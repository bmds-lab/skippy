/*
        MAIN.CPP
        --------
        indexReference

        Created by Shlomo Geva on 16/7/2023.
        Updated: runtime stripe pattern, metadata file output.
        Updated: -genome flag for directory-based genome management.
*/

#include <errno.h>
#include <sys/stat.h>
#include <unistd.h>

#include <chrono>
#include <climits>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <thread>
#include <vector>

#include "SkippyConstants.hpp"
#include "generate_seed_matrix.hpp"
#include "greedyOrder.hpp"
#include "indexGenome.hpp"
#include "protected_vector.hpp"
#include "serialiseKmersMap.hpp"
#include "stripe_pattern.h"

std::vector<std::string> sampleSequences;
std::vector<std::string> fileNames;

// Global run parameters
std::string REFERENCE = "";   // full path to FASTA (legacy) or derived from -genome
std::string GENOME_DIR = "";  // genome root directory (-genome flag)
std::string STRIPE_ARG = "";  // raw stripe string from command line
std::string STRIPE_ARG2 = ""; // optional complementary pattern
StripePattern gPattern;       // parsed at startup, used everywhere
StripePattern gPattern2;      // optional complementary pattern

// Technology → stripe mapping (for symlink maintenance)
static const std::vector<std::string> PLAIN_TECHS = {"ont", "nanopore", "pacbio"};
static const std::vector<std::string> STRIPED_TECHS = {"illumina", "hifi"};

/*
        WRITEMAPTOFILE()
        ----------------
*/
bool writeMapToFile(const std::string& filename,
                    const std::map<uint32_t, std::string>& referenceIDMap)
{
    std::ofstream outFile(filename);
    if (!outFile) {
        std::cerr << "Error opening the file: " << filename << std::endl;
        return false;
    }
    for (const auto& entry : referenceIDMap)
        outFile << entry.first << " " << entry.second;
    outFile.close();
    if (!outFile)
        std::cerr << "Error writing the file: " << filename << '\n';
    return static_cast<bool>(outFile);
}

/*
        RESOLVEREALPATH()
        -----------------
*/
std::string resolveRealPath(const std::string& path)
{
    char resolved[PATH_MAX];
    if (realpath(path.c_str(), resolved) != nullptr)
        return std::string(resolved);
    return path;
}

/*
        GETBASENAME()
        -------------
        Strips the file extension from a path.
*/
std::string getBaseName(const std::string& filePath)
{
    size_t lastSlash = filePath.find_last_of("/\\");
    size_t searchFrom = (lastSlash == std::string::npos) ? 0 : lastSlash + 1;
    size_t lastDot = filePath.find_last_of('.');
    if (lastDot != std::string::npos && lastDot > searchFrom)
        return filePath.substr(0, lastDot);
    return filePath;
}

/*
        DIREXISTS()
        -----------
*/
bool dirExists(const std::string& path)
{
    struct stat st;
    return (stat(path.c_str(), &st) == 0 && S_ISDIR(st.st_mode));
}

/*
        FILEEXISTS()
        ------------
*/
bool fileExists(const std::string& path)
{
    struct stat st;
    return (stat(path.c_str(), &st) == 0 && S_ISREG(st.st_mode));
}

/*
        UPDATESYMLINKS()
        ----------------
        Creates/updates technology symlinks in the genome directory.
        Symlink target is the stripe subdirectory name (relative, no leading
   slash).
*/
void updateSymlinks(const std::string& genomeDir, const std::string& stripeSubdir, bool isPlain)
{
    const auto& techs = isPlain ? PLAIN_TECHS : STRIPED_TECHS;
    for (const auto& tech : techs) {
        std::string linkPath = genomeDir + "/" + tech;
        // Remove existing symlink if present
        struct stat st;
        if (lstat(linkPath.c_str(), &st) == 0) {
            if (S_ISLNK(st.st_mode))
                unlink(linkPath.c_str());
            else {
                std::cerr << "Warning: " << linkPath << " exists and is not a symlink — skipping\n";
                continue;
            }
        }
        if (symlink(stripeSubdir.c_str(), linkPath.c_str()) != 0)
            std::cerr << "Warning: could not create symlink " << linkPath << " -> " << stripeSubdir
                      << ": " << strerror(errno) << "\n";
        else
            std::cout << "Symlink: " << tech << " -> " << stripeSubdir << "\n";
    }
}

/*
        WRITEMETADATA()
        ---------------
*/
bool writeMetadata(const std::string& inputFile, const StripePattern& pat,
                   const StripePattern& pat2, uint32_t genomeSize, int numBitsToKeep)
{
    std::string metaFilename = getBaseName(inputFile) + "_meta.idx";
    std::ofstream meta(metaFilename);
    if (!meta) {
        std::cerr << "Error: could not write metadata to " << metaFilename << "\n";
        return false;
    }
    meta << "SKIPPY_VERSION=2\n";
    meta << "STRIPE_PATTERN=" << pat.pattern_str << "\n";
    if (!pat2.pattern_str.empty())
        meta << "STRIPE_PATTERN2=" << pat2.pattern_str << "\n";
    meta << "KMER_WEIGHT=" << pat.weight << "\n";
    meta << "KMER_WINDOW=" << pat.window << "\n";
    meta << "GENOME_SIZE=" << genomeSize << "\n";
    meta << "KEFF=" << numBitsToKeep << "\n";
    meta.close();
    if (!meta) {
        std::cerr << "Error: could not write metadata to " << metaFilename << "\n";
        return false;
    }
    std::cout << "Metadata written to " << metaFilename << "\n";
    return true;
}

/*
        GETREFERENCE()
        --------------
        Build the index for REFERENCE (FASTA path) into outputDir.
        outputDir is the full path to the stripe subdirectory.
*/
bool buildIndex(const std::string& fastaFile, const std::string& outputDir,
                const StripePattern& pat, const StripePattern& pat2)
{
    std::map<uint32_t, std::string> referenceIDMap;
    char* genome = nullptr;
    uint32_t MASK = 0;

    // Load the genome
    auto start = std::chrono::steady_clock::now();
    uint32_t genomeSize;
    genome = load_genome_file(fastaFile, referenceIDMap, genomeSize);

    int numBitsToKeep = (int)::ceil(::log2(genomeSize));
    MASK = (numBitsToKeep == 32) ? UINT32_MAX : (1 << numBitsToKeep) - 1;
    std::cout << "Keeping " << numBitsToKeep << " bits in kmerHash\n";
    std::vector<protected_vector<uint32_t>> kmersMap(pow(2, numBitsToKeep));

    index_kmers(genome, genomeSize, kmersMap, MASK, pat, pat2);

    auto end = std::chrono::steady_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    int minutes = (int)duration.count() / (1000 * 60);
    float seconds = ((duration.count() - minutes * 1000 * (float)60)) / 1000;
    std::cout << "Building time: " << minutes << " min " << seconds << " sec\n";

    uint64_t kmerCount = genomeSize - pat.window;
    uint64_t kmersInMap = 0;
    for (size_t i = 0; i < kmersMap.size(); i++)
        if (!kmersMap[i].empty())
            kmersInMap++;
    std::cout << "Map size " << kmersMap.size() << ", kmersCount " << kmerCount << ", kmers in Map "
              << kmersInMap << "\n";

    // Output base path = outputDir/<genomeName>
    // Derive genome name from FASTA basename
    std::string genomeName = getBaseName(fastaFile);
    size_t slash = genomeName.find_last_of("/\\");
    if (slash != std::string::npos)
        genomeName = genomeName.substr(slash + 1);

    std::string basePath = outputDir + "/" + genomeName;

    start = std::chrono::steady_clock::now();
    std::string outerMapFilename = basePath + "_OuterBlob.idx";
    std::string innerMapFilename = basePath + "_InnerBlob.idx";
    std::string genomeFilename = basePath + "_genome.idx";
    std::string refIDFilename = basePath + "_refID.idx";

    std::cout << "Serialising genome to " << genomeFilename << "\n";
    if (!writeTextBlobToFile(genome, genomeSize, genomeFilename)) {
        std::cerr << "Error: could not write " << genomeFilename << '\n';
        return false;
    }

    std::cout << "Serialising map to " << outerMapFilename << " and " << innerMapFilename << "\n";
    if (!serializeMap(kmersMap, innerMapFilename, outerMapFilename))
        return false;

    end = std::chrono::steady_clock::now();
    duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    minutes = (int)duration.count() / (1000 * 60);
    seconds = ((duration.count() - minutes * 1000 * (float)60)) / 1000;
    std::cout << "Serialising time: " << minutes << " min " << seconds << " sec\n";

    std::cout << "Serialising ReferenceIDMap\n";
    if (!writeMapToFile(refIDFilename, referenceIDMap))
        return false;

    if (!writeMetadata(basePath, pat, pat2, genomeSize, numBitsToKeep))
        return false;

    std::string greedyFile = basePath + "_greedy_order.bin";
    std::cout << "Generating greedy order table...\n";
    initialiseGreedyOrders(pat);
    if (!saveGreedyOrders(greedyFile)) {
        std::cerr << "Error: failed to write " << greedyFile << std::endl;
        return false;
    }

    std::cout << "Generating seed matrices...\n";
    constexpr int MATRIX_TRIALS = 10000;
    constexpr double TARGET_PROB = 0.999;

    std::cout << "SUBSTITUTION / GREEDY seed matrix...\n";
    if (!generateSeedMatrix(pat, PlacementStrategy::GREEDY, ErrorType::SUBSTITUTION, basePath,
                            MATRIX_TRIALS, MAX_READ_LENGTH, MAX_ERROR_PCT, TARGET_PROB))
        return false;

    std::cout << "INDEL / GREEDY seed matrix...\n";
    return generateSeedMatrix(pat, PlacementStrategy::GREEDY, ErrorType::INDEL, basePath,
                              MATRIX_TRIALS, MAX_READ_LENGTH, MAX_ERROR_PCT, TARGET_PROB);
}

/*
        INITIALISE()
        ------------
*/
void initialise(int argc, char* argv[])
{
    if (argc <= 1) {
        std::cerr << "Error: specify -genome <dir> or -reference <fasta>; see skippy index --help\n";
        exit(2);
    }
    if (strcmp(argv[1], "-help") == 0 || strcmp(argv[1], "--help") == 0) {
        std::cout << "Usage:\n";
        std::cout << "  skippy index -genome <dir> [-stripe <pattern>]\n";
        std::cout << "  skippy index -reference <fasta> [-stripe <pattern>]  "
                     "(legacy)\n\n";
        std::cout << "  -genome <dir>    Genome root directory, e.g. "
                     "/data/genomes/GRCh38\n";
        std::cout << "                   FASTA: <dir>/<name>.fasta  (name = "
                     "basename of dir)\n";
        std::cout << "                   Index: <dir>/<stripe>/\n";
        std::cout << "  -stripe <pat>    Stripe pattern, e.g. 4XX4XXX1XX3X4 or "
                     "plain (16-mer)\n";
        std::cout << "                   Default: plain\n";
        std::cout << "  -stripe2 <pat>   Optional complementary stripe pattern\n";
        std::cout << "  -reference <f>   Legacy: full path to FASTA file\n\n";
        std::cout << "  --help           Show this help\n";
        std::cout << "  --version        Show version\n\n";
        std::cout << "Examples:\n";
        std::cout << "  skippy index -genome /data/genomes/GRCh38 -stripe plain\n";
        std::cout << "  skippy index -genome /data/genomes/GRCh38 -stripe "
                     "4XX4XXX1XX3X4\n";
        exit(0);
    }

    for (int i = 1; i < argc; i++) {
        std::string arg = argv[i];

        if (arg == "STRIPE") {
            if (STRIPE_ARG.empty())
                STRIPE_ARG = DEFAULT_STRIPE_PATTERN;
            continue;
        }

        if (arg != "-genome" && arg != "-reference" && arg != "-stripe" &&
            arg != "-stripe2") {
            std::cerr << "Error: unknown index option: " << arg << '\n';
            exit(2);
        }
        if (i + 1 >= argc || argv[i + 1][0] == '-') {
            std::cerr << "Error: missing value for " << arg << "; see skippy index --help\n";
            exit(2);
        }
        std::string value = argv[++i];

        if (arg == "-genome")
            GENOME_DIR = resolveRealPath(value);
        else if (arg == "-reference")
            REFERENCE = resolveRealPath(value);
        else if (arg == "-stripe")
            STRIPE_ARG = value;
        else if (arg == "-stripe2")
            STRIPE_ARG2 = value;
    }

    // Default stripe = plain
    if (STRIPE_ARG.empty())
        STRIPE_ARG = "plain";

    // Parse stripe pattern (plain = all-care 16-mer = 16 ones)
    bool isPlain = (STRIPE_ARG == "plain");
    if (!isPlain) {
        if (!gPattern.parse(STRIPE_ARG)) {
            std::cerr << "Fatal: invalid stripe pattern \"" << STRIPE_ARG << "\"\n";
            exit(2);
        }
    }
    else {
        // Plain 16-mer: pattern "1111111111111111" (16 care positions, no gaps)
        if (!gPattern.parse("16")) {
            std::cerr << "Fatal: could not parse plain 16-mer pattern\n";
            exit(1);
        }
    }

    if (!STRIPE_ARG2.empty()) {
        if (!gPattern2.parse(STRIPE_ARG2)) {
            std::cerr << "Fatal: invalid stripe pattern2 \"" << STRIPE_ARG2 << "\"\n";
            exit(2);
        }
        if (gPattern2.weight != gPattern.weight) {
            std::cerr << "Fatal: stripe2 weight (" << gPattern2.weight
                      << ") must match stripe weight (" << gPattern.weight << ")\n";
            exit(2);
        }
    }
}

/*
        MAIN()
        ------
*/
int main(int argc, char* argv[])
{
    auto startAll = std::chrono::steady_clock::now();

    initialise(argc, argv);

    bool isPlain = (STRIPE_ARG == "plain");

    if (!GENOME_DIR.empty()) {
        //--------------------------------------------------------------
        // New -genome mode
        //--------------------------------------------------------------

        // Derive genome name from directory basename
        std::string genomeName = GENOME_DIR;
        size_t slash = genomeName.find_last_of("/\\");
        if (slash != std::string::npos)
            genomeName = genomeName.substr(slash + 1);

        std::string fastaFile = GENOME_DIR + "/" + genomeName + ".fasta";

        // 1. Check FASTA exists
        if (!fileExists(fastaFile)) {
            std::cerr << "Error: Genome file not found: " << fastaFile << "\n";
            std::cerr << "Please place the genome FASTA at <genome_dir>/<name>.fasta\n";
            return 1;
        }

        // 2. Determine stripe subdirectory
        std::string stripeSubdir = STRIPE_ARG; // "plain" or "4XX4XXX1XX3X4" etc.
        std::string indexDir = GENOME_DIR + "/" + stripeSubdir;

        // 3. Error if index directory already exists
        if (dirExists(indexDir)) {
            std::cerr << "Error: Index directory already exists: " << indexDir << "\n";
            std::cerr << "To rebuild, delete it first:\n";
            std::cerr << "  rm -rf " << indexDir << "\n";
            return 1;
        }

        // 4. If building a striped index AND plain doesn't exist, build plain first
        if (!isPlain) {
            std::string plainDir = GENOME_DIR + "/plain";
            if (!dirExists(plainDir)) {
                std::cout << "\nPlain index not found — building plain index first...\n\n";
                if (mkdir(plainDir.c_str(), 0755) != 0) {
                    std::cerr << "Error: could not create " << plainDir << ": " << strerror(errno)
                              << "\n";
                    return 1;
                }
                StripePattern plainPat;
                plainPat.parse("16");
                if (!buildIndex(fastaFile, plainDir, plainPat, StripePattern{}))
                    return 1;
                updateSymlinks(GENOME_DIR, "plain", true);
                std::cout << "\nPlain index built. Now building " << STRIPE_ARG << " index...\n\n";
            }
        }

        // 5. Create stripe index directory
        if (mkdir(indexDir.c_str(), 0755) != 0) {
            std::cerr << "Error: could not create " << indexDir << ": " << strerror(errno) << "\n";
            return 1;
        }

        // 6. Build the requested index
        std::cout << "skippy index parameters\n";
        std::cout << "  Genome:  " << GENOME_DIR << "\n";
        std::cout << "  FASTA:   " << fastaFile << "\n";
        std::cout << "  Stripe:  " << STRIPE_ARG << "\n";
        std::cout << "  Index:   " << indexDir << "\n\n";
        gPattern.print();

        if (!buildIndex(fastaFile, indexDir, gPattern, gPattern2))
            return 1;

        // 7. Update symlinks
        updateSymlinks(GENOME_DIR, stripeSubdir, isPlain);
    }
    else if (!REFERENCE.empty()) {
        //--------------------------------------------------------------
        // Legacy -reference mode (unchanged behaviour)
        //--------------------------------------------------------------
        if (!fileExists(REFERENCE)) {
            std::cerr << "Error: reference FASTA not found: " << REFERENCE << '\n';
            return 1;
        }
        std::cout << "skippy index run parameters\n";
        std::cout << "reference: " << REFERENCE << "\n";
        if (!STRIPE_ARG.empty() && !isPlain)
            gPattern.print();
        else
            std::cout << "Stripe: plain 16-mer\n";

        // In legacy mode, output goes alongside the FASTA file
        std::string outputDir = REFERENCE.substr(0, REFERENCE.find_last_of("/\\"));
        if (!buildIndex(REFERENCE, outputDir, gPattern, gPattern2))
            return 1;
    }
    else {
        std::cerr << "Error: specify -genome <dir> or -reference <fasta>; see skippy index --help\n";
        return 2;
    }

    auto endAll = std::chrono::steady_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(endAll - startAll);
    int minutes = (int)duration.count() / (1000 * 60);
    float seconds = ((duration.count() - minutes * 1000 * (float)60)) / 1000;
    std::cout << "\nTotal time: " << minutes << " min " << seconds << " sec\n";

    return 0;
}
