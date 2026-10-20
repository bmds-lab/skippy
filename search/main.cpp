#include "main.hpp"

#include <dirent.h>
#include <fcntl.h>
#include <limits.h>
#include <limits.h> // For PATH_MAX
#include <stdio.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

#include <algorithm>
#include <atomic>
#include <cctype>
#include <chrono>
#include <cmath>
#include <condition_variable>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <functional>
#include <iomanip> // for std::setprecision
#include <iostream>
#include <map>
#include <memory> // for shared_ptr
#include <mutex>
#include <queue>
#include <set>
#include <span>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

#include "MonteCarlo.hpp"
#include "SafeQueue.hpp"
#include "SeedMatrix.hpp"
#include "SkippyConstants.hpp"
#include "encode_kmer_2bit.h"
#include "hash.hpp"
#include "kmerUtilities.hpp"
#include "seed_table.h"
#include "serialiseKmersMap.hpp"
#include "smithWaterman.hpp"
#include "ssw.h"
#include "stripe_pattern.hpp"
#include "threadPool.hpp"
using std::cerr;
using std::cout;
using std::endl;

// SeedTable gSeedTable;
SeedMatrix gSeedMatrix;

uint32_t ERROR_TARGET_PCT;
uint32_t GATE_COVERAGE;

// Runtime stripe pattern — loaded from index metadata at startup
StripePattern gPattern;
bool TECHNOLOGY_SPECIFIED = false;
std::string TECHNOLOGY_NAME = "";
bool MAX_SEEDS_SPECIFIED = false;
bool ERROR_TARGET_SPECIFIED = false;
std::string CONTRACT_TYPE = "subs"; // "subs" or "indels"
bool CONTRACT_SPECIFIED = false;
int SEED_STEP = 0;
bool SEED_STEP_SPECIFIED = false;
int NUM_SEEDS = 0;
bool NUM_SEEDS_SPECIFIED = false;
bool PAIRED_READS = false;                // if false, disable mate rescue
bool SAM_OUTPUT = false;
bool DEBUG_COVERAGE = false;              // if true, dump per-read seeds/coverage to stderr
bool BALANCED_MODE = false;               // if true, apply length-dependent secondary gate in getSW
float IDENTITY_THRESHOLD = 0.63f;         // SW identity filter threshold (eq/read_len)
uint32_t SUSPICIOUS_SPAN_THRESHOLD = 200; // span below this triggers SW even with cutoff=0
bool useGreedy = true;                    // greedy placement of seeds based on stripe pattern
std::string REFERENCE;                    // Reference file (genome) — full path to FASTA
std::string GENOME_DIR;                   // Genome root directory (-genome flag)
std::string STRIPE_NAME;                  // Stripe subdirectory name (-stripe flag)
std::vector<std::vector<uint32_t>> gGreedyOrders;
static constexpr uint32_t MAX_LIST_LEN = 250000; // tune

// =============================================================================
// SCORING VERSION CONTROL
// Must match the setting in smithWaterman.cpp!
// Uncomment ONE:
// =============================================================================
// #define USE_ABSOLUTE_MATCH_COUNT    // Original: count of '=' operations
#define USE_ALIGNMENT_DENSITY // New: density-based scoring
#ifdef USE_ALIGNMENT_DENSITY
// New density-based scoring thresholds
uint32_t MIN_ALIGNMENT_SPAN = 80;    // Minimum aligned region length (bp)
float MIN_ALIGNMENT_DENSITY = 0.70f; // Minimum match density (0.0-1.0)
uint32_t MIN_ABSOLUTE_MATCHES = 33;  // Minimum exact matches (bp)
#endif

// Logging globals
std::mutex logMutex;
std::condition_variable logCV;
std::map<size_t, std::vector<std::string>> fileLogMessages;
std::atomic<bool> loggingFinished{false};
size_t nextFileToPrint = 0;

void initialiseGreedyOrders(const StripePattern& pattern);
std::vector<int> buildGreedyOrder(int readLength, int span, const std::vector<int>& careOffsets,
                                  const std::vector<int>& protectOffsets);

// Holds exactly one chunk worth of ReadInfo plus shared ownership of the mmap'd
// file bytes.
struct ChunkData {
    std::vector<ReadInfo> readInfos;  // only for this chunk
    std::shared_ptr<char[]> data_ptr; // mmap backing store (shared across chunks)
};

struct FileContext {
    SafeQueue<Result> writerQueue;
    std::thread writerThread;
    std::vector<ReadInfo> samReads;
    std::shared_ptr<char[]> samData;
    bool writeFailed = false;

    size_t fileNum = 0;
    std::string outFilename;

    // Streaming accounting (we do NOT know totalReads upfront anymore)
    std::atomic<size_t> totalReads{0};

    // Completion signalling for this file
    std::atomic<uint64_t> inFlight{0};    // #chunks enqueued but not yet finished by workers
    std::atomic<bool> doneEnqueue{false}; // producer finished enqueueing all chunks for this file
};

struct WorkItem {
    std::shared_ptr<FileContext> ctx;
    std::shared_ptr<ChunkData> chunkData; // holds chunk readInfos + mmap ptr
    std::span<const ReadInfo> chunk;      // view into chunkData->readInfos
    size_t startIndex;                    // global read index within this file (for getSW)
};

std::vector<std::unique_ptr<FileContext>> contexts; // indexed by fileNum
std::mutex contextsMutex;                           // if you access contexts concurrently

/*
void buildGreedySeedTable();
*/

// Logging helper
// Appends " [LAST]" iff isLast==true, buffers per fileNum, and notifies logger.
void logForFile(size_t fileNum, const std::string& message, bool isLast)
{
    {
        std::lock_guard<std::mutex> lock(logMutex);
        std::string msg = message;
        if (isLast && msg.find("[LAST]") == std::string::npos)
            msg += " [LAST]";
        fileLogMessages[fileNum].push_back(std::move(msg));
    }
    logCV.notify_all();
}

// Drop-in overload: preserve existing call sites that pass only (fileNum,
// message)
inline void logForFile(size_t fileNum, const std::string& message)
{
    logForFile(fileNum, message, false);
}

std::string getSeedTableFilename()
{
    namespace fs = std::filesystem;

    fs::path p(REFERENCE);

    return (p.parent_path() / "seed_counts.bin").string();
}

std::vector<std::tuple<uint32_t, uint32_t, uint32_t, uint32_t>> results;
std::mutex resultsMutex;
std::vector<std::string> fileNames;
uint32_t query_length;
std::vector<std::tuple<uint32_t, uint32_t, uint32_t, double, std::string, uint32_t>> resultsSW;

// uint32_t KMERSIZE;
uint32_t MIN_MATCHES;
uint32_t MAX_SEEDS = 1000;
uint32_t THREADS;
int CUTOFF;
// std::string REFERENCE;
std::string READS;
std::string OUTPUT_DIR;
uint64_t kmer_encoding_table[256];
bool MATCH_ALL;
bool GET_SW;
uint32_t INDEL_SLACK;
QueryTable ref_ID_Table;

// ConcurrentQueue template (BOUNDED: blocks enqueue when full)
template <typename T> class ConcurrentQueue {
public:
    ConcurrentQueue() = default;

    // Set capacity after you know THREADS. 0 => unbounded (disabled).
    void setCapacity(size_t cap)
    {
        std::lock_guard<std::mutex> lock(mutex_);
        capacity_ = cap;
        // If capacity increased, wake producer(s)
        cv_not_full_.notify_all();
    }

    void enqueue(T item)
    {
        std::unique_lock<std::mutex> lock(mutex_);

        // If bounded, block when full.
        // (If capacity_ == 0 => unbounded.)
        if (capacity_ > 0) {
            cv_not_full_.wait(lock, [this] { return finished_ || queue_.size() < capacity_; });
        }

        // If finished, ignore further enqueues (producer should be done anyway)
        if (finished_)
            return;

        queue_.push(std::move(item));
        cv_not_empty_.notify_one();
    }

    bool dequeue(T& item)
    {
        std::unique_lock<std::mutex> lock(mutex_);
        cv_not_empty_.wait(lock, [this] { return finished_ || !queue_.empty(); });

        if (queue_.empty())
            return false;

        item = std::move(queue_.front());
        queue_.pop();

        // Wake producer if it was blocked on "full"
        if (capacity_ > 0)
            cv_not_full_.notify_one();

        return true;
    }

    void setFinished()
    {
        std::lock_guard<std::mutex> lock(mutex_);
        finished_ = true;
        cv_not_empty_.notify_all();
        cv_not_full_.notify_all();
    }

    bool isFinished() const
    {
        std::lock_guard<std::mutex> lock(mutex_);
        return finished_ && queue_.empty();
    }

private:
    std::queue<T> queue_;
    mutable std::mutex mutex_;
    std::condition_variable cv_not_empty_;
    std::condition_variable cv_not_full_;
    bool finished_ = false;

    size_t capacity_ = 0; // 0 => unbounded
};
ConcurrentQueue<WorkItem> workQueue;

std::string get_current_dir_name_custom()
{
    char cwd[PATH_MAX];
    if (getcwd(cwd, sizeof(cwd)) != nullptr) {
        return std::string(cwd);
    }
    else {
        perror("getcwd() error");
        return std::string();
    }
}
std::string getBaseFilename(const std::string& filePath)
{
    std::size_t lastSlash = filePath.find_last_of("/\\");
    std::size_t lastDot = filePath.find_last_of('.');
    if (lastDot != std::string::npos && (lastSlash == std::string::npos || lastDot > lastSlash)) {
        return filePath.substr(lastSlash + 1, lastDot - lastSlash - 1);
    }
    return filePath.substr(lastSlash + 1);
}

std::string getFirstToken(const std::string& str)
{
    std::istringstream iss(str);
    std::string firstToken;
    iss >> firstToken;
    return firstToken;
}

bool isRegularFile(const std::string& path)
{
    struct stat pathStat;
    stat(path.c_str(), &pathStat);
    return S_ISREG(pathStat.st_mode);
}

bool hasReadExtension(const std::string& fileName)
{
    std::vector<std::string> readExtensions = {".fasta", ".fa", ".fna", ".fastq", ".fq"};
    std::string::size_type idx = fileName.find_last_of('.');
    if (idx != std::string::npos) {
        std::string extension = fileName.substr(idx);
        return std::find(readExtensions.begin(), readExtensions.end(), extension) !=
               readExtensions.end();
    }
    return false;
}

std::vector<std::string> getFileNames(const std::string& rootDir)
{
    std::vector<std::string> fileNames;

    if (isRegularFile(rootDir) && hasReadExtension(rootDir)) {
        fileNames.push_back(rootDir);
        return fileNames;
    }

    DIR* directory = opendir(rootDir.c_str());
    if (directory == nullptr) {
        std::cerr << "Error: Could not open directory " << rootDir << std::endl;
        return fileNames;
    }

    dirent* entry;
    while ((entry = readdir(directory)) != nullptr) {
        std::string fileName = entry->d_name;
        std::string filePath = rootDir + "/" + fileName;

        if (fileName == "." || fileName == "..") {
            continue;
        }

        if (isRegularFile(filePath) && hasReadExtension(filePath)) {
            fileNames.push_back(filePath);
        }
        else if (entry->d_type == DT_DIR) {
            std::vector<std::string> subDirFiles = getFileNames(filePath);
            fileNames.insert(fileNames.end(), subDirFiles.begin(), subDirFiles.end());
        }
    }

    closedir(directory);
    return fileNames;
}

std::string extractRefName(const std::string& filePath)
{
    std::string fileName = filePath.substr(filePath.find_last_of('/') + 1);
    std::string baseName = fileName.substr(0, fileName.find_last_of('.'));
    return baseName;
}

std::string getBaseName(const std::string& filePath)
{
    std::stringstream ss(filePath);
    std::string baseName;
    std::getline(ss, baseName, '.');
    return baseName;
}

static bool samCigarSpans(const std::string& cigar, size_t& querySpan, size_t& referenceSpan)
{
    querySpan = referenceSpan = 0;
    size_t count = 0;
    for (char op : cigar) {
        if (op >= '0' && op <= '9') {
            count = count * 10 + static_cast<size_t>(op - '0');
            continue;
        }
        if (count == 0)
            return false;
        switch (op) {
        case 'M':
        case '=':
        case 'X':
            querySpan += count;
            referenceSpan += count;
            break;
        case 'I':
        case 'S':
            querySpan += count;
            break;
        case 'D':
        case 'N':
            referenceSpan += count;
            break;
        case 'H':
        case 'P':
            break;
        default:
            return false;
        }
        count = 0;
    }
    return count == 0 && referenceSpan > 0;
}

static std::string samReferenceName(const std::string& header)
{
    std::string name = getFirstToken(header);
    if (!name.empty() && name.front() == '>')
        name.erase(0, 1);
    return name;
}

void writer(FileContext& ctx, size_t genomeLength)
{
    auto& queue = ctx.writerQueue;
    const auto& filename = ctx.outFilename;
    auto& totalReads = ctx.totalReads;
    const size_t fileNum = ctx.fileNum;
    std::map<uint32_t, Result> resultsMap;
    Result result;

    while (true) {
        bool dequeued = queue.dequeue(result);

        if (dequeued) {
            resultsMap[result.sampleIndex] = result;
        }

        else if (queue.isFinished()) {
            break;
        }
    }

    std::set<uint32_t> missingMates;
    if (PAIRED_READS) {
        for (const auto& [index, res] : resultsMap) {
            uint32_t realIndex = res.sampleIndex + 1;

            if (realIndex % 2 == 1) {
                if (resultsMap.find(index + 1) == resultsMap.end()) {
                    missingMates.insert(index + 1);
                }
            }
            else {
                if (resultsMap.find(index - 1) == resultsMap.end()) {
                    missingMates.insert(index - 1);
                }
            }
        }

        logForFile(fileNum, "Inserted mates  : " + std::to_string(missingMates.size()));

        for (uint32_t idx : missingMates) {
            Result dummy;
            dummy.sampleIndex = idx;
            dummy.refBegin = 0;
            dummy.score = 0.0f;
            dummy.cigar = "*";
            dummy.revComp = 0;
            resultsMap[idx] = dummy;
        }
    }
    std::ofstream outFile(filename);
    if (!outFile.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        ctx.writeFailed = true;
        logForFile(fileNum, "Error opening file: " + filename, true);
        return;
    }

    char buffer[1024];
    uint32_t writtenCount = 0;

    if (SAM_OUTPUT) {
        outFile << "@HD\tVN:1.6\tSO:unsorted\n";
        for (auto it = ref_ID_Table.entries.begin(); it != ref_ID_Table.entries.end(); ++it) {
            auto next = std::next(it);
            size_t end = next == ref_ID_Table.entries.end() ? genomeLength : next->first;
            outFile << "@SQ\tSN:" << samReferenceName(it->second) << "\tLN:" << end - it->first
                    << '\n';
        }

        for (size_t index = 0; index < ctx.samReads.size(); ++index) {
            const ReadInfo& ri = ctx.samReads[index];
            const char* data = ctx.samData.get();
            std::string name(data + ri.nameStart, ri.nameLength);
            if (name.empty())
                name = std::to_string(index + 1);
            std::string sequence(data + ri.start, ri.length);
            std::string quality = ri.qualityLength == ri.length
                                      ? std::string(data + ri.qualityStart, ri.qualityLength)
                                      : "*";

            auto result = resultsMap.find(static_cast<uint32_t>(index));
            bool mapped = result != resultsMap.end() && !result->second.cigar.empty() &&
                          result->second.cigar != "*" && result->second.refBegin < genomeLength;
            auto reference = ref_ID_Table.entries.end();
            if (mapped) {
                reference = ref_ID_Table.entries.upper_bound(result->second.refBegin);
                if (reference == ref_ID_Table.entries.begin()) {
                    mapped = false;
                }
                else {
                    --reference;
                    auto next = std::next(reference);
                    const size_t end = next == ref_ID_Table.entries.end() ? genomeLength
                                                                           : next->first;
                    size_t querySpan = 0, referenceSpan = 0;
                    mapped = samCigarSpans(result->second.cigar, querySpan, referenceSpan) &&
                             querySpan == ri.length && referenceSpan <= end - result->second.refBegin;
                }
            }

            if (mapped) {
                const Result& res = result->second;
                if (res.revComp) {
                    sequence = reverseComplement(sequence);
                    if (quality != "*")
                        std::reverse(quality.begin(), quality.end());
                }
                outFile << name << '\t' << (res.revComp ? 16 : 0) << '\t'
                        << samReferenceName(reference->second) << '\t'
                        << (res.refBegin - reference->first + 1) << "\t255\t" << res.cigar
                        << "\t*\t0\t0\t" << sequence << '\t' << quality << '\n';
                ++writtenCount;
            }
            else {
                outFile << name << "\t4\t*\t0\t0\t*\t*\t0\t0\t" << sequence << '\t'
                        << quality << '\n';
            }
        }
    }
    else {
        for (const auto& [index, res] : resultsMap) {
            const bool missingMate = missingMates.count(index) != 0;
            std::string refID = missingMate
                                    ? "*"
                                    : getFirstToken(ref_ID_Table.getValue(ref_ID_Table.query(res.refBegin)));

            char statusChar = missingMates.find(index) != missingMates.end()
                                  ? '='
                                  : (res.revComp == 0 ? '>' : '<');

            snprintf(buffer, sizeof(buffer), "%c %8u\t%s\t%8u\t%5.2f\t%s\n", statusChar,
                     res.sampleIndex + 1, refID.c_str(),
                     missingMate ? 0 : res.refBegin - ref_ID_Table.query(res.refBegin) + 1,
                     res.score, res.cigar.c_str());

            outFile << buffer;
            ++writtenCount;
        }
    }
    outFile.close();
    if (!outFile) {
        std::cerr << "Error writing file: " << filename << std::endl;
        ctx.writeFailed = true;
        logForFile(fileNum, "Error writing file: " + filename);
    }

    const size_t tr = totalReads.load(std::memory_order_relaxed);
    logForFile(fileNum,
               std::to_string(tr) + " reads, " + std::to_string(writtenCount) + " aligned");
    float alignRate = tr > 0 ? (float)(writtenCount) / tr * 100.0f : 0.0f;

    std::ostringstream rateStream;
    rateStream << std::fixed << std::setprecision(2) << alignRate << "% alignment rate";
    logForFile(fileNum, rateStream.str(), true); // last message
}

bool loadGreedyOrders(const std::string& filename)
{
    std::ifstream in(filename, std::ios::binary);

    if (!in)
        return false;

    uint32_t version;
    uint32_t maxLength;

    in.read(reinterpret_cast<char*>(&version), sizeof(version));

    in.read(reinterpret_cast<char*>(&maxLength), sizeof(maxLength));

    if (version != 1)
        return false;

    gGreedyOrders.clear();
    gGreedyOrders.resize(maxLength);

    for (uint32_t i = 0; i < maxLength; ++i) {
        uint32_t count;

        in.read(reinterpret_cast<char*>(&count), sizeof(count));

        gGreedyOrders[i].resize(count);

        if (count) {
            in.read(reinterpret_cast<char*>(gGreedyOrders[i].data()), count * sizeof(uint32_t));
        }
    }

    return true;
}

bool saveGreedyOrders(const std::string& filename)
{
    std::ofstream out(filename, std::ios::binary);

    if (!out)
        return false;

    uint32_t version = 1;
    uint32_t maxLength = static_cast<uint32_t>(gGreedyOrders.size());

    out.write(reinterpret_cast<char*>(&version), sizeof(version));

    out.write(reinterpret_cast<char*>(&maxLength), sizeof(maxLength));

    for (const auto& vec : gGreedyOrders) {
        uint32_t count = static_cast<uint32_t>(vec.size());

        out.write(reinterpret_cast<char*>(&count), sizeof(count));

        if (count) {
            out.write(reinterpret_cast<const char*>(vec.data()), count * sizeof(uint32_t));
        }
    }

    return true;
}

void getReference(const std::string& inputFile, std::span<uint32_t>& innerMapBlob,
                  std::span<uint32_t>& outerMapBlob, std::span<char>& genomeStr, uint32_t& MASK,
                  QueryTable& ref_ID_Table)
{
    std::string outerMapFilename = getBaseName(inputFile) + "_OuterBlob.idx";
    std::string innerMapFilename = getBaseName(inputFile) + "_InnerBlob.idx";
    std::string genomeFilename = getBaseName(inputFile) + "_genome.idx";
    std::string refIDFilename = getBaseName(inputFile) + "_refID.idx";
    std::string metaFilename = getBaseName(inputFile) + "_meta.idx";

    for (const auto& path : {innerMapFilename, outerMapFilename, genomeFilename, refIDFilename}) {
        struct stat st;
        if (stat(path.c_str(), &st) != 0 || !S_ISREG(st.st_mode)) {
            cerr << "Error: index file not found: " << path << '\n';
            exit(1);
        }
    }

    // ------------------------------------------------------------
    // Load stripe pattern from metadata file
    // ------------------------------------------------------------
    std::ifstream meta(metaFilename);
    if (!meta) {
        std::cerr << "Warning: no metadata file found at " << metaFilename
                  << " — using default pattern " << DEFAULT_STRIPE_PATTERN << "\n";
        gPattern.parse(DEFAULT_STRIPE_PATTERN);
    }
    else {
        std::string line, patternStr, pattern2Str;
        while (std::getline(meta, line)) {
            if (line.rfind("STRIPE_PATTERN2=", 0) == 0)
                pattern2Str = line.substr(16);
            else if (line.rfind("STRIPE_PATTERN=", 0) == 0)
                patternStr = line.substr(15);
        }
        meta.close();

        if (patternStr.empty()) {
            std::cout << "No stripe pattern in metadata — using default " << DEFAULT_STRIPE_PATTERN
                      << "\n";
            gPattern.parse(DEFAULT_STRIPE_PATTERN);
        }
        else {
            if (!gPattern.parse(patternStr)) {
                std::cerr << "FATAL: could not parse stripe pattern \"" << patternStr
                          << "\" from metadata\n";
                exit(1);
            }
        }
        gPattern.print();
    }

    ref_ID_Table = readTableFromFile(refIDFilename);

    auto start = std::chrono::steady_clock::now();
    genomeStr = readTextBlobFromFile(genomeFilename);
    if (genomeStr.empty()) {
        std::cerr << "Error: empty genome index: " << genomeFilename << '\n';
        exit(1);
    }
    logForFile(0, "genome length: " + std::to_string(genomeStr.size()));
    deserializeMap(innerMapFilename, outerMapFilename, innerMapBlob, outerMapBlob);
    logForFile(0, "outerBlob size: " + std::to_string(outerMapBlob.size()));
    logForFile(0, "innerBlob size: " + std::to_string(innerMapBlob.size()));

    auto end = std::chrono::steady_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    int minutes = (int)duration.count() / (1000 * 60);
    float seconds = ((duration.count() - minutes * 1000 * (float)60)) / 1000;
    logForFile(0, "DeSerialising time: " + std::to_string(minutes) + " min " +
                      std::to_string(seconds) + " sec");

    int numBitsToKeep = std::ceil(std::log2(genomeStr.size()));
    if (numBitsToKeep >= 32) {
        MASK = UINT32_MAX;
    }
    else {
        MASK = (uint32_t{1} << numBitsToKeep) - 1;
    }

    // After loading index and printing parameters:

    std::cout << "REFERENCE = " << REFERENCE << std::endl;
    std::string greedyFile = getBaseName(inputFile) + "_greedy_order.bin";

    if (!loadGreedyOrders(greedyFile)) {
        std::cerr << "ERROR: unable to load " << greedyFile << std::endl;

        exit(1);
    }

    //    gSeedMatrix.printSummary();
    //    std::cerr
    //    << "After gSeedMatrix.printSummary()"
    //    << std::endl;
}

// mmap_entire_file: maps file read-only; returns shared_ptr<char[]> with custom
// deleter
static std::shared_ptr<char[]> mmap_entire_file(const std::string& path, uint64_t& fileSize)
{
    fileSize = 0;

    int fd = ::open(path.c_str(), O_RDONLY);
    if (fd < 0) {
        perror(("open() failed: " + path).c_str());
        return {};
    }

    struct stat st;
    if (::fstat(fd, &st) != 0) {
        perror(("fstat() failed: " + path).c_str());
        ::close(fd);
        return {};
    }

    if (st.st_size <= 0) {
        ::close(fd);
        return {};
    }

    fileSize = static_cast<uint64_t>(st.st_size);

    void* p = ::mmap(nullptr, static_cast<size_t>(fileSize), PROT_READ, MAP_PRIVATE, fd, 0);
    if (p == MAP_FAILED) {
        perror(("mmap() failed: " + path).c_str());
        ::close(fd);
        fileSize = 0;
        return {};
    }

    // shared_ptr takes ownership; unmap+close on last release
    char* data = static_cast<char*>(p);

    return std::shared_ptr<char[]>(data, [fileSize, fd](char* ptr) {
        if (ptr && fileSize > 0)
            ::munmap(ptr, static_cast<size_t>(fileSize));
        if (fd >= 0)
            ::close(fd);
    });
}

// enqueue_fasta_or_fastq_streaming_chunks:
//   - FASTA: expects records starting with '>' followed by a single sequence
//   line.
//   - FASTQ: expects records starting with '@' followed by a single sequence
//   line,
//            then a '+' line and a single quality line.
// SAM output also retains header and quality offsets in the mapped file.
static void enqueue_fasta_or_fastq_streaming_chunks(const std::string& readsFile,
                                                    const std::shared_ptr<FileContext>& ctx,
                                                    uint64_t& fileSize, size_t fileNum,
                                                    ConcurrentQueue<WorkItem>& workQueue)
{
    // mmap the file
    std::shared_ptr<char[]> data_ptr = mmap_entire_file(readsFile, fileSize);
    if (!data_ptr || fileSize == 0) {
        std::cerr << "Empty or failed file: " << readsFile << std::endl;
        return;
    }

    const char* data = data_ptr.get();
    if (SAM_OUTPUT)
        ctx->samData = data_ptr;

    logForFile(fileNum, "Streaming Samples: " + readsFile, false);

    constexpr size_t READS_PER_CHUNK = 500;

    std::vector<ReadInfo> batch;
    batch.reserve(READS_PER_CHUNK);

    size_t pos = 0;
    size_t globalStartIndex = 0;

    auto flush_batch = [&]() {
        if (batch.empty())
            return;

        auto cd = std::make_shared<ChunkData>();
        cd->data_ptr = data_ptr;
        cd->readInfos = std::move(batch);

        WorkItem item;
        item.ctx = ctx;
        item.chunkData = cd;
        item.startIndex = globalStartIndex;
        item.chunk = std::span<const ReadInfo>(cd->readInfos.data(), cd->readInfos.size());

        ctx->totalReads.fetch_add(cd->readInfos.size(), std::memory_order_relaxed);
        ctx->inFlight.fetch_add(1, std::memory_order_acq_rel);

        workQueue.enqueue(std::move(item));

        globalStartIndex += cd->readInfos.size();

        batch.clear();
        batch.reserve(READS_PER_CHUNK);
    };

    auto is_line_start = [&](size_t p) -> bool {
        return (p == 0) || (data[p - 1] == '\n') || (data[p - 1] == '\r');
    };

    auto consume_eol = [&]() {
        if (pos < fileSize && data[pos] == '\r')
            ++pos;
        if (pos < fileSize && data[pos] == '\n')
            ++pos;
    };

    auto skip_to_eol = [&]() {
        while (pos < fileSize && data[pos] != '\n' && data[pos] != '\r')
            ++pos;
        consume_eol();
    };

    while (pos < fileSize) {
        // Find next record start at beginning of a line.
        // (Important for FASTQ: '@' appears in quality lines, so we require
        // line-start.)
        while (pos < fileSize) {
            if (is_line_start(pos) && (data[pos] == '>' || data[pos] == '@'))
                break;
            ++pos;
        }

        if (pos >= fileSize)
            break;

        const char recType = data[pos];

        // Skip header line
        const size_t nameStart = pos + 1;
        size_t nameEnd = nameStart;
        while (nameEnd < fileSize && !std::isspace(static_cast<unsigned char>(data[nameEnd])))
            ++nameEnd;
        skip_to_eol();

        if (pos >= fileSize)
            break;

        // Sequence line start
        size_t seqStart = pos;
        while (pos < fileSize && data[pos] != '\n' && data[pos] != '\r')
            ++pos;

        size_t seqLen = pos - seqStart;
        consume_eol();

        // FASTQ: skip '+' line and quality line (single-line versions).
        size_t qualityStart = 0, qualityLength = 0;
        if (recType == '@') {
            // Skip '+' line
            if (pos < fileSize)
                skip_to_eol();

            // Skip quality line
            if (pos < fileSize) {
                qualityStart = pos;
                while (pos < fileSize && data[pos] != '\n' && data[pos] != '\r')
                    ++pos;
                qualityLength = pos - qualityStart;
                consume_eol();
            }
        }

        if (seqLen > 0) {
            ReadInfo info{seqStart, seqLen, nameStart, nameEnd - nameStart, qualityStart,
                          qualityLength};
            if (SAM_OUTPUT)
                ctx->samReads.push_back(info);
            batch.push_back(info);
            if (batch.size() >= READS_PER_CHUNK)
                flush_batch();
        }
    }

    flush_batch();
}

void initialiseGreedyOrders(const StripePattern& pattern)
{
    gGreedyOrders.resize(MAX_READ_LENGTH + 1);

    for (int L = 1; L <= 1000; ++L) {
        std::vector<int> order =
            buildGreedyOrder(L, pattern.window, pattern.positions, pattern.protect_positions);

        gGreedyOrders[L] = std::vector<uint32_t>(order.begin(), order.end());
    }
}

bool getInvertedLists(std::vector<ListStruct>& result, const std::span<uint32_t>& outerMapBlob,
                      uint32_t index)
{
    if (index >= outerMapBlob.size())
        return false;

    uint32_t startIndex = outerMapBlob[index];

    uint32_t endIndex = (index + 1 < outerMapBlob.size())
                            ? outerMapBlob[index + 1]
                            : static_cast<uint32_t>(outerMapBlob.size());

    if (endIndex <= startIndex) {
        result.push_back(ListStruct{startIndex, 0});
        return false;
    }

    const uint32_t len = endIndex - startIndex;

    // Gate: treat "too long" lists as not-found
    if (len > MAX_LIST_LEN) {
        result.push_back(ListStruct{startIndex, 0});
        return false;
    }

    result.push_back(ListStruct{startIndex, len});
    return true;
}

bool getInvertedListsPtr(std::vector<ListStructPtr>& result,
                         const std::span<uint32_t>& outerMapBlob,
                         const std::span<uint32_t>& innerMapBlob, uint32_t index)
{
    if (index >= outerMapBlob.size())
        return false;

    uint32_t startIndex = outerMapBlob[index];

    uint32_t endIndex = (index + 1 < outerMapBlob.size())
                            ? outerMapBlob[index + 1]
                            : static_cast<uint32_t>(innerMapBlob.size());

    if (endIndex <= startIndex) {
        result.push_back(ListStructPtr{innerMapBlob.data() + startIndex, 0});
        return false;
    }

    const uint32_t len = endIndex - startIndex;

    // Gate: treat "too long" lists as not-found
    if (len > MAX_LIST_LEN) {
        result.push_back(ListStructPtr{innerMapBlob.data() + startIndex, 0});
        return false;
    }

    result.push_back(ListStructPtr{innerMapBlob.data() + startIndex, len});
    return true;
}

void testInvertedList(const std::string& seed, const std::span<uint32_t>& outerMapBlob,
                      const std::span<uint32_t>& innerMapBlob, const std::span<char>& genomeStr,
                      uint32_t index)
{
    if (index < outerMapBlob.size()) {
        std::cout << "Seed: " << seed << " ";

        uint32_t startIndex = outerMapBlob[index];
        uint32_t endIndex = (index + 1 < outerMapBlob.size()) ? outerMapBlob[index + 1]
                                                              : (uint32_t)outerMapBlob.size();
        if (endIndex > startIndex) {
            constexpr size_t pre_length = 0;
            constexpr size_t post_length = 31;
            std::cout << std::endl;

            for (uint32_t i = startIndex; i <= endIndex; ++i) {
                if (i >= innerMapBlob.size()) {
                    std::cerr << "Index out of range in innerMapBlob: " << i << std::endl;
                    continue;
                }

                uint32_t pos = innerMapBlob[i];
                if (pos >= genomeStr.size()) {
                    std::cerr << "Index out of range in genomeStr: " << pos << std::endl;
                    continue;
                }

                size_t pre_start = (pos < pre_length) ? 0 : pos - pre_length;
                size_t post_end = (pos + post_length >= genomeStr.size()) ? genomeStr.size() - 1
                                                                          : pos + post_length;

                std::string refStr(genomeStr.begin() + pre_start, genomeStr.begin() + post_end + 1);
                if (seed == toUpperCase(refStr))
                    std::cout << "*ref: " << toUpperCase(refStr) << " Position " << pos
                              << std::endl;
                else
                    std::cout << " ref: " << toUpperCase(refStr) << " Position " << pos
                              << std::endl;
            }
        }
        else
            std::cout << "Empty list\n ";
    }
}

/*============= START GREEDY SEED PLACEMENT CODE BLOCK =============*/
struct GreedyState {
    int readLength;
    int maxStart;

    std::vector<int> coverCount;
    std::vector<int> protectCount;

    std::vector<bool> protectable;

    std::set<int> remaining;

    std::vector<int> ordering;
};

static GreedyState initialiseState(int readLength, int span, const std::vector<int>& protectOffsets)
{
    GreedyState state;

    state.readLength = readLength;
    state.maxStart = readLength - span;

    state.coverCount.assign(readLength, 0);
    state.protectCount.assign(readLength, 0);
    state.protectable.assign(readLength, false);

    for (int start = 0; start <= state.maxStart; start++) {
        for (int off : protectOffsets) {
            int pos = start + off;

            if (pos >= 0 && pos < readLength) {
                state.protectable[pos] = true;
            }
        }
    }

    for (int s = 0; s <= state.maxStart; s++) {
        state.remaining.insert(s);
    }

    return state;
}

static bool coverageIncomplete(const GreedyState& state)
{
    for (int x : state.coverCount) {
        if (x == 0) {
            return true;
        }
    }

    return false;
}

static int computeCoverageScore(const GreedyState& state, int start,
                                const std::vector<int>& careOffsets)
{
    int score = 0;

    for (int off : careOffsets) {
        int pos = start + off;

        if (state.coverCount[pos] == 0) {
            score++;
        }
    }

    return score;
}

static double computeProtectionScore(const GreedyState& state, int start,
                                     const std::vector<int>& protectOffsets)
{
    double score = 0.0;

    for (int off : protectOffsets) {
        int pos = start + off;

        score += log((double)state.protectCount[pos] + 1.0);
    }

    return score;
}

static int chooseBestCoverageSeed(const GreedyState& state, const std::vector<int>& careOffsets)
{
    int bestStart = -1;
    int bestScore = -1;

    for (int start : state.remaining) {
        int score = computeCoverageScore(state, start, careOffsets);

        if (score > bestScore) {
            bestScore = score;
            bestStart = start;
        }
    }

    return bestStart;
}

static int chooseBestProtectionSeed(const GreedyState& state,
                                    const std::vector<int>& protectOffsets)
{
    int bestStart = -1;

    double bestScore = std::numeric_limits<double>::max();

    for (int start : state.remaining) {
        double score = computeProtectionScore(state, start, protectOffsets);

        if (score < bestScore) {
            bestScore = score;
            bestStart = start;
        }
    }

    return bestStart;
}

static void applySeed(GreedyState& state, int start, const std::vector<int>& careOffsets,
                      const std::vector<int>& protectOffsets)
{
    state.ordering.push_back(start);

    for (int off : careOffsets) {
        state.coverCount[start + off]++;
    }

    for (int off : protectOffsets) {
        state.protectCount[start + off]++;
    }

    state.remaining.erase(start);
}

std::vector<int> buildGreedyOrder(int readLength, int span, const std::vector<int>& careOffsets,
                                  const std::vector<int>& protectOffsets)
{
    GreedyState state = initialiseState(readLength, span, protectOffsets);

    while (!state.remaining.empty()) {
        bool coveragePhase = coverageIncomplete(state);

        int start;

        if (coveragePhase) {
            start = chooseBestCoverageSeed(state, careOffsets);
        }
        else {
            start = chooseBestProtectionSeed(state, protectOffsets);
        }

        applySeed(state, start, careOffsets, protectOffsets);
    }

    return state.ordering;
}
/*============= END GREEDY SEED PLACEMENT CODE BLOCK =============*/

std::vector<uint32_t> calculateSeedPositions(int L, int N, int P)
{
    std::vector<uint32_t> positions;

    int totalLength = P * N;
    int buffer = 0;

    if (totalLength + 2 * buffer <= L) {
        int uncoveredLength = L - totalLength - 2 * buffer;
        int gap = uncoveredLength / (P - 1);

        for (int i = 0; i < P; ++i) {
            uint32_t pos = buffer + i * (N + gap);
            positions.push_back(pos);
        }
    }
    else {
        double step = static_cast<double>(L - 2 * buffer - N) / (P - 1);
        for (int i = 0; i < P; ++i) {
            uint32_t pos = static_cast<uint32_t>(buffer + i * step);
            positions.push_back(pos);
        }
    }
    return positions;
}

std::vector<uint32_t> calculateGreedySeedPositions(int L, const StripePattern& pattern, int P)
{
    const auto& order = gGreedyOrders[L];

    const size_t n = std::min(static_cast<size_t>(P), order.size());

    return std::vector<uint32_t>(order.begin(), order.begin() + n);
}

// ================ START GREEDY TABLE CALCULATION ===========
#if 0
void buildGreedySeedTable()
{
    std::cerr
        << "Testing MonteCarlo..."
        << std::endl;

    int seeds =
        MonteCarlo::findMinimumSeeds(
            gPattern,
            1000,
            0.15,
            0.999,
            PlacementStrategy::GREEDY,
            1000);

    std::cerr
        << "Greedy seeds @1000bp,25% = "
        << seeds
        << std::endl;
}
#endif
// ================ END GREEDY LLOKUP SEEDS TABLE CALCULATION ===========

// Forward packing - works for weight <=16 (32-bit) and >16 (64-bit)
inline bool packStripeKmer(std::string_view window, uint64_t& packed)
{
    if ((int)window.size() < gPattern.window)
        return false;

    if (gPattern.use64)
        return gPattern.encode64(window.data(), packed);
    else {
        uint32_t p32 = 0;
        bool ok = gPattern.encode32(window.data(), p32);
        packed = p32;
        return ok;
    }
}

// RC packing
inline bool packStripeKmerRC(const char* wptr, uint64_t& packed)
{
    if ((int)gPattern.window == 0)
        return false; // safety

    if (gPattern.use64)
        return gPattern.encode_rc64(wptr, packed);
    else {
        uint32_t p32 = 0;
        bool ok = gPattern.encode_rc32(wptr, p32);
        packed = p32;
        return ok;
    }
}

void findMatches(const std::span<char>& genomeStr, SafeQueue<Result>& queue, uint32_t minMatches,
                 std::span<const ReadInfo> chunk, const std::span<uint32_t>& innerMapBlob,
                 const std::span<uint32_t>& outerMapBlob, uint32_t MASK, size_t startIndex,
                 const std::shared_ptr<char[]>& data_ptr, const std::vector<ReadInfo>& readInfos)
{
    // We will pass at most ONE candidate to getSW() per read.
    std::vector<std::tuple<uint32_t, uint32_t, uint32_t, uint32_t, uint32_t>> oneCandidate;
    oneCandidate.reserve(1);

    // STRIDE ONLY if explicitly requested. There is NO stride mode in AUTO.
    const bool useStrideSeeding = (SEED_STEP_SPECIFIED);
    uint32_t seedsFound;
    uint32_t kmerHash;

    // POINTER-BASED seed lists for findSpansPtr()
    std::vector<ListStructPtr> kmerPositionsPtr;
    kmerPositionsPtr.reserve(MAX_SEEDS);

    bool debug = false;

    for (size_t idx = 0; idx < chunk.size(); ++idx) {
        const size_t IDX = idx + startIndex;
        (void)IDX;

        const ReadInfo& readInfo = chunk[idx];
        uint32_t sequenceLength = static_cast<uint32_t>(readInfo.length);

        int Min_Matches = MIN_MATCHES;

        if (sequenceLength < ((uint32_t)gPattern.window) * static_cast<uint32_t>(Min_Matches))
            continue;

        // ------------------------------------------------------------
        // SEED LOGIC: priority is
        //   1) stride if user specified -seed_step (handled later via
        //   useStrideSeeding) 2) fixed count if user specified -num_seeds 3) AUTO
        //   fixed-count based on ERROR_TARGET_PCT (no stride in auto)
        //      - <= 1%  : 10 seeds
        //      - <=10%  : 55 seeds
        //      - > 10%  : interpolate from empirical table (capped by what fits)
        //
        // NOTE: The +50% boost for short reads (<=200bp) applies ONLY in AUTO,
        //       exactly as per your current policy.
        // ------------------------------------------------------------
        uint32_t requestedSeeds = 0;

        const uint32_t maxPossibleSeeds = (sequenceLength >= ((uint32_t)gPattern.window))
                                              ? (sequenceLength - ((uint32_t)gPattern.window) + 1)
                                              : 0;

        if (MAX_SEEDS_SPECIFIED) {
            requestedSeeds = maxPossibleSeeds;
        }
        else if (NUM_SEEDS_SPECIFIED) {
            // User explicitly requested a fixed number of seeds (NOT stride).
            requestedSeeds = std::min<uint32_t>(maxPossibleSeeds, static_cast<uint32_t>(NUM_SEEDS));
        }

        else {
            // AUTO: use truncated seed table (version 2 .bin).
            // queryTrunc() returns optimal truncation length l*
            // and minimum seed count n* at that length.
            // Falls back to original query() for version 1 .bin files.
            auto ts = gSeedMatrix.queryTrunc(static_cast<double>(ERROR_TARGET_PCT) / 100.0,
                                             static_cast<int>(sequenceLength));

            // Apply truncation — never exceed actual read length
            const uint32_t effectiveLength =
                std::min(static_cast<uint32_t>(ts.truncLength), sequenceLength);

            const uint32_t maxEffectiveSeeds =
                (effectiveLength >= static_cast<uint32_t>(gPattern.window))
                    ? (effectiveLength - static_cast<uint32_t>(gPattern.window) + 1)
                    : 0u;

            uint32_t n = static_cast<uint32_t>(ts.seedCount);

            if (n < 5u)
                n = std::min(5u, maxEffectiveSeeds);
            if (n > maxEffectiveSeeds)
                n = maxEffectiveSeeds;
            if (n > static_cast<uint32_t>(MAX_SEEDS))
                n = static_cast<uint32_t>(MAX_SEEDS);

            requestedSeeds = n;
            sequenceLength = effectiveLength;
        }
        /*
                else
                {
                    // AUTO (fixed-count): use seed table computed at startup
                    uint32_t n =
                        gSeedMatrix.query(
                        static_cast<float>(
                        ERROR_TARGET_PCT / 100.0),
                        static_cast<uint32_t>(
                        sequenceLength));

                    if (n < 5u)
                        n = std::min(5u, maxPossibleSeeds);
                    if (n > maxPossibleSeeds)
                        n = maxPossibleSeeds;
                    if (n > static_cast<uint32_t>(MAX_SEEDS))
                        n = static_cast<uint32_t>(MAX_SEEDS);
                    requestedSeeds = n;
                }
        */

        std::vector<uint32_t> seedPos;
        seedPos.reserve(requestedSeeds);
        if (useStrideSeeding) {
            uint32_t step = SEED_STEP_SPECIFIED ? static_cast<uint32_t>(SEED_STEP) : 0;

            if (step == 0) {
                // (Should not happen when SEED_STEP_SPECIFIED is true, but preserved
                // behavior.)
                const uint32_t maxStart = (sequenceLength >= ((uint32_t)gPattern.window))
                                              ? (sequenceLength - ((uint32_t)gPattern.window))
                                              : 0;
                step = (requestedSeeds > 1) ? (maxStart / (requestedSeeds - 1)) : maxStart;
                if (step == 0)
                    step = 1;
            }

            const uint32_t maxPossibleSeeds =
                (sequenceLength >= ((uint32_t)gPattern.window))
                    ? (sequenceLength - ((uint32_t)gPattern.window) + 1)
                    : 0;

            const uint32_t seedCap = std::min<uint32_t>(MAX_SEEDS, maxPossibleSeeds);

            for (uint32_t pos = 0; pos + ((uint32_t)gPattern.window) <= sequenceLength;
                 pos += step) {
                seedPos.push_back(pos);
                if (seedPos.size() >= static_cast<size_t>(seedCap))
                    break;
            }

            if (seedPos.size() < static_cast<size_t>(Min_Matches))
                continue;
        }

        else {
            // Fixed-count mode: either explicit -num_seeds,
            // or AUTO fixed-count.

            if (useGreedy) {
                seedPos = calculateGreedySeedPositions(sequenceLength, gPattern,
                                                       static_cast<int>(requestedSeeds));
            }
            else {
                seedPos = calculateSeedPositions(sequenceLength, ((uint32_t)gPattern.window),
                                                 static_cast<int>(requestedSeeds));
            }

            if (seedPos.size() < static_cast<size_t>(Min_Matches)) {
                std::cerr << "Rejected by Min_Matches" << std::endl;

                continue;
            }
        }

        const uint32_t numSeeds = static_cast<uint32_t>(seedPos.size());

        std::string_view read(data_ptr.get() + readInfo.start, readInfo.length);

        // ------------------------------------------------------------
        // FORWARD: build pointer lists -> findSpansPtr -> take FIRST span
        // ------------------------------------------------------------
        bool haveF = false;
        uint32_t f_ref = 0;
        uint32_t f_cov = 0;
        uint32_t f_span = 0;

        kmerPositionsPtr.clear();
        seedsFound = 0;

        for (uint32_t seedIndex = 0; seedIndex < numSeeds; ++seedIndex) {
            std::string_view window(read.data() + seedPos[seedIndex], (uint32_t)gPattern.window);

            uint64_t packed = 0;
            if (!packStripeKmer(window, packed)) {
                kmerPositionsPtr.push_back(ListStructPtr{nullptr, 0});
                continue;
            }

            // Correct for any stripe width (16 or 24+ care bases)
            const bool needsHash = (gPattern.weight > 16) || (MASK != UINT32_MAX);
            kmerHash = needsHash ? murmurHash3(packed) & MASK : static_cast<uint32_t>(packed);

            bool ok = getInvertedListsPtr(kmerPositionsPtr, outerMapBlob, innerMapBlob, kmerHash);
            if (ok)
                ++seedsFound;
            (void)ok;
        }

        if (seedsFound > 2) // almost always; 16mers are common
        {
            auto spans_fwd =
                findSpansPtrWithSpan(kmerPositionsPtr, seedPos, static_cast<uint32_t>(Min_Matches),
                                     sequenceLength, INDEL_SLACK);

            if (!spans_fwd.empty()) {
                std::tie(f_ref, f_cov, f_span) = spans_fwd[0];
                haveF = true;
            }
        }

        // ------------------------------------------------------------
        // SIMPLE GATE: if forward best coverage >= GATE_COVERAGE -> skip RC
        // entirely
        // ------------------------------------------------------------
        bool haveR = false;
        uint32_t r_ref = 0;
        uint32_t r_cov = 0;
        uint32_t r_span = 0;

        if (!(haveF && f_cov >= GATE_COVERAGE)) {
            // ------------------------------------------------------------
            // REVERSE: pointer lists + findSpansPtr + FIRST span only
            // ------------------------------------------------------------
            kmerPositionsPtr.clear();
            seedsFound = 0;
            const uint32_t maxStart = (sequenceLength >= ((uint32_t)gPattern.window))
                                          ? (sequenceLength - ((uint32_t)gPattern.window))
                                          : 0;
            //////////////////////////////////////////

            for (uint32_t seedIndex = 0; seedIndex < numSeeds; ++seedIndex) {
                const uint32_t rc_pos = seedPos[seedIndex];

                if (rc_pos > maxStart) {
                    kmerPositionsPtr.push_back(ListStructPtr{nullptr, 0});
                    continue;
                }

                const uint32_t fwd_pos = maxStart - rc_pos;
                const char* wptr = read.data() + fwd_pos;

                uint64_t packed = 0;
                if (!packStripeKmerRC(wptr, packed)) {
                    if (debug) {
                        std::string_view dbg_window(wptr, (uint32_t)gPattern.window);
                        std::cout << "Rejected RC-striped seed window (ambiguous?) " << dbg_window
                                  << "\n";
                    }
                    kmerPositionsPtr.push_back(ListStructPtr{nullptr, 0});
                    continue;
                }

                const bool needsHash = (gPattern.weight > 16) || (MASK != UINT32_MAX);
                kmerHash = needsHash ? murmurHash3(packed) & MASK : static_cast<uint32_t>(packed);

                bool ok =
                    getInvertedListsPtr(kmerPositionsPtr, outerMapBlob, innerMapBlob, kmerHash);
                if (ok)
                    ++seedsFound;
                (void)ok;
            }

            //////////////////////////////////////////////
            if (seedsFound > 2) // almost always; 16mers are common
            {
                auto spans_rev = findSpansPtrWithSpan(kmerPositionsPtr, seedPos,
                                                      static_cast<uint32_t>(Min_Matches),
                                                      sequenceLength, INDEL_SLACK);

                if (!spans_rev.empty()) {
                    std::tie(r_ref, r_cov, r_span) = spans_rev[0];
                    haveR = true;
                }
            }
        }
        // ------------------------------------------------------------
        // PICK BEST (coverage first, then smaller refBegin)
        // ------------------------------------------------------------
        bool haveBest = false;
        uint32_t bestRefBegin = 0;
        uint32_t bestCoverage = 0;
        uint32_t bestRevComp = 0;
        uint32_t bestSpan = 0;

        if (haveF) {
            haveBest = true;
            bestRefBegin = f_ref;
            bestCoverage = f_cov;
            bestRevComp = 0;
            bestSpan = f_span;
        }

        if (haveR) {
            if (!haveBest || (r_cov > bestCoverage) ||
                (r_cov == bestCoverage && r_ref < bestRefBegin)) {
                haveBest = true;
                bestRefBegin = r_ref;
                bestCoverage = r_cov;
                bestRevComp = 1;
                bestSpan = r_span;
            }
        }

        // FINAL GATE (same threshold)
        if (DEBUG_COVERAGE) {
            std::lock_guard<std::mutex> lock(logMutex);
            std::cerr << "SEEDS\t" << (startIndex + idx) << "\t" << seedsFound << "\t"
                      << (haveBest ? (int)bestCoverage : -1) << "\t"
                      << ((haveBest && bestCoverage >= GATE_COVERAGE) ? 1 : 0) << "\n";
        }
        if (!haveBest || bestCoverage < GATE_COVERAGE)
            continue;

        // ------------------------------------------------------------
        // SINGLE getSW CALL (exactly once per read)
        // ------------------------------------------------------------
        oneCandidate.clear();
        oneCandidate.emplace_back(static_cast<uint32_t>(idx), bestRefBegin, bestCoverage, bestSpan,
                                  bestRevComp);

        int matchScore = 1;
        int mismatchScore = -1;
        int gapScore = -1;

        getSW(data_ptr.get(), readInfos, startIndex, genomeStr, queue, oneCandidate, matchScore,
              mismatchScore, gapScore, CUTOFF);
    }
}

// Strip non-ASCII and trailing whitespace from path/string arguments.
// Guards against invisible Unicode characters (e.g. non-breaking space U+00A0)
// that can appear when commands are copied from terminal logs or web pages.
/*
static std::string cleanArg(const std::string &s)
{
    std::string out;
    out.reserve(s.size());
    for (unsigned char c : s)
        if (c >= 0x20 && c <= 0x7e)   // printable ASCII only
            out += (char)c;
    while (!out.empty() && out.back() == ' ')
        out.pop_back();
    return out;
}
*/
// Replace non-ASCII / non-printable characters with a single space,
// collapse repeated whitespace, and trim leading/trailing spaces.
// Prevents invisible Unicode characters from merging tokens together.

static std::string cleanArg(const std::string& s)
{
    std::string out;
    out.reserve(s.size());

    bool prevSpace = false;

    for (unsigned char c : s) {
        if (c >= 0x21 && c <= 0x7e) // printable ASCII excluding space
        {
            out += (char)c;
            prevSpace = false;
        }
        else {
            // Convert ANY non-printable / non-ASCII / whitespace
            // into a single ordinary ASCII space.
            if (!prevSpace) {
                out += ' ';
                prevSpace = true;
            }
        }
    }

    // Trim leading space
    if (!out.empty() && out.front() == ' ')
        out.erase(out.begin());

    // Trim trailing space
    if (!out.empty() && out.back() == ' ')
        out.pop_back();

    return out;
}

void initialiseSKIPPY(int argc, char* argv[])
{
    std::string arg;
    if (argc == 1) {
        cerr << "Error: missing search options; see skippy search --help\n";
        exit(2);
    }
    if (argc > 1) {
        arg = argv[1];
    }
    if (arg == "-help" || arg == "--help") {
        bool fullHelp = (argc > 2 && std::string(argv[2]) == "all");

        cout << "\nSkippy -- fast probabilistic reference sequence detection and "
                "depletion\n";
        cout << "\nUsage: skippy search -genome <dir> -technology <tech> -data "
                "<reads.fasta> -out_dir <dir> [options]\n";
        cout << "\nRequired:\n";
        cout << "  -genome <dir>          reference reference genome root "
                "directory\n";
        cout << "  -data <file>           input reads in FASTA or FASTQ format\n";
        cout << "  -out_dir <dir>         output directory\n";
        cout << "  -format <out|sam>      output format (default: out)\n";
        cout << "\nTechnology (sets index, contract and error rate):\n";
        cout << "  -technology <tech>     illumina | hifi | pacbio | nanopore | "
                "ont\n";
        cout << "\nOperating mode (choose one):\n";
        cout << "  --sensitivity          maximum sensitivity (default)\n";
        cout << "  --balanced             balanced sensitivity and specificity\n";
        cout << "  --specificity          maximum specificity, minimum false "
                "positives\n";
        cout << "\nCommon options:\n";
        cout << "  -error_target <pct>    override error rate set by -technology "
                "(%)\n";
        cout << "  -threads <n>           number of threads (default: all "
                "available)\n";
        cout << "  --help                 show this help\n";
        cout << "  --version              show version\n";

        if (fullHelp) {
            cout << "\nIndex and contract:\n";
            cout << "  -stripe <name>              index stripe pattern (plain | "
                    "4XX4XXX1XX3X4)\n";
            cout << "  -contract <type>            seed contract type (subs | "
                    "indels)\n";
            cout << "\nTuning:\n";
            cout << "  -gate <n>                   minimum seed coverage gate "
                    "(default: 33)\n";
            cout << "  -cut_off <n>                minimum SW matching bases "
                    "(default: 0)\n";
            cout << "  -slack <n>                  indel chain slack 0-10 (default: "
                    "10)\n";
            cout << "  -identity_threshold <f>     SW identity filter (default: "
                    "0.63)\n";
            cout << "  -suspicious_span <n>        span threshold for SW (default: "
                    "200)\n";
            cout << "  -num_seeds <n>              fix seeds per read (minimum 5)\n";
            cout << "  -sw <true|false>            compute SW alignment (default: "
                    "false; enabled for SAM)\n";
            cout << "\nPower user flags (not shown in basic help):\n";
            cout << "  -@gate, -@cut_off, -@identity_threshold, -@suspicious_span\n";
            cout << "\nExamples:\n";
            cout << "  skippy search -genome /data/GRCh38 -technology nanopore -data "
                    "reads.fasta -out_dir results\n";
            cout << "  skippy search -genome /data/GRCh38 -technology nanopore -data "
                    "reads.fasta -out_dir results --specificity\n";
        }
        else {
            cout << "  For full options: skippy search --help all\n";
        }
        cout << endl;
        exit(0);
    }

    ERROR_TARGET_PCT = 10; // default is 10% base error
                           //    KMERSIZE = 32;
    MIN_MATCHES = 2;
    THREADS = std::thread::hardware_concurrency();
    CUTOFF = 0;
    REFERENCE = "";
    READS = "";
    OUTPUT_DIR = get_current_dir_name_custom() + "/save";
    MATCH_ALL = false;
    GET_SW = false;
    INDEL_SLACK = 10;
    GATE_COVERAGE = 33; // 40?
    for (int i = 1; i < argc; ++i) {
        arg = argv[i];
        /* --------------------------------------*/
        // Boolean flags with no value
        if (arg == "-maxseeds") {
            MAX_SEEDS_SPECIFIED = true;
            if (NUM_SEEDS_SPECIFIED || ERROR_TARGET_SPECIFIED || TECHNOLOGY_SPECIFIED)
                std::cerr << "WARNING: -maxseeds overrides all other seed count settings.\n";
            continue;
        }
        if (arg == "-greedy") {
            useGreedy = true;
            continue;
        }
        if (arg == "-skippy") {
            useGreedy = false;
            continue;
        }
        if (arg == "--sensitivity") {
            GATE_COVERAGE = 33;
            CUTOFF = 0;
            IDENTITY_THRESHOLD = 0.50f;
            SUSPICIOUS_SPAN_THRESHOLD = 200;
            continue;
        }
        if (arg == "--specificity") {
            GATE_COVERAGE = 49;
            CUTOFF = 45;
            IDENTITY_THRESHOLD = 0.65f;
            SUSPICIOUS_SPAN_THRESHOLD = 200;
            continue;
        }
        if (arg == "--balanced") {
            GATE_COVERAGE = 33;
            CUTOFF = 0;
            IDENTITY_THRESHOLD = 0.50f;
            SUSPICIOUS_SPAN_THRESHOLD = 200;
            BALANCED_MODE = true;
            continue;
        }
        static const std::set<std::string> valueOptions = {
            "-min_matches", "-contract", "-num_seeds", "-seed_step", "-slack", "-gate",
            "-technology", "-error_target", "-cut_off", "-@gate", "-@cut_off",
            "-@identity_threshold", "-@suspicious_span", "-identity_threshold",
            "-suspicious_span", "-threads", "-reference", "-genome", "-stripe", "-data",
            "-out_dir", "-format", "-match_all", "-sw", "-paired", "-debug_coverage", "-min_span",
            "-min_density", "-min_matches_abs"};
        if (!valueOptions.count(arg)) {
            cerr << "Error: unknown search option: " << arg << '\n';
            exit(2);
        }
        /* --------------------------------------*/
        if (i + 1 >= argc ||
            (argv[i + 1][0] == '-' && (argv[i + 1][1] < '0' || argv[i + 1][1] > '9'))) {
            cerr << "Error: missing value for " << arg << "; see skippy search --help\n";
            exit(2);
        }

        std::string value = argv[i + 1];
        ++i;

        try {
            if (arg == "-min_matches") {
                MIN_MATCHES = std::stoi(value);
                if (MIN_MATCHES < 2)
                    MIN_MATCHES = 2;
            }
            else if (arg == "-contract") {
                CONTRACT_TYPE = cleanArg(value);
                CONTRACT_SPECIFIED = true;
                std::cerr << "DEBUG: -contract set to '" << CONTRACT_TYPE << "'\n";
            }
            else if (arg == "-num_seeds") {
                NUM_SEEDS = std::stoi(value);
                if (NUM_SEEDS < 5) {
                    std::cerr << "WARNING: -num_seeds " << NUM_SEEDS
                              << " is below minimum of 5. Using 5.\n";
                    NUM_SEEDS = 5;
                }
                NUM_SEEDS_SPECIFIED = true;
            }
            else if (arg == "-seed_step") {
                SEED_STEP = std::stoi(value);
                if (SEED_STEP < 1)
                    SEED_STEP = 1;
                SEED_STEP_SPECIFIED = true;
            }
            else if (arg == "-slack") {
                const int slack = std::stoi(value);
                INDEL_SLACK = (slack < 0 || slack > 30) ? 3u : static_cast<uint32_t>(slack);
            }
            else if (arg == "-gate") {
                const int gate = std::stoi(value);
                if (gate < 0) {
                    cerr << "Error: -gate must be nonnegative\n";
                    exit(2);
                }
                GATE_COVERAGE = static_cast<uint32_t>(gate);
            }

            else if (arg == "-technology") {
                TECHNOLOGY_NAME = cleanArg(value);
                std::transform(TECHNOLOGY_NAME.begin(), TECHNOLOGY_NAME.end(), TECHNOLOGY_NAME.begin(),
                               ::tolower);
                if (TECHNOLOGY_NAME == "illumina")
                    ERROR_TARGET_PCT = 6;
                else if (TECHNOLOGY_NAME == "hifi")
                    ERROR_TARGET_PCT = 6;
                else if (TECHNOLOGY_NAME == "pacbio")
                    ERROR_TARGET_PCT = 15;
                else if (TECHNOLOGY_NAME == "nanopore")
                    ERROR_TARGET_PCT = 15;
                else if (TECHNOLOGY_NAME == "ont")
                    ERROR_TARGET_PCT = 15;
                else {
                    cerr << "Error: unknown technology '" << TECHNOLOGY_NAME
                         << "'; expected illumina, hifi, pacbio, nanopore, or ont\n";
                    exit(2);
                }
                TECHNOLOGY_SPECIFIED = true;

                if (TECHNOLOGY_NAME == "illumina") {
                    CONTRACT_TYPE = "subs";
                }
                else if (TECHNOLOGY_NAME == "hifi") {
                    CONTRACT_TYPE = "subs";
                }
                else if (TECHNOLOGY_NAME == "pacbio") {
                    CONTRACT_TYPE = "indels";
                }
                else if (TECHNOLOGY_NAME == "nanopore") {
                    CONTRACT_TYPE = "indels";
                }
                else if (TECHNOLOGY_NAME == "ont") {
                    CONTRACT_TYPE = "indels";
                }
            }
            else if (arg == "-error_target") {
                ERROR_TARGET_PCT = std::stoi(value);
                if (TECHNOLOGY_SPECIFIED) {
                    std::cerr << "WARNING: -error_target " << value
                              << " overrides the default error target set by -technology "
                              << TECHNOLOGY_NAME << " (" << ERROR_TARGET_PCT << "%).\n";
                }
                ERROR_TARGET_SPECIFIED = true;
            }
            else if (arg == "-cut_off") {
                CUTOFF = std::stoi(value);
                if (CUTOFF > 100)
                    CUTOFF = 100;
                else if (CUTOFF < 0)
                    CUTOFF = 0;
            }
            else if (arg == "-@gate") {
                GATE_COVERAGE = std::stoi(value);
            }
            else if (arg == "-@cut_off") {
                CUTOFF = std::stoi(value);
            }
            else if (arg == "-@identity_threshold") {
                IDENTITY_THRESHOLD = std::stof(value);
            }
            else if (arg == "-@suspicious_span") {
                SUSPICIOUS_SPAN_THRESHOLD = static_cast<uint32_t>(std::stoi(value));
            }
            else if (arg == "-identity_threshold") {
                IDENTITY_THRESHOLD = std::stof(value);
                if (IDENTITY_THRESHOLD > 1.0f)
                    IDENTITY_THRESHOLD = 1.0f;
                else if (IDENTITY_THRESHOLD < 0.0f)
                    IDENTITY_THRESHOLD = 0.0f;
            }
            else if (arg == "-suspicious_span") {
                SUSPICIOUS_SPAN_THRESHOLD = static_cast<uint32_t>(std::stoi(value));
            }
            else if (arg == "-@gate") {
                GATE_COVERAGE = std::stoi(value);
            }
            else if (arg == "-@cut_off") {
                CUTOFF = std::stoi(value);
            }
            else if (arg == "-@identity_threshold") {
                IDENTITY_THRESHOLD = std::stof(value);
            }
            else if (arg == "-@suspicious_span") {
                SUSPICIOUS_SPAN_THRESHOLD = static_cast<uint32_t>(std::stoi(value));
            }
            else if (arg == "-identity_threshold") {
                IDENTITY_THRESHOLD = std::stof(value);
                if (IDENTITY_THRESHOLD > 1.0f)
                    IDENTITY_THRESHOLD = 1.0f;
                else if (IDENTITY_THRESHOLD < 0.0f)
                    IDENTITY_THRESHOLD = 0.0f;
            }
            else if (arg == "-suspicious_span") {
                SUSPICIOUS_SPAN_THRESHOLD = static_cast<uint32_t>(std::stoi(value));
            }
            else if (arg == "-threads") {
                int threads = std::stoi(value);
                if (threads < 1) {
                    cerr << "Error: -threads must be at least 1\n";
                    exit(2);
                }
                THREADS = threads;
            }
            else if (arg == "-reference") {
                REFERENCE = cleanArg(value);
            }
            else if (arg == "-genome") {
                GENOME_DIR = cleanArg(value);
            }
            else if (arg == "-stripe") {
                STRIPE_NAME = cleanArg(value);
            }
            else if (arg == "-data") {
                READS = cleanArg(value);
            }
            else if (arg == "-out_dir") {
                OUTPUT_DIR = cleanArg(value);
            }
            else if (arg == "-format") {
                if (value != "out" && value != "sam") {
                    cerr << "Error: -format must be out or sam\n";
                    exit(2);
                }
                SAM_OUTPUT = (value == "sam");
            }
            else if (arg == "-match_all") {
                MATCH_ALL = (value == "true");
            }
            else if (arg == "-sw") {
                GET_SW = (value == "true");
            }
            else if (arg == "-paired") {
                PAIRED_READS = (value == "true");
            }
            else if (arg == "-debug_coverage") {
                DEBUG_COVERAGE = (value == "true");
            }
#ifdef USE_ALIGNMENT_DENSITY
            else if (arg == "-min_span") {
                MIN_ALIGNMENT_SPAN = std::stoi(value);
                if (MIN_ALIGNMENT_SPAN < 1)
                    MIN_ALIGNMENT_SPAN = 1;
            }
            else if (arg == "-min_density") {
                MIN_ALIGNMENT_DENSITY = std::stof(value);
                if (MIN_ALIGNMENT_DENSITY > 1.0f)
                    MIN_ALIGNMENT_DENSITY = 1.0f;
                else if (MIN_ALIGNMENT_DENSITY < 0.0f)
                    MIN_ALIGNMENT_DENSITY = 0.0f;
            }
            else if (arg == "-min_matches_abs") {
                MIN_ABSOLUTE_MATCHES = std::stoi(value);
                if (MIN_ABSOLUTE_MATCHES < 1)
                    MIN_ABSOLUTE_MATCHES = 1;
            }
#endif

            else {
                cerr << "Error: unknown search option: " << arg << '\n';
                exit(2);
            }
        }
        catch (const std::invalid_argument&) {
            cerr << "Error: invalid value for " << arg << ": " << value << '\n';
            exit(2);
        }
        catch (const std::out_of_range&) {
            cerr << "Error: value out of range for " << arg << ": " << value << '\n';
            exit(2);
        }
    }

    if (READS.empty()) {
        cerr << "Error: missing -data <file>; see skippy search --help\n";
        exit(2);
    }
    if (GENOME_DIR.empty() && REFERENCE.empty()) {
        cerr << "Error: missing -genome <dir> or -reference <fasta>; see skippy search --help\n";
        exit(2);
    }
    if (CONTRACT_TYPE != "subs" && CONTRACT_TYPE != "indels") {
        cerr << "Error: invalid -contract value: " << CONTRACT_TYPE << " (expected subs or indels)\n";
        exit(2);
    }
    if (SAM_OUTPUT && PAIRED_READS) {
        cerr << "Error: paired SAM output is not supported\n";
        exit(2);
    }
    if (SAM_OUTPUT)
        GET_SW = true;

    kmer_encoding_table[(size_t)'A'] = 0;
    kmer_encoding_table[(size_t)'C'] = 1;
    kmer_encoding_table[(size_t)'G'] = 2;
    kmer_encoding_table[(size_t)'T'] = 3;

    kmer_encoding_table[(size_t)'a'] = 0;
    kmer_encoding_table[(size_t)'c'] = 1;
    kmer_encoding_table[(size_t)'g'] = 2;
    kmer_encoding_table[(size_t)'t'] = 3;

    std::string rootDir = READS;
    fileNames = getFileNames(rootDir);

    if (fileNames.empty()) {
        cerr << "Sample File: " + READS + ", does not exist in this directory " << endl;
        exit(1);
    }

    std::cout << "Processing " << fileNames.size() << " files" << std::endl;

    cout << "SKIPPY run parameters" << endl;
    cout << "min_matches: "
         << (MIN_MATCHES > 0 ? std::to_string(MIN_MATCHES) : "Dynamic (3 + read_length / 100")
         << endl;
    cout << "num-seeds: ";
    if (NUM_SEEDS_SPECIFIED)
        cout << std::to_string(NUM_SEEDS) << endl;
    else if (SEED_STEP_SPECIFIED)
        cout << "Dynamic, step size : " << std::to_string(SEED_STEP) << endl;
    else
        cout << "Auto" << endl;
    cout << "threads: " << THREADS << endl;
    cout << "coverage gate: " << GATE_COVERAGE << endl;
    cout << "cut_off: " << CUTOFF << endl;
    cout << "reference: " << REFERENCE << endl;
    cout << "compute Smith Waterman: " << (GET_SW == true ? "TRUE" : "FALSE") << endl;
    cout << "indel_slack: " << INDEL_SLACK << endl;
    cout << "error_target: " << ERROR_TARGET_PCT << endl;
    cout << "data: " << READS << endl;
    cout << "output directory: " << OUTPUT_DIR << endl << endl;
    cout << "output format: " << (SAM_OUTPUT ? "sam" : "out") << endl;
#ifdef USE_ALIGNMENT_DENSITY
    cout << "min_alignment_span: " << MIN_ALIGNMENT_SPAN << endl;
    cout << "min_alignment_density: " << MIN_ALIGNMENT_DENSITY << endl;
    cout << "min_absolute_matches: " << MIN_ABSOLUTE_MATCHES << endl;
#else
    cout << "cut_off: " << CUTOFF << endl;
#endif
    if (TECHNOLOGY_SPECIFIED)
        cout << "technology:              " << TECHNOLOGY_NAME << "\n"
             << "sequencing error rate:   " << ERROR_TARGET_PCT << "% (set by -technology)\n";
    else if (ERROR_TARGET_SPECIFIED)
        cout << "sequencing error rate:   " << ERROR_TARGET_PCT << "% (set by -error_target)\n";
    else
        cout << "sequencing error rate:   " << ERROR_TARGET_PCT
             << "% (default — use -technology or -error_target to override)\n";
    cout << "paired reads : " << PAIRED_READS << "\n";
    //------------------------------------------------------------------
    // Resolve -genome / -stripe / -technology into REFERENCE path
    //------------------------------------------------------------------
    if (!GENOME_DIR.empty()) {
        // Derive genome name from directory basename
        std::string genomeName = GENOME_DIR;
        size_t slash = genomeName.find_last_of("/\\");
        if (slash != std::string::npos)
            genomeName = genomeName.substr(slash + 1);

        std::string fastaFile = GENOME_DIR + "/" + genomeName + ".fasta";

        // Determine stripe subdirectory
        // Priority: -stripe > -technology > default (illumina symlink)
        std::string stripeSubdir;
        if (!STRIPE_NAME.empty()) {
            stripeSubdir = STRIPE_NAME;
        }
        else if (TECHNOLOGY_SPECIFIED) {
            stripeSubdir = TECHNOLOGY_NAME; // symlink name, e.g. "nanopore"
        }
        else {
            stripeSubdir = "illumina"; // default
        }

        std::string indexDir = GENOME_DIR + "/" + stripeSubdir;

        // Check FASTA exists
        struct stat fst;
        if (stat(fastaFile.c_str(), &fst) != 0) {
            cerr << "Error: Genome FASTA not found: " << fastaFile << "\n";
            cerr << "Please place the genome at <genome_dir>/<name>.fasta\n";
            exit(1);
        }

        // Check index directory exists
        struct stat dst;
        if (stat(indexDir.c_str(), &dst) != 0 || !S_ISDIR(dst.st_mode)) {
            // Try falling back to plain if requested stripe is a technology symlink
            std::string plainDir = GENOME_DIR + "/plain";
            struct stat pst;
            bool plainExists = (stat(plainDir.c_str(), &pst) == 0 && S_ISDIR(pst.st_mode));

            cerr << "Error: Index '" << stripeSubdir << "' not found for " << genomeName << ".\n";
            cerr << "Expected: " << indexDir << "\n\n";
            if (!STRIPE_NAME.empty()) {
                cerr << "To build this index, run:\n";
                cerr << "  skippy index -genome " << GENOME_DIR << " -stripe " << STRIPE_NAME
                     << "\n";
            }
            else {
                if (plainExists) {
                    cerr << "Note: A plain index exists. To use it, add: -stripe plain\n";
                    cerr << "Or build the " << stripeSubdir << " index:\n";
                    cerr << "  skippy index -genome " << GENOME_DIR << " -stripe 4XX4XXX1XX3X4\n";
                }
                else {
                    cerr << "To build a plain index (works for all technologies), run:\n";
                    cerr << "  skippy index -genome " << GENOME_DIR << " -stripe plain\n";
                }
            }
            exit(1);
        }

        // Resolve symlink to get actual stripe directory name for index file paths
        char resolved[PATH_MAX];
        std::string resolvedIndexDir = indexDir;
        if (realpath(indexDir.c_str(), resolved) != nullptr)
            resolvedIndexDir = std::string(resolved);

        // Get the actual genome name from the resolved index dir
        // Index files are: <indexDir>/<genomeName>_OuterBlob.idx etc.
        REFERENCE = resolvedIndexDir + "/" + genomeName + ".fasta";

        cout << "Genome:    " << GENOME_DIR << "\n";
        cout << "Index:     " << resolvedIndexDir << "\n";
        cout << "FASTA:     " << fastaFile << "\n\n";

        // Override REFERENCE to point to the FASTA, but index files come from
        // resolvedIndexDir The existing code uses REFERENCE's parent_path() to find
        // index files, so set REFERENCE to a path whose parent is resolvedIndexDir
        REFERENCE = resolvedIndexDir + "/" + genomeName + ".fasta";
    }

    if (REFERENCE.empty() || READS.empty()) {
        cerr << "Reference file name and Data file name must be provided" << endl;
        exit(1);
    }

    // Reload contract with final CONTRACT_TYPE after all args parsed
    {
        bool ui = (CONTRACT_TYPE == "indels");
        std::filesystem::path rp(REFERENCE);
        std::string sf = (rp.parent_path() / (rp.stem().string() + (ui ? "_indels" : "") +
                                              (useGreedy ? "_greedy" : "_skippy") + ".bin"))
                             .string();
        gSeedMatrix.clear();
        if (!gSeedMatrix.load(sf)) {
            std::cerr << "Failed to load " << sf << std::endl;
            exit(1);
        }
        else
            std::cerr << "Successful load " << sf << std::endl;
        cerr << "contract type: " << CONTRACT_TYPE << endl;
    }
    cerr << "error_target: " << ERROR_TARGET_PCT << endl;
}

// main (GLOBAL WORK QUEUE + GLOBAL WORKERS + PER-FILE WRITER THREADS)
int main(int argc, char* argv[])
{
    std::span<uint32_t> innerMapBlob;
    std::span<uint32_t> outerMapBlob;
    std::span<char> genomeStr;
    uint32_t MASK;

    auto startAll = std::chrono::steady_clock::now();

    initialiseSKIPPY(argc, argv);
    cout << "\n********************* initialise complete *******************\n";
    // Bound the global work queue to enforce back-pressure.
    // Capacity is "number of chunks", not reads.
    // A good conservative default is 8 chunks per worker.
    workQueue.setCapacity(static_cast<size_t>(THREADS) * 8);

    getReference(REFERENCE, innerMapBlob, outerMapBlob, genomeStr, MASK, ref_ID_Table);
    if (SAM_OUTPUT && ref_ID_Table.entries.empty()) {
        std::cerr << "Error: SAM output requires reference names in the index\n";
        return 1;
    }

    // -------------------------
    // Go to work now
    // -------------------------
    cout << "starting threads\n";
    std::thread loggingThread([&] {
        while (true) {
            std::unique_lock<std::mutex> lock(logMutex);

            logCV.wait(lock, [&] {
                return loggingFinished || fileLogMessages.count(nextFileToPrint) > 0 ||
                       !fileLogMessages.empty();
            });

            while (fileLogMessages.count(nextFileToPrint) > 0) {
                auto& messages = fileLogMessages[nextFileToPrint];
                auto it = std::find_if(messages.begin(), messages.end(), [](const std::string& s) {
                    return s.find("[LAST]") != std::string::npos;
                });
                if (it == messages.end())
                    break;

                for (const auto& msg : messages)
                    std::cout << msg << std::endl;

                fileLogMessages.erase(nextFileToPrint);
                ++nextFileToPrint;
            }

            if (loggingFinished && nextFileToPrint >= fileNames.size()) {
                for (auto& kv : fileLogMessages) {
                    for (const auto& msg : kv.second)
                        std::cout << msg << std::endl;
                }
                fileLogMessages.clear();
                break;
            }
        }
    });

    // -------------------------
    // GLOBAL WORKERS
    // -------------------------

    auto workerLoop = [&]() {
        WorkItem item;
        while (workQueue.dequeue(item)) {
            findMatches(genomeStr, item.ctx->writerQueue, MIN_MATCHES, item.chunk, innerMapBlob,
                        outerMapBlob, MASK, item.startIndex, item.chunkData->data_ptr,
                        item.chunkData->readInfos);

            // Completion signalling for streaming:
            // decrement inFlight; if producer is done and this was the last in-flight
            // item, finish writer queue.
            if (item.ctx->inFlight.fetch_sub(1, std::memory_order_acq_rel) == 1) {
                if (item.ctx->doneEnqueue.load(std::memory_order_acquire)) {
                    item.ctx->writerQueue.setFinished();
                }
            }
        }
    };

    // NOTE: This is your "global worker pool". We do NOT use ThreadPool here,
    // because these workers are blocking on a queue.
    std::vector<std::thread> workers;
    workers.reserve(THREADS);
    for (uint32_t t = 0; t < THREADS; ++t)
        workers.emplace_back(workerLoop);

    // -------------------------
    // PRODUCER: loads files, starts writer threads, enqueues chunks to GLOBAL
    // queue
    // -------------------------
    std::vector<std::shared_ptr<FileContext>> allFileContexts;
    allFileContexts.reserve(fileNames.size());

    std::thread producer([&] {
        for (size_t fileNum = 0; fileNum < fileNames.size(); ++fileNum) {
            const std::string& filename = fileNames[fileNum];
            auto ctx = std::make_shared<FileContext>();
            allFileContexts.push_back(ctx);

            ctx->fileNum = fileNum;

            logForFile(fileNum,
                       std::to_string(fileNum) + " ******** processing " + filename + " ********");

            std::string baseFilename = getBaseFilename(filename);
            ctx->outFilename = OUTPUT_DIR + "/" + baseFilename + (SAM_OUTPUT ? ".sam" : ".out");

            logForFile(fileNum, "Results written to " + ctx->outFilename);

            // Start writer thread immediately (it blocks on writerQueue)
            ctx->writerThread = std::thread(writer, std::ref(*ctx), genomeStr.size());

            // Stream chunks into the global workQueue (mmap + incremental FASTA scan)
            uint64_t fileSize = 0;
            enqueue_fasta_or_fastq_streaming_chunks(filename, ctx, fileSize, fileNum, workQueue);

            // Producer is done for this file.
            ctx->doneEnqueue.store(true, std::memory_order_release);

            // If nothing is in flight, release writer immediately (e.g., empty/no
            // valid reads)
            if (ctx->inFlight.load(std::memory_order_acquire) == 0)
                ctx->writerQueue.setFinished();

            logForFile(fileNum, std::to_string(ctx->totalReads.load(std::memory_order_relaxed)) +
                                    " samples enqueued");
        }
    });

    // Wait for producer to finish enqueuing everything
    producer.join();

    // Now no more work items will be enqueued
    workQueue.setFinished();

    // Join global workers
    for (auto& th : workers)
        th.join();

    // Join all writers (they should have been released by remainingChunks hitting
    // 0)
    bool writeFailed = false;
    for (auto& ctx : allFileContexts) {
        if (ctx && ctx->writerThread.joinable())
            ctx->writerThread.join();
        if (ctx)
            writeFailed |= ctx->writeFailed;
    }

    // Shutdown logging
    {
        std::lock_guard<std::mutex> lock(logMutex);
        loggingFinished = true;
    }
    logCV.notify_all();
    loggingThread.join();

    std::cout << "Slack = " << INDEL_SLACK << endl;
    auto endAll = std::chrono::steady_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(endAll - startAll);
    int minutes = (int)duration.count() / (1000 * 60);
    float seconds = (duration.count() - minutes * 1000 * 60) / 1000.0f;

    std::cout << "Total time: " << minutes << " min " << seconds << " sec" << std::endl;

    return writeFailed ? 1 : 0;
}
