/*
        INDEXGENOME.CPP
        ---------------
        indexSKIPPY

        Created by Shlomo Geva on 13/7/2023.
        Updated: runtime stripe pattern (StripePattern), 32/64-bit packed
   k-mers.
*/

#include "indexGenome.hpp"

#include <chrono>
#include <iomanip>
#include <iostream>
#include <thread>

#include "encode_kmer_2bit.h"
#include "hash.hpp"
#include "packGenomeBlob.hpp"
#include "stripe_pattern.h"

/*
        DISPLAYPROGRESS()
        -----------------
*/
void displayProgress(std::chrono::time_point<std::chrono::steady_clock>& start,
                     uint64_t& lastDisplayedPercent, uint64_t current, uint64_t total,
                     int desiredUpdateInterval)
{
    uint64_t percent = (current * 100) / total;
    if (percent - lastDisplayedPercent >= (uint64_t)desiredUpdateInterval) {
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(
                            std::chrono::steady_clock::now() - start)
                            .count();
        std::cout << "progress:" << percent << " (" << duration << " milliseconds)\n";
        lastDisplayedPercent = percent;
    }
}

/*
    READ_ENTIRE_FILE()
    ------------------
*/
char* read_entire_file(const char* filename, uint64_t& fileSize)
{
    FILE* fp;
    struct stat details;
    char* contents = NULL;
    fileSize = 0;
    if ((fp = fopen(filename, "rb")) != NULL) {
        if (fstat(fileno(fp), &details) == 0) {
            if (details.st_size != 0 || details.st_size > UINT32_MAX) {
                contents = (char*)malloc(details.st_size + 1);
                if (fread(contents, details.st_size, 1, fp) != 1) {
                    free(contents);
                    contents = NULL;
                }
                else {
                    fileSize = details.st_size;
                    contents[details.st_size] = '\0';
                }
            }
        }
        fclose(fp);
    }
    return contents;
}

/*
        LOAD_GENOME_FILE()
        ------------------
*/
char* load_genome_file(const std::string& fastaFile,
                       std::map<uint32_t, std::string>& referenceIDMap, uint32_t& genomeSize)
{
    std::cout << "\nLoading References: " << fastaFile << "\n";
    uint64_t fileSize;
    char* genome = read_entire_file(fastaFile.c_str(), fileSize);
    if (genome == NULL) {
        std::cerr << "Failed to read " << fastaFile << " - Either missing, or larger than 4GB\n";
        exit(1);
    }
    std::cout << "Reference file size on disk " << fileSize << "\n";

    genomeSize = packGenome(genome, fileSize, referenceIDMap);
    std::cout << "        Reference blob size " << genomeSize << "\n";

    return genome;
}

/*
        INDEX_KMERS_THREAD()
        --------------------
        Dispatches on pattern.use64 to handle both <=32-mer (uint32_t packed)
        and >32-mer (uint64_t packed) cases.  Both hash down to uint32_t for
        the index, so outerBlob/innerBlob are unchanged.
*/
void index_kmers_thread(char* genome, uint64_t offset, uint64_t genomeSize,
                        std::vector<protected_vector<uint32_t>>& kmersMap, uint32_t MASK,
                        const StripePattern& pat, const StripePattern& pat2)
{
    genome += offset;
    auto start = std::chrono::steady_clock::now();
    uint64_t lastDisplayedPercent = (uint64_t)-10;

    const bool needsHashing = (MASK != UINT32_MAX);
    const bool useStripe = !pat.pattern_str.empty();
    const bool useDual = useStripe && !pat2.pattern_str.empty();

    if (useStripe) {
        if (!pat.use64) {
            // ---- 32-bit packed (weight <= 32) ----
            for (uint32_t pos = 0; pos + (uint32_t)pat.window <= genomeSize; ++pos) {
                uint32_t packed;
                if (pat.encode32(genome + pos, packed)) {
                    uint32_t kmerHash =
                        needsHashing ? murmurHash3(static_cast<uint64_t>(packed)) & MASK : packed;
                    kmersMap[kmerHash].push_back(pos + (uint32_t)offset);
                }

                // Second pattern — same position, independent hash
                if (useDual) {
                    uint32_t packed2;
                    if (pat2.encode32(genome + pos, packed2)) {
                        uint32_t kmerHash2 =
                            needsHashing ? murmurHash3(static_cast<uint64_t>(packed2)) & MASK
                                         : packed2;
                        kmersMap[kmerHash2].push_back(pos + (uint32_t)offset);
                    }
                }

                displayProgress(start, lastDisplayedPercent, pos, genomeSize, 10);
            }
        }
        else {
            // ---- 64-bit packed (weight 33..64) ----
            for (uint32_t pos = 0; pos + (uint32_t)pat.window <= genomeSize; ++pos) {
                uint64_t packed;
                if (pat.encode64(genome + pos, packed)) {
                    uint32_t kmerHash = murmurHash3(packed) & MASK;
                    kmersMap[kmerHash].push_back(pos + (uint32_t)offset);
                }

                if (useDual) {
                    uint64_t packed2;
                    if (pat2.encode64(genome + pos, packed2)) {
                        uint32_t kmerHash2 = murmurHash3(packed2) & MASK;
                        kmersMap[kmerHash2].push_back(pos + (uint32_t)offset);
                    }
                }

                displayProgress(start, lastDisplayedPercent, pos, genomeSize, 10);
            }
        }
    }
    else {
        // ---- Plain 32-mer (no stripe, legacy non-stripe mode) ----
        int ncountdown = 0;
        for (size_t pos = 0; pos < 31; pos++)
            if (toupper(genome[pos]) == 'N')
                ncountdown = pos + 1;

        uint64_t pkmer = encode_kmer_2bit::pack_32mer(genome);
        uint64_t remkp = encode_kmer_2bit::reverse_complement_32mer(pkmer);
        pkmer >>= 2;
        remkp <<= 2;
        char* encode_pos = genome + 31;

        for (uint32_t pos = 0; pos < genomeSize; pos++) {
            if (toupper(*encode_pos) == 'N')
                ncountdown = 32;
            uint64_t new_base = encode_kmer_2bit::pack_1mer(*encode_pos++);
            pkmer = (pkmer << 2) | new_base;
            uint32_t kmerHash = needsHashing ? murmurHash3(pkmer) & MASK
                                             : static_cast<uint32_t>(pkmer >> 32) & MASK;
            if (ncountdown == 0)
                kmersMap[kmerHash].push_back((uint32_t)(pos + offset));
            else
                ncountdown--;
            displayProgress(start, lastDisplayedPercent, pos, genomeSize, 10);
        }
    }
}

/*
        INDEX_KMERS()
        -------------
*/
void index_kmers(char* genome, uint64_t genomeSize,
                 std::vector<protected_vector<uint32_t>>& kmersMap, uint32_t MASK,
                 const StripePattern& pat, const StripePattern& pat2)
{
    size_t thread_count = std::thread::hardware_concurrency();
    uint64_t chunk_size = genomeSize / thread_count;
    uint64_t start = 0;

    std::vector<std::thread> threads;

    std::cout << "Launching " << thread_count << " threads each with " << chunk_size << " pieces\n";

    for (size_t i = 0; i < thread_count - 1; i++) {
        threads.push_back(std::thread(index_kmers_thread, genome, start, chunk_size,
                                      std::ref(kmersMap), MASK, std::cref(pat), std::cref(pat2)));
        start += chunk_size;
    }
    index_kmers_thread(genome, start, genomeSize - start - 31, kmersMap, MASK, pat, pat2);

    for (auto& thread : threads)
        thread.join();
}
