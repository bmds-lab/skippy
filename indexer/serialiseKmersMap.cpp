/*
        SERIALISEKMERSMAP.CPP
        ---------------------
        indexReference

        Created by Shlomo Geva on 19/7/2023.
        Parallel serialisation added Feb 2026.
        33-bit hi-bit bitmap added Feb 2026.

        Phase 1 - Parallel build (all hardware threads):
            Each thread owns a contiguous range of kmersMap slots.
            Sorts each slot, copies positions into local innerChunk,
            records zero-based local offsets in local outerChunk,
            frees each slot immediately. No synchronisation.

        Phase 2 - Parallel offset fixup + hi-bit bitmap:
            Compute uint64_t prefix sums across chunks.
            Each thread independently fixes up its own outerChunk and
            builds its own local bitmap slice. No shared writes, no
            contention. Local bitmaps merged by OR into global hiBitmap
            after all threads complete.

        Phase 3 - Serial write:
            Walk chunks in order writing outerChunk and innerChunk.
            Then write hiBitmap as a separate file.
*/
#include "serialiseKmersMap.hpp"

#include <algorithm>
#include <cstring>
#include <fstream>
#include <iostream>
#include <limits>
#include <thread>
#include <vector>

bool serializeMap(std::vector<protected_vector<uint32_t>>& kmersMap,
                  const std::string& innerMapFilename, const std::string& outerMapFilename)
{
    constexpr uint32_t M = 10000;

    const size_t totalSlots = kmersMap.size();
    const size_t numThreads = std::max(size_t(1), (size_t)std::thread::hardware_concurrency());
    const size_t slotsPerThread = (totalSlots + numThreads - 1) / numThreads;

    std::cout << "Serialising " << totalSlots << " slots using " << numThreads << " threads...\n";

    struct ChunkResult {
        std::vector<uint32_t> outerChunk;
        std::vector<uint32_t> innerChunk;
        uint32_t longest = 0;
        uint32_t longCount = 0;
    };

    std::vector<ChunkResult> chunks(numThreads);

    // ------------------------------------------------------------------
    // Phase 1: parallel build
    // ------------------------------------------------------------------
    {
        std::vector<std::thread> threads;
        threads.reserve(numThreads);

        for (size_t t = 0; t < numThreads; ++t) {
            const size_t slotStart = t * slotsPerThread;
            const size_t slotEnd = std::min(slotStart + slotsPerThread, totalSlots);
            if (slotStart >= totalSlots)
                break;

            threads.emplace_back([&, t, slotStart, slotEnd]() {
                ChunkResult& cr = chunks[t];
                const size_t slotCount = slotEnd - slotStart;
                cr.outerChunk.reserve(slotCount);

                uint32_t localOffset = 0;

                for (size_t i = slotStart; i < slotEnd; ++i) {
                    auto& vec = kmersMap[i];
                    const uint32_t len = static_cast<uint32_t>(vec.size());

                    if (len > M) {
                        ++cr.longCount;
                        if (len > cr.longest)
                            cr.longest = len;
                    }

                    cr.outerChunk.push_back(localOffset);

                    if (len > 1)
                        std::sort(vec.begin(), vec.end());

                    cr.innerChunk.insert(cr.innerChunk.end(), vec.begin(), vec.end());

                    localOffset += len;

                    {
                        protected_vector<uint32_t> tmp;
                        vec.swap(tmp);
                    }
                }
            });
        }

        for (auto& th : threads)
            th.join();
    }

    std::cout << "Phase 1 complete. Fixing up offsets...\n";

    // ------------------------------------------------------------------
    // Phase 2: parallel offset fixup + hi-bit bitmap
    // Each thread fixes up its own outerChunk and builds its own local
    // bitmap slice. No shared writes, no contention.
    // ------------------------------------------------------------------

    std::vector<uint64_t> prefixSum(numThreads, 0);
    for (size_t t = 1; t < numThreads; ++t)
        prefixSum[t] = prefixSum[t - 1] + chunks[t - 1].innerChunk.size();

    const uint64_t totalEntries =
        prefixSum[numThreads - 1] + chunks[numThreads - 1].innerChunk.size();
    if (totalEntries > 0xFFFFFFFFULL)
        std::cout << "Note: innerBlob has " << totalEntries
                  << " entries (>4G) - hi-bit bitmap active\n";

    const size_t bitmapWords = (totalSlots + 31) / 32;
    std::vector<uint32_t> hiBitmap(bitmapWords, 0);
    std::vector<std::vector<uint32_t>> localBitmaps(numThreads);

    {
        std::vector<std::thread> threads;
        threads.reserve(numThreads);

        for (size_t t = 0; t < numThreads; ++t) {
            threads.emplace_back([&, t]() {
                const uint64_t base = prefixSum[t];
                const size_t slotStart = t * slotsPerThread;
                const size_t slotEnd = std::min(slotStart + slotsPerThread, totalSlots);
                const size_t slotCount = slotEnd - slotStart;

                const size_t localWords = (slotCount + 31) / 32;
                localBitmaps[t].assign(localWords, 0);

                size_t localSlot = 0;
                for (uint32_t& entry : chunks[t].outerChunk) {
                    const uint64_t absOffset = base + entry;
                    entry = static_cast<uint32_t>(absOffset & 0xFFFFFFFF);
                    if (absOffset >> 32)
                        localBitmaps[t][localSlot / 32] |= (1u << (localSlot & 31));
                    ++localSlot;
                }
            });
        }

        for (auto& th : threads)
            th.join();
    }

    // Merge local bitmaps into global hiBitmap
    for (size_t t = 0; t < numThreads; ++t) {
        const size_t slotStart = t * slotsPerThread;
        for (size_t i = 0; i < localBitmaps[t].size(); ++i)
            if (localBitmaps[t][i])
                hiBitmap[(slotStart / 32) + i] |= localBitmaps[t][i];
        {
            std::vector<uint32_t> tmp;
            localBitmaps[t].swap(tmp);
        }
    }

    std::cout << "Phase 2 complete. Writing to disk...\n";

    // ------------------------------------------------------------------
    // Phase 3: serial write
    // ------------------------------------------------------------------
    {
        std::ofstream innerFile(innerMapFilename, std::ios::binary);
        std::ofstream outerFile(outerMapFilename, std::ios::binary);
        if (!innerFile || !outerFile) {
            std::cerr << "Error: could not open " << innerMapFilename << " or " << outerMapFilename
                      << " for writing\n";
            return false;
        }

        for (size_t t = 0; t < numThreads; ++t) {
            ChunkResult& cr = chunks[t];

            if (!cr.outerChunk.empty())
                outerFile.write(reinterpret_cast<const char*>(cr.outerChunk.data()),
                                cr.outerChunk.size() * sizeof(uint32_t));

            if (!cr.innerChunk.empty())
                innerFile.write(reinterpret_cast<const char*>(cr.innerChunk.data()),
                                cr.innerChunk.size() * sizeof(uint32_t));

            {
                std::vector<uint32_t> tmp;
                cr.outerChunk.swap(tmp);
            }
            {
                std::vector<uint32_t> tmp;
                cr.innerChunk.swap(tmp);
            }
        }

        innerFile.close();
        outerFile.close();
        if (!innerFile || !outerFile) {
            std::cerr << "Error: could not write " << innerMapFilename << " or " << outerMapFilename
                      << '\n';
            return false;
        }
    }

    // Write hi-bit bitmap
    {
        std::string hiBitFilename = outerMapFilename;
        size_t pos = hiBitFilename.rfind("_OuterBlob.idx");
        if (pos != std::string::npos)
            hiBitFilename.replace(pos, 14, "_OuterBlobHiBit.idx");
        else
            hiBitFilename += "_HiBit";

        std::ofstream hiBitFile(hiBitFilename, std::ios::binary);
        if (!hiBitFile) {
            std::cerr << "Error: could not open " << hiBitFilename << " for writing\n";
            return false;
        }
        hiBitFile.write(reinterpret_cast<const char*>(hiBitmap.data()),
                        hiBitmap.size() * sizeof(uint32_t));
        hiBitFile.close();
        if (!hiBitFile) {
            std::cerr << "Error: could not write " << hiBitFilename << '\n';
            return false;
        }
        std::cout << "Hi-bit bitmap written (" << hiBitmap.size() * 4 / (1024 * 1024)
                  << " MB)\n";
    }

    // Reduce stats
    uint32_t longest = 0, longCount = 0;
    for (const auto& cr : chunks) {
        if (cr.longest > longest)
            longest = cr.longest;
        longCount += cr.longCount;
    }
    std::cout << "Longest inner vector: " << longest << "\n";
    std::cout << "Inner vectors length > " << M << ": " << longCount << "\n";
    return true;
}

void deserializeMap(const std::string& innerMapFilename, const std::string& outerMapFilename,
                    std::vector<uint32_t>& innerMapBlob, std::vector<uint32_t>& outerMapBlob,
                    std::vector<uint32_t>& hiBitmap)
{
    std::ifstream innerMapFile(innerMapFilename, std::ios::binary);
    std::ifstream outerMapFile(outerMapFilename, std::ios::binary);

    innerMapFile.seekg(0, std::ios::end);
    size_t innerBlobSize = innerMapFile.tellg();
    innerMapFile.seekg(0, std::ios::beg);

    outerMapFile.seekg(0, std::ios::end);
    size_t outerBlobSize = outerMapFile.tellg();
    outerMapFile.seekg(0, std::ios::beg);

    innerMapBlob.resize(innerBlobSize / sizeof(uint32_t));
    innerMapFile.read(reinterpret_cast<char*>(innerMapBlob.data()), innerBlobSize);

    outerMapBlob.resize(outerBlobSize / sizeof(uint32_t));
    outerMapFile.read(reinterpret_cast<char*>(outerMapBlob.data()), outerBlobSize);

    innerMapFile.close();
    outerMapFile.close();

    // Load hi-bit bitmap if present
    std::string hiBitFilename = outerMapFilename;
    size_t pos = hiBitFilename.rfind("_OuterBlob.idx");
    if (pos != std::string::npos)
        hiBitFilename.replace(pos, 14, "_OuterBlobHiBit.idx");
    else
        hiBitFilename += "_HiBit";

    std::ifstream hiBitFile(hiBitFilename, std::ios::binary);
    if (hiBitFile.is_open()) {
        hiBitFile.seekg(0, std::ios::end);
        size_t hiBitSize = hiBitFile.tellg();
        hiBitFile.seekg(0, std::ios::beg);
        hiBitmap.resize(hiBitSize / sizeof(uint32_t));
        hiBitFile.read(reinterpret_cast<char*>(hiBitmap.data()), hiBitSize);
        hiBitFile.close();
        std::cerr << "Loaded hi-bit bitmap (" << hiBitSize / (1024 * 1024) << " MB)\n";
    }
    else
        hiBitmap.clear();
}

std::vector<uint32_t> getInnerVector(const std::vector<uint32_t>& innerMapBlob,
                                     const std::vector<uint32_t>& outerMapBlob, size_t index)
{
    std::vector<uint32_t> innerVector;
    if (index < outerMapBlob.size()) {
        size_t startOffset = outerMapBlob[index];
        size_t endOffset =
            (index + 1 < outerMapBlob.size()) ? outerMapBlob[index + 1] : innerMapBlob.size();
        innerVector = std::vector<uint32_t>(innerMapBlob.begin() + startOffset,
                                            innerMapBlob.begin() + endOffset);
    }
    return innerVector;
}

bool writeTextBlobToFile(const char* text, std::size_t length, const std::string& filename)
{
    std::ofstream outputFile(filename, std::ios::binary);
    if (!outputFile.is_open())
        return false;
    outputFile.write(text, length);
    outputFile.close();
    return static_cast<bool>(outputFile);
}

std::pair<char*, std::size_t> readTextBlobFromFile(const std::string& filename)
{
    std::ifstream inputFile(filename, std::ios::binary | std::ios::ate);
    if (!inputFile.is_open())
        return {nullptr, 0};

    std::size_t length = inputFile.tellg();
    inputFile.seekg(0);

    std::vector<char> buffer(length);
    if (!inputFile.read(buffer.data(), length))
        return {nullptr, 0};

    inputFile.close();
    char* textBlob = new char[length];
    std::memcpy(textBlob, buffer.data(), length);
    return {textBlob, length};
}
