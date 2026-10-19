//
//  kmerUtilities.cpp
//  indexKISS
//
//  Created by Shlomo Geva on 13/7/2023.
//

#include "kmerUtilities.hpp"
#include <algorithm> // for std::min and std::max
#include <climits>
#include <iostream>
#include <limits>
#include <utility>
#include <vector>

//extern uint32_t KMERSIZE;
extern uint32_t NUM_SEEDS; // seed skipping distancee
extern int CUTOFF;
extern uint64_t kmer_encoding_table[];

uint32_t getSpan(const std::vector<uint32_t> &pos, uint32_t span)
{
    // calculate the span of matches covered by seeds
    if (pos.empty())
    {
        return 0;
    }

    std::vector<Span> spans;

    for (const uint32_t p : pos)
    {
        spans.push_back({p, p + span - 1});
    }

    std::vector<Span> mergedSpans;
    std::sort(spans.begin(), spans.end(),
              [](const Span &a, const Span &b) { return a.start < b.start; });

    for (const Span &span : spans)
    {
        if (mergedSpans.empty() || span.start > mergedSpans.back().end)
        {
            mergedSpans.push_back(span);
        }
        else
        {
            mergedSpans.back().end = std::max(mergedSpans.back().end, span.end);
        }
    }

    uint32_t totalSpan = 0;
    for (const Span &span : mergedSpans)
    {
        totalSpan += (span.end - span.start + 1);
    }

    return totalSpan;
}

/*
    PACK KMER
    -----------
*/
/*
uint64_t packKmer(const char *sequence)
{
    uint64_t packed = 0;

    for (int pos = 0; pos < KMERSIZE; pos++)
        packed = (packed << 2) | kmer_encoding_table[(size_t)sequence[pos]];

    return packed;
}
*/
/*
    UNPACK KMER
    -------------
*/
/*
std::string unpackKmer(uint64_t packed_sequence)
{
    std::string sequence;
    for (int32_t pos = KMERSIZE - 1; pos >= 0; pos--)
    {
        switch ((packed_sequence >> (pos * 2)) & 3)
        {
        case 0:
            sequence += 'A';
            break;
        case 1:
            sequence += 'C';
            break;
        case 2:
            sequence += 'G';
            break;
        case 3:
            sequence += 'T';
            break;
        }
    }
    return sequence;
}
*/
/*
    READ_ENTIRE_FILE
 */
char *read_entire_file(const char *filename, uint64_t &fileSize)
{
    FILE *fp;
    struct stat details;
    char *contents = NULL;
    fileSize = 0;
    if ((fp = fopen(filename, "rb")) != NULL)
    {
        if (fstat(fileno(fp), &details) == 0)
        {
            if (details.st_size != 0 || details.st_size > UINT32_MAX)
            {
                contents = (char *)malloc(details.st_size);
                if (fread(contents, details.st_size, 1, fp) != 1)
                {
                    free(contents);
                    contents = NULL;
                }
                else
                {
                    fileSize = details.st_size;
                }
            }
        }
        fclose(fp);
    }
    return contents;
}

void displayProgress(uint64_t current, uint64_t total, int desiredUpdateInterval)
{
    uint64_t percent = (current * 100) / total;
    static uint64_t lastDisplayedPercent =
        -desiredUpdateInterval; // Initialize to a value that will trigger the first update
    if (percent - lastDisplayedPercent >= desiredUpdateInterval)
    {
        std::cout << "Progress: " << std::setw(3) << percent << "%" << std::endl;
        std::cout.flush();
        lastDisplayedPercent = percent;
    }
}
