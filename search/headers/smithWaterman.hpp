//
//  smithWaterman.hpp
//  indexKISS
//
//  Created by Shlomo Geva on 15/7/2023.
//

#ifndef smithWaterman_hpp
#define smithWaterman_hpp

#include <stdio.h>
#include <string.h>

#include <algorithm>
#include <climits>
#include <iostream>
#include <map>
#include <span>
#include <string>
#include <tuple>
#include <unordered_map>
#include <vector>

#include "SafeQueue.hpp"
#include "main.hpp" // includes ReadInfo struct
#include "threadPool.hpp"

// Utility: convert string to uppercase
std::string toUpperCase(const std::string& str);

// Reverse complement a DNA sequence (handles mixed case and non-ACGT)
std::string reverseComplement(const std::string& sequence);

// Perform Smith-Waterman alignment using StripedSmithWaterman library
void alignSequences(std::string_view reference, std::string_view query, std::string& cigar,
                    int& sw_score, uint32_t& ref_begin, uint32_t& ref_end);

// Updated getSW with zero-copy blob + ReadInfo
bool getSW(const char* blob, const std::vector<ReadInfo>& readInfos, size_t startIndex,
           const std::span<char>& genome, SafeQueue<Result>& queue,
           const std::vector<std::tuple<uint32_t, uint32_t, uint32_t, uint32_t, uint32_t>>&
               candidate_results,
           int matchScore, int mismatchScore, int gapScore, int cutoff);

uint64_t reverse_complement(uint64_t kmer, uint64_t bases);

struct SortByScore {
    bool operator()(
        const std::tuple<uint32_t, uint32_t, uint32_t, double, int, std::string, int>& a,
        const std::tuple<uint32_t, uint32_t, uint32_t, double, int, std::string, int>& b) const
    {
        if (std::get<3>(a) != std::get<3>(b)) {
            return std::get<3>(a) > std::get<3>(b); // Sort by swScore
        }
        else if (std::get<4>(a) != std::get<4>(b)) {
            return std::get<4>(a) < std::get<4>(b); // Sort by swLength
        }
        else if (std::get<0>(a) != std::get<0>(b)) {
            return std::get<0>(a) < std::get<0>(b); // Sort by read number
        }
        else {
            return std::get<1>(a) < std::get<1>(b); // Sort by reference number
        }
    }
};

#endif /* smithWaterman_hpp */
