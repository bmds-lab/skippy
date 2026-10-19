//
//  findSpans.hpp
//  KISSfilter
//
//  Created by Shlomo Geva on 3/7/2024.
//

#ifndef findSpans_hpp
#define findSpans_hpp

#include <stdio.h>

#include <algorithm>
#include <cstdint>
#include <iostream>
#include <set>
#include <span>
#include <tuple>
#include <vector>
#include <cstddef>
// #include "main.hpp"

struct ListStruct
{
    uint32_t start;
    uint32_t length;
};

struct ListStructPtr
{
    const uint32_t *start; // pointer to the start of the sub-array in innerMapBlob
    uint32_t length;       // length of the sub-array
};

std::vector<std::tuple<uint32_t, uint32_t>>
findSpansPtr(std::vector<ListStructPtr> &seedLists,
             std::vector<uint32_t> &pos,
             uint32_t minMatches,
             uint32_t range,
             uint32_t slack = 0);

std::vector<std::tuple<uint32_t, uint32_t, uint32_t>>
findSpansPtrWithSpan(std::vector<ListStructPtr> &seedLists,
                     std::vector<uint32_t> &pos,
                     uint32_t minMatches,
                     uint32_t range,
                     uint32_t slack = 0);

//template <typename T> class SafeQueue;

#endif /* findSpans_hpp */

