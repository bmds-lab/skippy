//
//  kmerUtilities.hpp
//  indexKISS
//
//  Created by Shlomo Geva on 13/7/2023.
//
#pragma once

#include <stdio.h>
#include <sys/stat.h>

#include <algorithm>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <set>
#include <string>
#include <vector>

struct Span {
    uint32_t start;
    uint32_t end;
};

struct store_position {
    const uint32_t* position;
    std::size_t size;
};

uint64_t packKmer(const char* sequence);
std::string unpackKmer(uint64_t packed_sequence);
char* read_entire_file(const char* filename, uint64_t& fileSize);
