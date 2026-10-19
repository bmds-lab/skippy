//
//  main.hpp
//  indexKISS
//
//  Created by Shlomo Geva on 13/7/2023.
//

#ifndef main_h
#define main_h

#include "dirent.h"
#include <algorithm>
#include <array>
#include <chrono>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <map>
#include <mutex>
#include <random>
#include <regex>
#include <set>
#include <sstream>
#include <string>
#include <sys/stat.h>
#include <thread>
#include <unordered_map>
#include <vector>

#include "findSpans.hpp"
#include "serialiseKmersMap.hpp"

// New: ReadInfo struct for offset-based zero-copy reads
struct ReadInfo {
    size_t start;   // byte offset in blob where sequence begins
    size_t length;  // number of bases
};

struct SpanStruct
{
    uint32_t readID;
    uint32_t start;
    uint32_t end;
    uint32_t nSeed;
    uint32_t reverse;
};

#endif /* main_h */
