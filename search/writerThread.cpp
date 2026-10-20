//
//  writerThread.cpp
//  KISSfilter
// Producer-Cnsumer writeer thread
// KISS search workers queue matches for writing to file
//
//  Created by Shlomo Geva on 30/6/2024.
//

#include "writerThread.hpp"

#include <condition_variable>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <mutex>
#include <queue>
#include <thread>
#include <vector>

/* EXAMPLE of a worker thread */
//
// void worker(SafeQueue& queue, int threadId) {
//    // Simulate some processing
//    for (int i = 0; i < 1000; ++i) {
//        Result result;
//        result.data = "Thread " + std::to_string(threadId) + " result " +
//        std::to_string(i); queue.push(result);
//    }
//}

extern QueryTable ref_ID_Table;

static std::string getFirstToken(const std::string& str)
{
    std::istringstream iss(str);
    std::string firstToken;
    iss >> firstToken;
    return firstToken;
}

void writer(SafeQueue& queue, const std::string& filename)
{
    std::ofstream outFile(filename);
    Result result;
    std::string refID = "";
    uint32_t key = 0;
    char buffer[1024];
    uint32_t resultsCount = 0;
    while (queue.pop(result)) {
        // here format the output result and write it out

        // key is reference# in fasta file
        // note: the queryTable struct is defined in serialiseKmersMap (for some
        // unknown reason...
        // :-) )

        key = ref_ID_Table.query(result.refBegin);
        if (key == -1) {
            std::cout << "referenceIndex = " << result.refBegin << " Ref_ID not found."
                      << std::endl;
        }
        else {
            refID = getFirstToken(ref_ID_Table.getValue(key));
        }

        snprintf(buffer, 1024, "%c %8u\t%s\t%6u\t%5.1f\t%s\n", (result.revComp == 0 ? ' ' : '*'),
                 result.sampleIndex + 1, refID.c_str(), result.refBegin - key + 1, result.score,
                 result.cigar.c_str());
        outFile << std::string(buffer);
        ++resultsCount;
    }
    std::cout << resultsCount << " matches written\n";
}
