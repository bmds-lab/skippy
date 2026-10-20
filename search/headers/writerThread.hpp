//
//  writerThread.hpp
//  KISSfilter
//
//  Created by Shlomo Geva on 30/6/2024.
//

#ifndef writerThread_hpp
#define writerThread_hpp

#include <condition_variable>
#include <fstream>
#include <iostream>
#include <mutex>
#include <queue>
#include <thread>
#include <vector>

#include "serialiseKmersMap.hpp"

// Assuming Result is a structure to hold the result
struct Result {
    // Add your data members here
    uint32_t sampleIndex;
    uint32_t refBegin;
    float score;
    std::string cigar;
    uint32_t revComp;
};

// A thread-safe queue to hold results
class SafeQueue {
private:
    std::queue<Result> queue;
    std::mutex mtx;
    std::condition_variable cv;
    bool finished = false;

public:
    void push(const Result& result)
    {
        std::lock_guard<std::mutex> lock(mtx);
        queue.push(result);
        cv.notify_one();
    }

    bool pop(Result& result)
    {
        std::unique_lock<std::mutex> lock(mtx);
        cv.wait(lock, [this] { return !queue.empty() || finished; });

        if (!queue.empty()) {
            result = queue.front();
            queue.pop();
            return true;
        }

        return false;
    }

    void setFinished()
    {
        std::lock_guard<std::mutex> lock(mtx);
        finished = true;
        cv.notify_all();
    }
};

void writer(SafeQueue& queue, const std::string& filename);

#endif /* writerThread_hpp */
