//
//  SafeQueue.h
//  KISSfilter
//
//  Created by Shlomo Geva on 10/7/2024.
//

#ifndef SafeQueue_h
#define SafeQueue_h
#include <condition_variable>
#include <mutex>
#include <queue>

struct Result {
    // Add your data members here
    uint32_t sampleIndex;
    uint32_t refBegin;
    float score;
    std::string cigar;
    uint32_t revComp;
};

template <typename T> class SafeQueue {
public:
    void enqueue(const T& item)
    {
        std::lock_guard<std::mutex> lock(mutex);
        queue.push(item);
        condition.notify_one(); // Notify one waiting thread
    }

    bool dequeue(T& item)
    {
        std::unique_lock<std::mutex> lock(mutex);
        condition.wait(lock, [this] { return !queue.empty() || finished; });
        if (!queue.empty()) {
            item = queue.front();
            queue.pop();
            return true;
        }
        return false;
    }

    void setFinished()
    {
        std::lock_guard<std::mutex> lock(mutex);
        finished = true;
        condition.notify_all(); // Notify all waiting threads
    }

    bool isFinished() const
    {
        std::lock_guard<std::mutex> lock(mutex);
        return finished && queue.empty();
    }

private:
    std::queue<T> queue;
    mutable std::mutex mutex;
    std::condition_variable condition;
    bool finished = false;
};

#endif /* SafeQueue_h */
