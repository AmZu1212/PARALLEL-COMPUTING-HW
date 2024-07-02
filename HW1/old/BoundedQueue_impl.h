#include "BoundedQueue.h"
#include <mutex>
#include <condition_variable>
#include <deque>
#include <thread>
#include <iostream>
class BoundedQueue :  public BoundedQueueAbstract {

    private:
    int capacity;
    int count;
    std::deque<int> BQueue;
    std::mutex QMutex;
    std::condition_variable notEmpty; // for pop
    std::condition_variable notFull; // for push

    public:
    // Constructor to initialize the capacity
    BoundedQueue(int cap) : capacity(cap), count(0) {}

    // Return the number of elements in the queue
    int size() override {
        std::unique_lock<std::mutex> lock(QMutex);
        return count;
    }

    // Pop the next element (integer value) from the queue.
    // if the buffer is empty, the calling thread waits until being notified of new elements in the queue
    int pop() override {
        std::unique_lock<std::mutex> lock(QMutex);
        notEmpty.wait(lock, [this]() {return count > 0;});
        int value = BQueue.front();
        BQueue.pop_front();
        --count;
        notFull.notify_one();
        return value;
    }

    // Push a new integer to the queue.
    // if the buffer is full, the calling thread should wait until being notified that the queue is not full anymore
    // v - the new integer to push into the queue.
    void push(int v) override {
        std::unique_lock<std::mutex> lock(QMutex);
        notFull.wait(lock, [this]() {return count < capacity;});
        BQueue.push_back(v);
        ++count;
        notEmpty.notify_one();
    }
};
