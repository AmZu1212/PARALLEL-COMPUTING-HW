#include "BoundedQueue1p1c.h"
#include <atomic>
#include <thread>
#include <iostream>
#include <chrono>
#include <vector>

namespace chrono = std::chrono;

// trying again the rotating array.
// i forgot to use atomics damn
class BoundedQueue1p1c : public BoundedQueueAbstract_1p1c
{
    private:
        // changed everything to atomics, buffer has to be atomic
        std::vector<std::atomic<int>> buffer;
        int head;
        int tail;
        std::atomic<int> count;
        int capacity;
        // not using memory flags.
        // load is read, store is write.

        // i think head and tail dont need to be atomic because
        // there is only 1 producer and 1 consumer, and 
        // each one updates the different index.
    public:
        BoundedQueue1p1c(int capacity)
            : buffer(capacity), head(0), tail(0), count(0), capacity(capacity) {}

        int size() override
        {
            return count.load();
        }

        bool pop(int &val) override
        {
            if (count.load() == 0)
            {
                return false; // Queue is empty
            }

            val = buffer[head].load();
            head = (head + 1) % capacity; // update head mod n
            count.fetch_sub(1);// this is -- in atomics

            return true;
        }

        bool push(int v) override
        {
            if (count.load() == capacity)
            {
                return false; // Queue is full
            }

            buffer[tail].store(v);
            tail = (tail + 1) % capacity; // update tail mod n
            count.fetch_add(1);// this is ++ in atomics

            return true;
        }
};
