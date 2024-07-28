#include "tree_barrier_impl.h"
#include <iostream>
#include <thread>
#include <vector>

//const int N_THREADS = 16; // Number of threads to test with
std::atomic<int> threadLoopCount(0);
std::atomic<int> globalThreadCount(0);
std::atomic<int> globalPrintingPing(0);
std::mutex LOOPING_MUTEX;


void threadFunction(BinaryTreeBarrier &barrier, int threadID) {
    thread_id = threadID;
    
    for (size_t i = 0; i < threadLoopCount; i++)
    {
        barrier.barrier();
        PRINTING_MUTEX.lock();
            std::cout << "threadFunction(): Thread "<< threadID <<" has exited the barrier().\n";
            globalPrintingPing++;
            if(globalPrintingPing == globalThreadCount)
            {
                std::cout << "threadFunction(): Thread " << threadID << " was the last one to leave barrier " << i << "." << endl;
                std::cout << "=========================================================================================================" << endl;
                globalPrintingPing = 0;
            }
        PRINTING_MUTEX.unlock();
    }
}

int main(int argc, char** argv) {
    int threadCount = atoi(argv[1]);
    globalThreadCount = threadCount;
    threadLoopCount = atoi(argv[2]);
    PRINTING_MUTEX.lock();
        std::cout << "main(): Thread count is: "<< threadCount << "." << endl;
    PRINTING_MUTEX.unlock();
    BinaryTreeBarrier barrier(threadCount); // Initialize barrier with N_THREADS

    std::vector<std::thread> threads;
    for (int threadID = 0; threadID < threadCount; ++threadID) {
        threads.emplace_back(threadFunction, std::ref(barrier), threadID);
    }

    // Gather all of the threads before terminating
    for (auto &t : threads) {
        
        t.join(); // Wait for all threads to complete
    }

    std::cout << "main(): Done running and joined the threads.\n";

    return 0;
}

