#include "tree_barrier_impl.h"
#include <iostream>
#include <thread>
#include <vector>

//const int N_THREADS = 16; // Number of threads to test with

void threadFunction(BinaryTreeBarrier &barrier, int threadID) {
    thread_id = threadID;
    //std::cout << "main(): A thread is about to enter the barrier(). It's ID is: "<< threadID <<"\n";
    barrier.barrier();
    //PRINTING_MUTEX.lock();
    //    std::cout << "main(): Thread "<< threadID <<" has exited the barrier().\n";
    //PRINTING_MUTEX.unlock();
    //std::cout << "main(): A thread has exited the barrier().\n";
}

int main(int argc, char** argv) {
    int threadCount = atoi(argv[1]);
    //PRINTING_MUTEX.lock();
    //    std::cout << "main(): Thread count is: "<< threadCount << "." << endl;
    //PRINTING_MUTEX.unlock();
    BinaryTreeBarrier barrier(threadCount); // Initialize barrier with N_THREADS

    std::vector<std::thread> threads;
    for (int threadID = 0; threadID < threadCount; ++threadID) {
        threads.emplace_back(threadFunction, std::ref(barrier), threadID);
    }

    // Gather all of the threads before terminating
    for (auto &t : threads) {
        
        t.join(); // Wait for all threads to complete
    }

    //std::cout << "main(): Done running and joined the threads.\n";

    return 0;
}

