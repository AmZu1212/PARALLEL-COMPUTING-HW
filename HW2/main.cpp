#include "tree_barrier_impl.h"
#include <iostream>
#include <thread>
#include <vector>

const int N_THREADS = 8; // Number of threads to test with

void threadFunction(BinaryTreeBarrier &barrier, int thread_id) {
    std::cout << "Thread " << thread_id << " is waiting at the barrier.\n";
    barrier.barrier();
    std::cout << "Thread " << thread_id << " has passed the barrier.\n";
}

int main() {
    BinaryTreeBarrier barrier(N_THREADS); // Initialize barrier with N_THREADS

    std::vector<std::thread> threads;
    for (int i = 0; i < N_THREADS; ++i) {
        threads.emplace_back(threadFunction, std::ref(barrier), i);
    }

    for (auto &t : threads) {
        t.join(); // Wait for all threads to complete
    }

    std::cout << "All threads have passed the barrier.\n";

    return 0;
}
