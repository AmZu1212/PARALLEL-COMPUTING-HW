#include "tree_barrier.h"
#include <vector>
#include <thread>
#include <atomic>
#include <condition_variable>
#include <mutex>
#include <cmath>
#include <iostream>

class BinaryTreeBarrier : public BinaryTreeBarrierAbstract {
public:
    BinaryTreeBarrier(int num_threads);
    void barrier() override;

private:
    struct BarrierNode {
        std::mutex mtx;
        std::condition_variable cv;
        int count = 0;
    };

    int num_threads;
    int levels;
    std::vector<BarrierNode> nodes;

    int calculate_levels(int num_threads);
    int get_node_index(int thread_id, int level);
    int get_thread_id();
};

BinaryTreeBarrier::BinaryTreeBarrier(int num_threads)
    : num_threads(num_threads), levels(calculate_levels(num_threads)), nodes((1 << (levels + 1)) - 1) {}

void BinaryTreeBarrier::barrier() {
    static thread_local int thread_id = std::hash<std::thread::id>()(std::this_thread::get_id()) % num_threads;
    int node_index = thread_id;

    // Traverse the binary tree from leaf to root
    for (int level = 0; level < levels; ++level) {
        node_index = get_node_index(thread_id, level);
        BarrierNode &node = nodes[node_index];

        std::unique_lock<std::mutex> lock(node.mtx);
        node.count++;

        // If we are the last thread to arrive, notify others and reset count
        if (node.count == 2) { // Adjust this based on how many threads you expect at each node
            node.count = 0; // Reset count for next use
            lock.unlock();
            node.cv.notify_all(); // Notify waiting threads
        } else {
            // Wait for the others to arrive
            node.cv.wait(lock);
        }
    }

    // Synchronize at the root node
    if (thread_id == 0) {
        BarrierNode &root = nodes[0];
        std::unique_lock<std::mutex> lock(root.mtx);
        while (root.count < num_threads) {
            root.cv.wait(lock);
        }
        root.count = 0; // Reset count for next use
        root.cv.notify_all(); // Notify all threads waiting at the root
    } else {
        BarrierNode &root = nodes[0];
        std::unique_lock<std::mutex> lock(root.mtx);
        if (++root.count == num_threads) {
            root.cv.notify_all(); // Notify all waiting threads if this is the last thread
        } else {
            root.cv.wait(lock); // Wait for others to reach the root
        }
    }
}

int BinaryTreeBarrier::calculate_levels(int num_threads) {
    return std::ceil(std::log2(num_threads));
}

int BinaryTreeBarrier::get_node_index(int thread_id, int level) {
    return (1 << level) - 1 + (thread_id >> (levels - level));
}