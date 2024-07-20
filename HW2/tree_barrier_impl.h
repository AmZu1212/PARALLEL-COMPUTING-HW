#include "tree_barrier.h"
#include <vector>
#include <thread>
#include <atomic>
#include <condition_variable>
#include <mutex>
#include <cmath>
#include <iostream>

/*
// figure out how to run this
class BarrierNode
{
    public:
        std::atomic<int> count;
        std::condition_variable cv;
        std::mutex mtx;

        BarrierNode() : count(0) {}
};

class BinaryTreeBarrier : BinaryTreeBarrierAbstract
{
    private:
        std::vector<BarrierNode> nodes;
        int num_threads;
        int num_nodes;

    public:
        BinaryTreeBarrier(int num_threads) : num_threads(num_threads)
        {
            // The number of nodes in a binary tree is N-1. where N = number of threads.
            num_nodes = num_threads - 1;
            // Set the vector to the correct size
            nodes.resize(num_nodes);
        }

        void barrier() override
        {
            // get the entering thread's id, in order to know where to stick it.
            size_t thread_id = std::hash<std::thread::id>()(std::this_thread::get_id()) % num_threads;
            // calc the barrier id where the thread will wait.
            int node_idx = num_nodes + thread_id / 2;

            while (node_idx > 0)
            {
                // reference to the node?
                auto &node = nodes[node_idx];
                //lock on it
                std::unique_lock<std::mutex> lock(node.mtx);
                // increment 1/2
                node.count++;

                // Calculate the number of threads that need to synchronize at this level, each level doubles
                int threads_at_this_level = 2 << (int)log2(num_threads / (node_idx + 1));

                // If the required number of threads have not reached this node, wait
                if (node.count < threads_at_this_level) {
                    node.cv.wait(lock, [&node, threads_at_this_level] { return node.count == threads_at_this_level; });
                } else {
                    node.cv.notify_all(); // If everyone has arrived, notify all threads waiting on this node
                }

                lock.unlock();// free , node id update
                node_idx = (node_idx - 1) / 2;
            }

            // if a thread is here it reached the root
            auto &root = nodes[0];
            std::unique_lock<std::mutex> lock(root.mtx); // last lock waiting for everyone else...
            root.count++;

            if (root.count < num_threads)// wait for the others
            {
                root.cv.wait(lock, [&root, this]
                            { return root.count == num_threads; });
            }
            else
            {
                root.cv.notify_all();// everyone arrive, start to reset the tree
            }

            // resetting the barriers
            while (node_idx < num_nodes)
            {
                auto &node = nodes[node_idx];
                std::unique_lock<std::mutex> lock(node.mtx);// lock to reset 
                node.count = 0;
                node.cv.notify_all();// is this needed even? not sure

                lock.unlock(); // unlock and move ot the next one 
                node_idx = node_idx * 2 + 1;
            }
        }
};
*/

// Define a maximum number of nodes for the binary tree
constexpr size_t MAX_NODES = 128; // Adjust this number based on your maximum thread count

class BarrierNode
{
    public:
        std::atomic<int> count;
        std::condition_variable cv;
        std::mutex mtx;

        BarrierNode() : count(0) {}
};

class BinaryTreeBarrier : BinaryTreeBarrierAbstract
{
    private:
        std::array<BarrierNode, MAX_NODES> nodes;
        int num_threads;
        int num_nodes;

    public:
        BinaryTreeBarrier(int num_threads) : num_threads(num_threads)
        {
            // The number of nodes in a binary tree is N-1, where N = number of threads.
            num_nodes = num_threads - 1;
            // Ensure num_nodes does not exceed MAX_NODES
            if (num_nodes > MAX_NODES)
                throw std::runtime_error("Number of nodes exceeds maximum limit.");
        }

        void barrier() override
        {
            size_t thread_id = std::hash<std::thread::id>()(std::this_thread::get_id()) % num_threads;
            int node_idx = num_nodes + thread_id / 2;

            std::cout << "Thread " << thread_id << " is waiting at node " << node_idx << std::endl;

            while (node_idx > 0)
            {
                auto &node = nodes[node_idx];
                std::unique_lock<std::mutex> lock(node.mtx);
                node.count++;

                int threads_at_this_level = 2 << (int)log2(num_threads / (node_idx + 1));

                if (node.count < threads_at_this_level) {
                    std::cout << "Thread " << thread_id << " waiting at node " << node_idx << std::endl;
                    node.cv.wait(lock, [&node, threads_at_this_level] { return node.count == threads_at_this_level; });
                    std::cout << "Thread " << thread_id << " done waiting at node " << node_idx << std::endl;
                } else {
                    std::cout << "Thread " << thread_id << " notifying at node " << node_idx << std::endl;
                    node.cv.notify_all();
                }

                lock.unlock();
                node_idx = (node_idx - 1) / 2;
            }

            auto &root = nodes[0];
            std::unique_lock<std::mutex> lock(root.mtx);
            root.count++;

            if (root.count < num_threads) {
                std::cout << "Thread " << thread_id << " waiting at root" << std::endl;
                root.cv.wait(lock, [this, &root] { return root.count == num_threads; });
                std::cout << "Thread " << thread_id << " done waiting at root" << std::endl;
            } else {
                std::cout << "Thread " << thread_id << " notifying root" << std::endl;
                root.cv.notify_all();
            }

            // Reset the barrier
            while (node_idx < num_nodes)
            {
                auto &node = nodes[node_idx];
                std::unique_lock<std::mutex> lock(node.mtx);
                node.count = 0;
                node.cv.notify_all();

                lock.unlock();
                node_idx = node_idx * 2 + 1;
            }
        }

};
