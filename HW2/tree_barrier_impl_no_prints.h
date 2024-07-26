#include "tree_barrier.h"
#include <vector>
#include <thread>
#include <atomic>
#include <condition_variable>
#include <mutex>
#include <cmath>
#include <iostream>
#define DEBUG 1
using namespace std;

// Define a maximum number of nodes for the binary tree
constexpr size_t MAX_NODES = 64;//128;//256;//512;//1024;//16384; // Adjust this number based on your maximum thread count 
                                    //(future note: 16384 is the maximum limit before crashing.)
                                    //(Also: test only runs 64 threads.)

int next_id = 0;
std::mutex PRINTING_MUTEX;
std::mutex id_mutex;
thread_local int currentHeight = 1;

class BarrierNode
{
    public:
        std::atomic<int> count; // if this works, we can implement a phase lock
        std::condition_variable cv;
        std::mutex mtx;

        BarrierNode() : count(0) {}
};

class BinaryTreeBarrier : BinaryTreeBarrierAbstract
{
    private:
        std::array<BarrierNode, MAX_NODES> nodes;
        int numOfThreads;
        int numOfNodes;
        int treeHeight;

    public:
        int get_thread_id() {
            std::lock_guard<std::mutex> lock(id_mutex);
            return next_id++; // Assign the current ID and increment for the next thread
        }


        BinaryTreeBarrier(int numOfThreads) : numOfThreads(numOfThreads)
        {
            treeHeight = (int)log2(numOfThreads - 1);
            
            // The number of nodes in a binary tree is N-1, where N = number of threads.
            numOfNodes = numOfThreads - 1;

            // Ensure num_nodes does not exceed MAX_NODES. REVISE THIS LATER
            if (numOfNodes > MAX_NODES)
                throw std::runtime_error("Number of nodes exceeds maximum limit.");
        }

        void barrier() override
        {
            int threadID = thread_id;
            int nodeID = ((numOfNodes + threadID - 1) / 2);

            while (nodeID > 0)
            {
                auto &node = nodes[nodeID];
                std::unique_lock<std::mutex> lock(node.mtx);
                node.count++;
                int threads_at_this_level = (int)pow(2, currentHeight);
                if (node.count < threads_at_this_level) {
                    node.cv.wait(lock, [&node, threads_at_this_level] { return node.count == threads_at_this_level; });
                } else {
                    node.cv.notify_all();
                }

                lock.unlock();
                currentHeight++;
                nodeID = (nodeID - 1) / 2;
            }

            auto &root = nodes[0];
            std::unique_lock<std::mutex> lock(root.mtx);
            root.count++;
            
            if (root.count < numOfThreads) {
                root.cv.wait(lock, [this, &root] { return root.count == numOfThreads; });
            } else {
                root.cv.notify_all();
            }
            /*
            // Reset the barrier (add later)
            while (nodeID < numOfNodes)
            {
                auto &node = nodes[nodeID];
                std::unique_lock<std::mutex> lock(node.mtx);
                node.count = 0;
                node.cv.notify_all();

                lock.unlock();
                nodeID = nodeID * 2 + 1;
            }*/
        }
};

