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
//int thread_ids[MAX_NODES + 1];
std::mutex PRINTING_MUTEX;
std::mutex id_mutex;
thread_local int currentHeight = 1;
thread_local int currentPhase = 0;

class BarrierNode
{
    public:
        std::atomic<int> count[2]; // if this works, we can implement a phase lock
        std::condition_variable cv[2];
        std::mutex mtx[2];

        BarrierNode(){ count[0] = 0; count[1] = 0;}
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
            cout << "==============================================================" << endl;
            treeHeight = (int)log2(numOfThreads - 1);
            cout << "BinaryTreeBarrier(): tree height is: "<< treeHeight << "."<< endl;
            cout << "BinaryTreeBarrier(): number of threads is: "<< numOfThreads << "."<< endl;
            // The number of nodes in a binary tree is N-1, where N = number of threads.
            numOfNodes = numOfThreads - 1;
            cout << "BinaryTreeBarrier(): number of nodes is: "<< numOfNodes << "."<< endl;
            // Ensure num_nodes does not exceed MAX_NODES. REVISE THIS LATER
            if (numOfNodes > MAX_NODES)
                throw std::runtime_error("Number of nodes exceeds maximum limit.");
            cout << "=============================================================="<< endl;
        }

        void barrier() override
        {
            PRINTING_MUTEX.lock();
                //get the thread's id so we can map it to the correct barrier node.
                int threadID = thread_id;
                if(DEBUG) cout << "barrier(): entering threads ID is " << threadID << "."<< endl;
                int nodeID = ((numOfNodes + threadID - 1) / 2);//(threadID / 2) + treeHeight;
                if(DEBUG) cout << "barrier(): node ID assigned is " << nodeID << "." << endl;
            PRINTING_MUTEX.unlock();


            while (nodeID > 0)
            {
                auto &node = nodes[nodeID];
                std::unique_lock<std::mutex> lock(node.mtx[currentPhase]);
                node.count[currentPhase]++;

                int threads_at_this_level = (int)pow(2, currentHeight);
                PRINTING_MUTEX.lock();
                    cout << "barrier(): Thread " << threadID << " is checking " << nodeID << " (Phase = "<< currentPhase <<")." << endl;
                    cout << "barrier(): The Number of threads on this level is: " << threads_at_this_level << endl;
                PRINTING_MUTEX.unlock();

                if (node.count[currentPhase] < threads_at_this_level) {
                    PRINTING_MUTEX.lock();
                       cout << "barrier(): Thread " << threadID << " decided to wait at node " << nodeID << "."
                             <<" [Currently Waiting ("<< node.count[currentPhase] <<"/"<< threads_at_this_level << ")]"<< endl;
                    PRINTING_MUTEX.unlock();
                    node.cv[currentPhase].wait(lock/*, [&node, threads_at_this_level] { return node.count == threads_at_this_level; }*/);
                    PRINTING_MUTEX.lock();
                        cout << "barrier(): Thread " << threadID << " has stopped waiting at node " << nodeID << " (Phase = "<< currentPhase <<")." << endl;
                    PRINTING_MUTEX.unlock();
                } else {
                    PRINTING_MUTEX.lock();
                    cout << "barrier(): Thread " << threadID << " is notifying at node " << nodeID << " (Phase = "<< currentPhase <<").";
                    cout <<"..... [  Capacity was ("<< (node.count[currentPhase]) <<"/"<< threads_at_this_level << ")  ]"<< endl;
                    PRINTING_MUTEX.unlock();
                    node.cv[currentPhase].notify_all();
                }

                lock.unlock();
                PRINTING_MUTEX.lock();
                cout << "barrier(): Thread "<< threadID <<" is now moving to the next node number " <<  ((nodeID - 1) / 2) << "."<< endl;
                currentHeight++;
                nodeID = (nodeID - 1) / 2;
                PRINTING_MUTEX.unlock();
            }


            PRINTING_MUTEX.lock();
               cout << "barrier(): Thread "<< threadID <<" is passing the generic loop."<< endl;
            PRINTING_MUTEX.unlock();
            
            auto &root = nodes[0];
            std::unique_lock<std::mutex> lock(root.mtx[currentPhase]);
            root.count[currentPhase]++;
            
            if (root.count[currentPhase] < numOfThreads) {
                PRINTING_MUTEX.lock();
                    cout << "barrier(): Thread " << threadID << " decided to wait at the root (Phase = "<< currentPhase <<")."
                            <<" [Currently Waiting ("<< root.count[currentPhase] <<"/"<< numOfThreads << ")]"<< endl;
                PRINTING_MUTEX.unlock();
                root.cv[currentPhase].wait(lock/*, [this, &root] { return root.count == numOfThreads; }*/);
                PRINTING_MUTEX.lock();
                    cout << "barrier(): Thread " << threadID << " done waiting at root (Phase = "<< currentPhase <<")."<< endl;
                PRINTING_MUTEX.unlock();
            } else {
                PRINTING_MUTEX.lock();
                    cout << "barrier(): Thread " << threadID << " notifying root (Phase = "<< currentPhase <<")...............";
                    cout <<" [Capacity was ("<< root.count[currentPhase] <<"/"<< numOfThreads << ")]"<< endl;
                PRINTING_MUTEX.unlock();
                

                for (auto currentNodeID = 0; currentNodeID < numOfNodes; ++currentNodeID)
                {
                    auto &currentNode = nodes[currentNodeID];
                    currentNode.count[currentPhase] = 0;      // reset all counts,
                    //currentNode.mtx.unlock();   // reset all mutexes 
                }

                root.cv[currentPhase].notify_all(); // wake everyone up and move on.
            }
            currentPhase = (currentPhase == 0) ? 1 : 0; // Toggle phase
            currentHeight = 1;
        }
};

