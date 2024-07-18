#include "tree_barrier.h"
#include <mutex>
#include <condition_variable>
#include <iostream>

class BinaryTreeBarrier : BinaryTreeBarrierAbstract
{
    public:
         void barrier(int numberOfThreads)
         {
            // int arrived[]; <- the pings.
            // a node pings itself when their sons ping it

            // how do we put a binary tree into an array
         }


         void ResetTree(){}
};




class Node
{
    public:
    bool ready;

    // initializes a node with the "initial" value (true/false)
    Node(bool initial = false) : ready(initial){}
    
    void SetReady(bool value = true)
    {
        ready = value;
    }
};


class BarrierNode
{
    public:
    Node right;
    Node left;
    bool barrierReady = false;
    
    BarrierNode(bool rightValue = false, bool leftValue = false)
    {
        right = new Node(rightValue);
        left = new Node(leftValue);
    }

    bool UpdateReadiness()
    {
        if(right.ready && left.ready)
        {
            barrierReady = true;
        } 
        return barrierReady;
    }
};





class Barrier
{
    private:
        int numberOfThreads;
        int threadsWaiting;
        std::mutex mutex;
        std::condition_variable CV;
    public:
        Barrier(int numberOfThreads) : numberOfThreads(numberOfThreads), threadsWaiting(0){}
        void Wait(){
            std::unique_lock<std::mutex> lock(mutex);
            ++threadsWaiting;
            if(threadsWaiting < numberOfThreads)
            {
                CV.wait(lock, [this] {return threadsWaiting >= numberOfThreads;}); // return when everyone arrives
            }
            else
            {
                threadsWaiting = 0;
                CV.notify_all();// release the waiting threads
            }
        }
};
