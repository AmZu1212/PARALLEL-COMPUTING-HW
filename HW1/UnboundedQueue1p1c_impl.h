#include "UnboundedQueue1p1c.h"
#include <atomic>
#include <thread>
#include <iostream>
#include <chrono>

namespace chrono = std::chrono;
// need to speed this up
struct Node {
    int value;
    Node* next;
    Node(int val = 0) : value(val), next(nullptr) {}
};


class UnboundedQueue1p1c : public UnboundedQueue1p1cAbstract {
private:
    std::atomic<Node*> head;
    std::atomic<Node*> tail;
    std::atomic<int> count;

public:
    UnboundedQueue1p1c() {// make the dummy
        Node* dummy = new Node();
        head.store(dummy);
        tail.store(dummy);
        count.store(0);
    }

    ~UnboundedQueue1p1c() {// kill every node
        while (head.load() != nullptr) {
            Node* tmp = head.load();
            head.store(head.load()->next);
            delete tmp;
        }
    }

    int size() const override {
        return count.load();
    }

    bool pop(int &val) override {// remove head, and update pointer
        Node* head_node = head.load();
        Node* next_node = head_node->next;

        if (next_node == nullptr) {
            return false; // Queue is empty
        }

        val = next_node->value;
        head.store(next_node);
        delete head_node;
        count.fetch_sub(1); // this is -- in atomics
        return true;
    }

    void push(int value) override {
        Node* new_node = new Node(value);
        Node* old_tail = tail.exchange(new_node);
        old_tail->next = new_node;
        count.fetch_add(1); // this is ++ in atomics
    }  
};
