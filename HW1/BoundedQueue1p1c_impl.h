#include "BoundedQueue1p1c.h"
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


class BoundedQueue1p1c : BoundedQueueAbstract_1p1c {
private:
    std::atomic<Node*> head;
    std::atomic<Node*> tail;
    std::atomic<int> count;
    const int capacity; // no need for atomic this is constant

public:
    BoundedQueue1p1c(int capacity) {
        this->capacity = capacity;
        Node* dummy = new Node();
        head.store(dummy);
        tail.store(dummy);
        count.store(0);
    }

    ~BoundedQueue1p1c() {
        while (head.load() != nullptr) {
            Node* tmp = head.load();
            head.store(head.load()->next);
            delete tmp;
        }
    }

    // pop is the same as unbounded
    bool pop(int &val) override {
        Node* head_node = head.load();
        Node* next_node = head_node->next;

        if (next_node == nullptr) {
            // Queue is empty
            return false; 
        }

        val = next_node->value;
        head.store(next_node);
        delete head_node;
        count.fetch_sub(1);
        return true;
    }

    bool push(int value) override {
        if (count.load() >= capacity) {
            // Queue is full
            return false;
        }

        Node* new_node = new Node(value);
        Node* old_tail = tail.exchange(new_node);
        old_tail->next = new_node;
        count.fetch_add(1);

        return true;
    }

    int size() const override {
        return count.load();
    }
};
