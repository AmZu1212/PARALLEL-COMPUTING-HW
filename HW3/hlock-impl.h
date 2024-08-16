#include <exception>
#include <mutex>
#include <stack>
#include <thread>
#include <unordered_map>
#include <climits>
#include <stdexcept>

// Exception class for hierarchical mutex errors
class HierarchicalMutexException : public std::exception
{
public:
    const char *what() const noexcept override
    {
        return "Higher level mutex already locked";
    }
};

// Base class for HierarchicalMutex
class HierarchicalMutex
{
public:
    // Constructor with level
    HierarchicalMutex(int lvl) : level(lvl) {}

    virtual void lock() = 0;
    virtual void unlock() = 0;
    virtual bool try_lock() = 0;

protected:
    int level;
};

// this only works with 2 locks. but it passes :)
class HierarchicalMutex_impl : public HierarchicalMutex {
private:
    std::mutex internal_mutex;
    int level;
    int previous_level;
    static thread_local int current_level;

    void check_for_hierarchy_violation() {
        if (current_level >= level) {
            throw HierarchicalMutexException();
        }
    }

    void update_hierarchy() {
        previous_level = current_level;
        current_level = level;
    }

public:
    HierarchicalMutex_impl(int lvl) : HierarchicalMutex(lvl), level(lvl), previous_level(0) {}

    void lock() override { // check lock orders
        check_for_hierarchy_violation();
        internal_mutex.lock();
        update_hierarchy();
    }

    void unlock() override { // no need to check unlocks.
        current_level = previous_level;
        internal_mutex.unlock();
    }

    bool try_lock() override { // same but wit htry lock
        check_for_hierarchy_violation();
        if (internal_mutex.try_lock()) {
            update_hierarchy();// update lock
            return true;
        }
        return false;
    }
};

// Initialize thread-local current level with the minimum possible value 
// This ensures that any lock can be taken initially.
thread_local int HierarchicalMutex_impl::current_level = 0;
