#include <iostream>
#include <vector>
#include <chrono>
#include <memory>

class Enemy {
public:
    Enemy() {}
    ~Enemy() {}
};

int main() {
    const int amount = 1000000;

    // Measure time with raw pointers
    auto start = std::chrono::steady_clock::now();
    {
        std::vector<Enemy*> enemies;
        for (size_t i = 0; i < amount; ++i) {
            enemies.push_back(new Enemy());
        }

        // we clear the vecotr for the next use
        for (Enemy* itr : enemies) {
            delete itr;
        }
        enemies.clear();
    }
    auto end = std::chrono::steady_clock::now();
    std::cout << "Time with raw pointers: " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << " ms\n";
    
    // ========================================================================================================================================================

    // Measure time with smart pointers
    start = std::chrono::steady_clock::now();
    {
        std::vector<std::unique_ptr<Enemy>> enemies;
        for (int i = 0; i < amount; ++i) {
            enemies.push_back(std::make_unique<Enemy>());
        }
        // no need to delete since smart pointers destruct on leaving this scope
    }

    end = std::chrono::steady_clock::now();
    std::cout << "Time with smart pointers: " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << " ms\n";

    return 0;
}
