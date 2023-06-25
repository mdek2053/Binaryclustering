#pragma once

#include <chrono>
#include <iostream>

class StopWatch {
    std::chrono::steady_clock::time_point startTime;

public:
    void start() { startTime = std::chrono::steady_clock::now(); }

    void print() {
        auto endTime = std::chrono::steady_clock::now();
        auto elapsedTime = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime);
        std::cout << "Time elapsed: " << elapsedTime.count() / 1000.0 << " seconds." << std::endl;
    }
};