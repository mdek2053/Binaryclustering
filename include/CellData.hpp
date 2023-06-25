#pragma once

#include <cstdint>
#include <vector>

class CellData {
public:
    int size;
    std::vector<uint64_t> chunks;
    // I have to rename this because it conflicts with the annoylib definitions.
    int count;
    CellData() {}
    CellData(int size, const std::vector<uint64_t> &chunks, int count) : size(size), chunks(chunks), count(count) {}
};