#pragma once

#include <cstdint>
#include <vector>

class DynamicBitset {
public:
    const int CHUNK_SIZE = 64;
    std::vector<uint64_t> data;
    int bit_size;

    DynamicBitset() {}

    DynamicBitset(int bit_size_) : bit_size(bit_size_) {
        int size = bit_size / CHUNK_SIZE;
        if (bit_size % CHUNK_SIZE != 0) {
            ++size;
        }
        data = std::vector<uint64_t>(size);
    }

    DynamicBitset& operator=(const DynamicBitset& other) {
        if (this != &other) {
            // Perform deep copy of the data vector
            data = other.data;
            bit_size = other.bit_size;
        }
        return *this;
    }

    int size() {
        return bit_size;
    }

    bool getBit(int i) {
        int chunk = i / CHUNK_SIZE;
        int offset = i % CHUNK_SIZE;

        return data[chunk] & (1 << offset);
    }

    void setBit(int i, bool v) {
        int chunk = i / CHUNK_SIZE;
        int offset = i % CHUNK_SIZE;

        if (v) {
            // 1
            data[chunk] |= (1 << offset);
        } else {
            // 0
            data[chunk] &= ~(1 << offset);
        }
    }
};