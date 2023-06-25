#pragma once

#include <cstdint>
#include <vector>

// const uint64_t m1 = 0x5555555555555555;
// const uint64_t m2 = 0x3333333333333333;
// const uint64_t m4 = 0x0f0f0f0f0f0f0f0f;
// const uint64_t m8 = 0x00ff00ff00ff00ff;
// const uint64_t m16 = 0x0000ffff0000ffff;
// const uint64_t m32 = 0x00000000ffffffff;
// const uint64_t h01 = 0x0101010101010101;

// int popcount(uint64_t x) {
// 	x -= (x >> 1) & m1;
// 	x = (x & m2) + ((x >> 2) & m2);
// 	x = (x + (x >> 4)) & m4;
// 	return (x * h01) >> 56;
// }

// Fast
int popcount(uint64_t x) {
    return __builtin_popcount(x);
}

int popcount(const std::vector<uint64_t> &bits) {
    int count = 0;
    for (int k = 0; k < (int)bits.size(); k++) {
        count += popcount(bits[k]);
    }
    return count;
}

int ANDCount(const std::vector<uint64_t> &a, const std::vector<uint64_t> &b) {
    int count = 0;
    for (int k = 0; k < (int)a.size(); k++) {
        count += popcount(a[k] & b[k]);
    }
    return count;
}

int XORCount(const std::vector<uint64_t> &a, const std::vector<uint64_t> &b) {
    int count = 0;
    for (int k = 0; k < (int)a.size(); k++) {
        count += popcount(a[k] ^ b[k]);
    }
    return count;
}