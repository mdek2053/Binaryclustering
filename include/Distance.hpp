#pragma once

#include <chrono>
#include <cmath>
#include <functional>
#include <iostream>
#include <stdexcept>

#include "BinaryMath.hpp"
#include "CellData.hpp"

using DistanceFunction = std::function<float(const CellData&, const CellData&)>;

enum DistanceMetric {
    Jaccard_Distance,
    Euclidean_Distance,
    Cosine_Distance,
    Ochiai_Distance,
    Manhattan_Distance,
    Kulsinski_Distance,
    Zero_Distance
};

std::string getName(DistanceMetric metric) {
    switch (metric) {
        case Jaccard_Distance:
            return "Jaccard distance";
        case Cosine_Distance:
            return "Cosine distance";
        case Ochiai_Distance:
            return "Ochiai distance";
        case Kulsinski_Distance:
            return "Kulsinski distance";
        case Manhattan_Distance:
            return "Manhattan distance";
        case Euclidean_Distance:
            return "Euclidean distance";
        case Zero_Distance:
            return "Zero distance";
        default:
            return "unknown";
    }
}

// ABOUT int64_t and integer bounds
// popcount and others limited by number of genes
// this will always smaller than 100mil we should not use int64_t etc

auto jaccardDistance = [](const CellData& left, const CellData& right) {
    int diff = ANDCount(left.chunks, right.chunks);

    return 1.f - (diff / (float)(left.count + right.count - diff));
};

auto euclideanDistance = [](const CellData& left, const CellData& right) {
    int diff = XORCount(left.chunks, right.chunks);
    return std::sqrt((float)diff);
};

auto ochiaiDistance = [](const CellData& left, const CellData& right) {
    int diff = ANDCount(left.chunks, right.chunks);
    return 1.f - (diff / std::sqrt((float)left.count * (float)right.count));
};

auto manhattanDistance = [](const CellData& a, const CellData& b) {
    return (float)XORCount(a.chunks, b.chunks);
};

auto kulsinskiDistance = [](const CellData& left, const CellData& right) {
    int a = ANDCount(left.chunks, right.chunks);

    int lcnt = left.count;
    int rcnt = right.count;
    int N = left.size;

    int numerator = lcnt + rcnt - 3 * a - N;
    int denominator = lcnt + rcnt - 2 * a - N;

    return numerator / (float)denominator;
};

auto zeroDistance = [](const CellData& left, const CellData& right) {
    return 0;
};

DistanceFunction getDistanceFunction(DistanceMetric metric) {
    switch (metric) {
        case Jaccard_Distance:
            return jaccardDistance;
        case Cosine_Distance:
        case Ochiai_Distance:
            return ochiaiDistance;
        case Kulsinski_Distance:
            return kulsinskiDistance;
        case Manhattan_Distance:
            return manhattanDistance;
        case Euclidean_Distance:
            return euclideanDistance;
        case Zero_Distance:
            return zeroDistance;
        default:
            throw std::runtime_error("Unknown distance metric");
    }
}