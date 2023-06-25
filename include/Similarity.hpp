#pragma once

#include <chrono>
#include <cmath>
#include <functional>
#include <iostream>
#include <stdexcept>

#include "BinaryMath.hpp"
#include "CellData.hpp"

using SimilarityFunction = std::function<float(const CellData&, const CellData&)>;

enum SimilarityMetric {
    Jaccard_Similarity,
    Pearson_Similarity,
    Matthews_Similarity,
    Phi_Similarity,
    Cosine_Similarity,
    Ochiai_Similarity,
    YuleQ_Similarity,
    Sorensen_Similarity,
    Dice_Similarity,
    Szymkiewicz_Simpson_Similarity,
    Sokal_Michener_Similarity,
    Rogers_Tanimoto_Similarity,
    Russel_Rao_Similarity,
    Sokal_Sneath_I_Similarity,
    Sokal_Sneath_II_Similarity
};

SimilarityMetric metricFromName(const std::string& name) {
    if (name == "Jaccard") return Jaccard_Similarity;
    if (name == "Pearson") return Pearson_Similarity;
    if (name == "Matthews") return Matthews_Similarity;
    if (name == "Phi") return Phi_Similarity;
    if (name == "Cosine") return Cosine_Similarity;
    if (name == "Ochiai") return Ochiai_Similarity;
    if (name == "YuleQ") return YuleQ_Similarity;
    if (name == "Sorensen") return Sorensen_Similarity;
    if (name == "Dice") return Dice_Similarity;
    if (name == "Szymkiewicz_Simpson") return Szymkiewicz_Simpson_Similarity;
    if (name == "Sokal_Michener") return Sokal_Michener_Similarity;
    if (name == "Rogers_Tanimoto") return Rogers_Tanimoto_Similarity;
    if (name == "Russel_Rao") return Russel_Rao_Similarity;
    if (name == "Sokal_Sneath_I") return Sokal_Sneath_I_Similarity;
    if (name == "Sokal_Sneath_II") return Sokal_Sneath_II_Similarity;

    throw std::runtime_error("Unknown name for similarity metric");
}

std::string getName(SimilarityMetric metric) {
    switch (metric) {
        case Jaccard_Similarity:
            return "Jaccard similarity";
        case Pearson_Similarity:
            return "Pearson similarity";
        case Matthews_Similarity:
            return "Matthews similarity";
        case Phi_Similarity:
            return "Phi similarity";
        case Cosine_Similarity:
            return "Cosine similarity";
        case Ochiai_Similarity:
            return "Ochiai similarity";
        case YuleQ_Similarity:
            return "Yule's Q similarity";
        case Sorensen_Similarity:
            return "Sorensen similarity";
        case Dice_Similarity:
            return "Dice similarity";
        case Szymkiewicz_Simpson_Similarity:
            return "Szymkiewicz Simpson similarity";
        case Sokal_Michener_Similarity:
            return "Sokal Michener similarity";
        case Rogers_Tanimoto_Similarity:
            return "Rogers Tanimoto similarity";
        case Russel_Rao_Similarity:
            return "Russel Rao similarity";
        case Sokal_Sneath_I_Similarity:
            return "Sokal Sneath I similarity";
        case Sokal_Sneath_II_Similarity:
            return "Sokal Sneath II similarity";
        default:
            throw std::runtime_error("Similarity metric unknown...");
    }
}

// ABOUT int64_t and integer bounds
// popcount and others limited by number of genes
// this will always smaller than 100mil we should not use int64_t etc

auto jaccardSimilarity = [](const CellData& left, const CellData& right) {
    int diff = ANDCount(left.chunks, right.chunks);
    return diff / (float)(left.count + right.count - diff);
};

auto ochiaiSimilarity = [](const CellData& left, const CellData& right) {
    int diff = ANDCount(left.chunks, right.chunks);
    return diff / std::sqrt((float)left.count * (float)right.count);
};

// using float will not result in big loss of information but is more robust and faster

auto phiSimilarity = [](const CellData& left, const CellData& right) {
    float a = ANDCount(left.chunks, right.chunks);
    float b = left.count - a;
    float c = right.count - a;
    float d = left.size - a - b - c;
    return (a * d - b * c) / std::sqrt((a + b) * (c + d) * (a + c) * (b + d));
};

auto pearsonSimilarity = [](const CellData& left, const CellData& right) {
    float a = ANDCount(left.chunks, right.chunks);
    float b = left.count - a;
    float c = right.count - a;
    float d = left.size - a - b - c;
    return (a * d - b * c) / std::sqrt(left.count * right.count * (b + d) * (c + d));
};

auto yuleQSimilarity = [](const CellData& left, const CellData& right) {
    float a = ANDCount(left.chunks, right.chunks);
    float b = left.count - a;
    float c = right.count - a;
    float d = left.size - a - b - c;
    return (a * d - b * c) / (a * d + b * c);
};

auto sorensenSimilarity = [](const CellData& a, const CellData& b) {
    int symdif = XORCount(a.chunks, b.chunks);
    return 1.0f - (symdif / (float)(a.count + b.count));
};

auto diceSimilarity = [](const CellData& left, const CellData& right) {
    int a = ANDCount(left.chunks, right.chunks);
    // int b = left.popcount - a;
    // int c = right.popcount - a;

    // return 2 * a / (float)(2 * a + b + c);
    // THIS two formulas are equal just expand b and c and 2*a is canceled?
    return 2 * a / (float)(left.count + right.count);
};

auto szymkiewiczSimpson = [](const CellData& left, const CellData& right) {
    int a = ANDCount(left.chunks, right.chunks);
    // int b = left.popcount - a;
    // int c = right.popcount - a;

    // return a / ((float)std::min(a + b, a + c));
    // AGAIN THIS two formulas are equal just expand b and c and aa is canceled?
    return a / (float)std::min(left.count, right.count);
};

auto sokalMichener = [](const CellData& left, const CellData& right) {
    int a = ANDCount(left.chunks, right.chunks);
    int b = left.count - a;
    int c = right.count - a;
    int d = left.size - a - b - c;
    return (float)(a + d) / (a + b + c + d);

    return (float)(a + d) / (left.size);
};

auto rogersTanimoto = [](const CellData& left, const CellData& right) {
    int a = ANDCount(left.chunks, right.chunks);
    int b = left.count - a;
    int c = right.count - a;
    int d = left.size - a - b - c;

    return (float)(a + d) / (a + 2 * b + 2 * c + d);
};

auto russellRao = [](const CellData& left, const CellData& right) {
    int a = ANDCount(left.chunks, right.chunks);
    int b = left.count - a;
    int c = right.count - a;
    int d = left.size - a - b - c;

    return (float)(a) / (a + b + c + d);
};

auto sokalSneathI = [](const CellData& left, const CellData& right) {
    int a = ANDCount(left.chunks, right.chunks);
    int b = left.count - a;
    int c = right.count - a;

    return (float)(a) / (a + 2 * b + 2 * c);
};

auto sokalSneathII = [](const CellData& left, const CellData& right) {
    int a = ANDCount(left.chunks, right.chunks);
    int b = left.count - a;
    int c = right.count - a;
    int d = left.size - a - b - c;

    return (float)(2 * a + 2 * d) / (2 * a + b + c + 2 * d);
};

SimilarityFunction getSimilarityFunction(SimilarityMetric metric) {
    switch (metric) {
        case Jaccard_Similarity:
            return jaccardSimilarity;
        case Cosine_Similarity:
        case Ochiai_Similarity:
            return ochiaiSimilarity;
        case Matthews_Similarity:
        case Phi_Similarity:
            return phiSimilarity;
        case Pearson_Similarity:
            return pearsonSimilarity;
        case YuleQ_Similarity:
            return yuleQSimilarity;
        case Sorensen_Similarity:
            return sorensenSimilarity;
        case Dice_Similarity:
            return diceSimilarity;
        case Szymkiewicz_Simpson_Similarity:
            return szymkiewiczSimpson;
        case Sokal_Michener_Similarity:
            return sokalMichener;
        case Rogers_Tanimoto_Similarity:
            return rogersTanimoto;
        case Russel_Rao_Similarity:
            return russellRao;
        case Sokal_Sneath_I_Similarity:
            return sokalSneathI;
        case Sokal_Sneath_II_Similarity:
            return sokalSneathII;
        default:
            throw std::runtime_error("Similarity metric unknown...");
    }
}