#pragma once

#include <cstdint>
#include <functional>
#include <thread>
#include <vector>

#include "CellData.hpp"
#include "Distance.hpp"
#include "Similarity.hpp"

// TODO: convert this to a ((n * (n - 1)) / 2) 1D array for optiomal spacing.
class Matrix {
public:
    int n;
    int m;
    float** data;

    Matrix(int n, int m) : n(n), m(m) {
        data = new float*[n];
        for (int i = 0; i < n; i++) {
            data[i] = new float[m];
        }
    }

    // destructor
    ~Matrix() {
        for (int i = 0; i < n; i++) {
            delete[] data[i];
        }
        delete[] data;
    }

    // copy constructor
    Matrix(const Matrix& other) : n(other.n), m(other.m) {
        data = new float*[n];
        for (int i = 0; i < n; i++) {
            data[i] = new float[m];
            std::copy(other.data[i], other.data[i] + m, data[i]);
        }
    }

    // copy assignment operator
    Matrix& operator=(const Matrix& other) {
        if (this != &other) {
            for (int i = 0; i < n; i++) {
                delete[] data[i];
            }
            delete[] data;

            n = other.n;
            m = other.m;
            data = new float*[n];
            for (int i = 0; i < n; i++) {
                data[i] = new float[m];
                std::copy(other.data[i], other.data[i] + m, data[i]);
            }
        }
        return *this;
    }
};

class DistanceMatrix : public Matrix {
public:
    DistanceMatrix(int n) : Matrix(n, n){};
};

class SimilarityMatrix : public Matrix {
public:
    SimilarityMatrix(int n) : Matrix(n, n){};
};

void fillMatrix(Matrix& matrix, const std::vector<CellData>& cellData, std::function<float(const CellData&, const CellData&)> comparisonFunction) {
    for (int i = 0; i < matrix.n; i++) {
        for (int j = i + 1; j < matrix.m; j++) {
            float result = comparisonFunction(cellData[i], cellData[j]);
            matrix.data[i][j] = result;
            matrix.data[j][i] = result;
        }
    }
}

SimilarityMatrix similarityMatrix(int nCells, const std::vector<CellData>& cellData, SimilarityMetric metric) {
    SimilarityFunction function = getSimilarityFunction(metric);
    SimilarityMatrix matrix(nCells);
    fillMatrix(matrix, cellData, function);
    return matrix;
}

DistanceMatrix distanceMatrix(int nCells, const std::vector<CellData>& cellData, DistanceMetric metric) {
    auto function = getDistanceFunction(metric);
    DistanceMatrix matrix(nCells);
    fillMatrix(matrix, cellData, function);
    return matrix;
}