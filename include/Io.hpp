#pragma once

#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "DynamicBitset.hpp"
#include "Evaluation.hpp"
#include "Preprocess.hpp"

std::vector<std::string> readCells(const std::string& filepath) {
    std::ifstream infile(filepath);
    std::vector<std::string> cells;
    for (std::string line; getline(infile, line);) {
        cells.push_back(line);
    }
    infile.close();
    return cells;
}

std::vector<std::string> readGenes(const std::string& filepath) {
    std::ifstream infile(filepath);
    std::vector<std::string> genes;
    for (std::string line; std::getline(infile, line);) {
        genes.push_back(line);
    }
    infile.close();
    return genes;
}

std::vector<DynamicBitset> readBits(int nCells, int nGenes, const std::string& filepath) {
    std::ifstream infile(filepath);
    std::vector<DynamicBitset> bits(nCells);
    for (int i = 0; i < nCells; i++) {
        bits[i] = DynamicBitset(nGenes);
    }

    int gene = 0;
    for (std::string line; std::getline(infile, line);) {
        int cell = 0;
        // it is faster just read char not reference (&) to it
        for (char c : line) {
            if (c == '0') {
                bits[cell].setBit(gene, 0);
            } else if (c == '1') {
                bits[cell].setBit(gene, 1);
            } else if (c == ',') {
                cell++;
            }
        }

        if (cell != nCells - 1) {
            throw std::runtime_error("The amount of cell data rows doesn't coincide with the required cell count for gene row: " + gene);
        }

        gene++;
    }

    if (gene != nGenes) {
        throw std::runtime_error("The amount of gene data rows doesn't coincide with the required gene count.");
    }

    return bits;
}

int countGenes(std::ifstream& file) {
    int genes = 0;
    std::string line;
    while (std::getline(file, line)) {
        genes++;
    }

    file.clear();
    file.seekg(0, std::ios::beg);
    return genes;
};

int countCells(std::ifstream& file) {
    int cells = 1;
    std::string line;
    std::getline(file, line);

    for (char& c : line) {
        if (c == ',') cells++;
    }
    file.clear();
    file.seekg(0, std::ios::beg);
    return cells;
};

std::vector<DynamicBitset> readCSV(const std::string& filepath) {
    std::ifstream file(filepath);

    if (!file.is_open()) {
        throw std::runtime_error("Cannot open file!");
    }

    int nGenes = countGenes(file);
    int nCells = countCells(file);

    std::vector<DynamicBitset> bits(nCells);

    for (int i = 0; i < nCells; i++) {
        bits[i] = DynamicBitset(nGenes);
    }

    int gene = 0;
    for (std::string line; getline(file, line);) {
        int cell = 0;
        std::istringstream iss(line);
        std::string value;
        while (getline(iss, value, ',')) {
            float floatValue = stof(value);
            bits[cell].setBit(gene, (floatValue > 0.0f) ? 1 : 0);
            cell++;
        }

        gene++;
    }

    return bits;
}

std::vector<std::vector<double>> readNonbinary(int nCells, int nGenes, const std::string& file) {
    std::ifstream infile(file);
    std::vector<std::vector<double>> expressions(nCells);
    for (int i = 0; i < nCells; i++) {
        expressions[i] = std::vector<double>(nGenes);
    }

    int gene = 0;
    for (std::string line; std::getline(infile, line);) {
        int cell = 0;
        std::string expr = "";
        for (char& c : line) {
            if (c == ';') {
                expressions[cell++][gene] = std::stod(expr);
                expr = "";
            } else {
                expr.push_back(c);
            }
        }

        if (cell != nCells - 1) {
            throw std::runtime_error("The amount of cell data rows doesn't coincide with the required cell count for gene row: " + gene);
        }

        expressions[cell][gene++] = std::stod(expr);
    }

    if (gene != nGenes) {
        throw std::runtime_error("The amount of gene data rows doesn't coincide with the required gene count.");
    }

    return expressions;
}

void outputClusterSizes(const std::vector<int>& clusters) {
    std::map<int, int> counts;
    for (int i = 0; i < (int)clusters.size(); i++) {
        if (counts.find(clusters[i]) == counts.end()) {
            counts[clusters[i]] = 1;
        } else {
            counts[clusters[i]]++;
        }
    }

    std::cout << "Cluster sizes:" << std::endl;
    for (std::pair<const int, int> cnt : counts) {
        std::cout << "Cluster " << cnt.first << " contains " << cnt.second << " cells." << std::endl;
    }
}

void outputARI(int nCells, const std::vector<std::string>& cells, const std::vector<int>& clustersPred) {
    std::vector<int> clustersTrue = groundTruth(cells);
    double ari = calculateARI(nCells, clustersPred, clustersTrue);
    std::cout << "ARI: " << ari << std::endl;
}

void outputARI(int nCells, const std::vector<int>& clustersTrue, const std::vector<int>& clustersPred) {
    double ari = calculateARI(nCells, clustersPred, clustersTrue);
    std::cout << "ARI: " << ari << std::endl;
}