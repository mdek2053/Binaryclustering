#pragma once

#include <cstdint>
#include <vector>

#include "BinaryMath.hpp"
#include "CellData.hpp"
#include "DynamicBitset.hpp"

std::vector<uint64_t> chunkToLongs(int /* nGenes */, const DynamicBitset& bits) {
    // int chunkSize = 64;
    // int start = (nGenes % chunkSize) - chunkSize;
    // auto chunks = std::vector<uint64_t>((nGenes - start) / chunkSize);

    // int cur = 0;
    // for (int i = start; i < nGenes; i += chunkSize) {
    // 	uint64_t chunk = 0L;
    // 	for (int j = std::max(i, 0); j < i + chunkSize; j++) {
    // 		chunk = (chunk << 1) + bits[j];
    // 	}
    // 	chunks[cur++] = chunk;
    // }

    return bits.data;
}

std::vector<CellData> preprocess(int nCells, int nGenes, std::vector<DynamicBitset>& bits) {
    auto data = std::vector<CellData>(nCells);

    for (int i = 0; i < nCells; i++) {
        auto chunk = chunkToLongs(nGenes, bits[i]);
        int popcnt = popcount(chunk);
        data[i] = CellData(nGenes, chunk, popcnt);
    }

    bits.clear();
    return data;
}

std::vector<int> groundTruth(const std::vector<std::string>& cells) {
    std::vector<int> groundTruth(cells.size());
    for (int i = 0; i < (int)cells.size(); i++) {
        const auto& gt = cells[i];
        size_t first = gt.find(",");
        size_t second = gt.find(",", first + 1);
        auto id = gt.substr(second + 1);
        groundTruth[i] = atoi(id.c_str());
    }
    return groundTruth;
}