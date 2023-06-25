#include <Rcpp.h>
#include "experiments/Common.hpp"

int nCells;
int nGenes;
vector<DynamicBitset> bits;
vector<CellData> cells;
Graph knn;

// [[Rcpp::export]]
void read(string dataPath) {
    bits = readCSV(dataPath);
    nCells = bits.size();
    nGenes = bits[0].size();
    cout << "Read csv with " << nCells << " Cells and " << nGenes << " Genes." << endl;
}

// [[Rcpp::export]]
void preprocess() {
    cells = preprocess(nCells, nGenes, bits);
}

// [[Rcpp::export]]
void findNeighboursExact(string metric, int k) {
    // FIX
	if (metric == "Bray-Curtis") {
		knn = knnGraph(nCells, cells, k, Sorensen_Similarity, thread::hardware_concurrency());
		return;
	}
	if (metric == "Ochiai") {
        knn = knnGraph(nCells, cells, k, Ochiai_Distance, thread::hardware_concurrency());
		return;
	}
	if (metric == "Dice") {
        knn = knnGraph(nCells, cells, k, Dice_Similarity, thread::hardware_concurrency());
		return;
	}
	if (metric == "Euclidean") {
        knn = knnGraph(nCells, cells, k, Euclidean_Distance, thread::hardware_concurrency());
		return;
	}
	if (metric == "Manhattan") {
        knn = knnGraph(nCells, cells, k, Manhattan_Distance, thread::hardware_concurrency());
		return;
	}
	if (metric == "Jaccard") {
        knn = knnGraph(nCells, cells, k, Jaccard_Distance, thread::hardware_concurrency());
		return;
	}
	if (metric == "Kulsinski") {
        knn = knnGraph(nCells, cells, k, Kulsinski_Distance, thread::hardware_concurrency());
		return;
	}
	if (metric == "Pearson") {
        knn = knnGraph(nCells, cells, k, Pearson_Similarity, thread::hardware_concurrency());
		return;
	}
	if (metric == "YuleQ") {
        knn = knnGraph(nCells, cells, k, YuleQ_Similarity, thread::hardware_concurrency());
		return;
	}
	if (metric == "Zero") {
        knn = knnGraph(nCells, cells, k, Zero_Distance, thread::hardware_concurrency());
		return;
	}

}

// [[Rcpp::export]]
void findNeighboursApprox(int k) {
    cout << "Create knn graph."  << endl;
    knn = approxKnnGraph(nCells, nGenes, cells, k, thread::hardware_concurrency());
    cout << "Created knn graph."  << endl;
}

// [[Rcpp::export]]
vector<int> cluster() {
	return clusterize(knn, Louvain, 0.2);
}