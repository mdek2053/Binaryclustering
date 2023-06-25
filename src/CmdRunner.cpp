#include <sstream>

#include "include/BinaryMath.hpp"
#include "include/CellData.hpp"
#include "include/Clustering.hpp"
#include "include/Distance.hpp"
#include "include/DynamicBitset.hpp"
#include "include/Evaluation.hpp"
#include "include/Graph.hpp"
#include "include/Io.hpp"
#include "include/KNNGraph.hpp"
#include "include/Matrix.hpp"
#include "include/Preprocess.hpp"
#include "include/Similarity.hpp"
#include "include/StopWatch.hpp"

void runClustering(int threads, std::string dataFilepath, SimilarityMetric metric, CommunityDetection cd, int k, float resolution) {
    // Stopwatch to time the algorithm
    StopWatch sw = StopWatch();
    StopWatch fw = StopWatch();
    fw.start();

    std::cout << "Start reading data..." << std::endl;
    sw.start();
    auto binaryData = readCSV(dataFilepath);
    std::cout << "Finished preprocessing data..." << std::endl;
    sw.print();
    std::cout << std::endl;

    std::cout << "Starting clustering pipeline..." << std::endl;

    // Reading code moved to another method to avoid repeated reading of data.
    int nCells = binaryData.size();
    int nGenes = binaryData[0].size();
    std::cout << "Input data consists of " << nCells << " cells and " << nGenes << " genes." << std::endl;

    // Preprocess data
    std::cout << "Preprocessing data..." << std::endl;
    sw.start();
    std::vector<CellData> data = preprocess(nCells, nGenes, binaryData);
    std::cout << "Finished preprocessing data..." << std::endl;
    sw.print();
    std::cout << std::endl;

    // Create KNN graph
    std::cout << "Creating KNN graph..." << std::endl;
    sw.start();
    Graph g = knnGraph(nCells, data, k, metric, threads);
    std::cout << "Finished creating KNN graph." << std::endl;
    sw.print();
    std::cout << std::endl;

    // Create clusters
    std::cout << "Detecting communities..." << std::endl;
    sw.start();
    std::vector<int> clusters = clusterize(g, cd, resolution);
    std::cout << "Finished community detection." << std::endl;
    sw.print();
    std::cout << std::endl;

    // Output results
    std::cout << "Results:" << std::endl;
    outputClusterSizes(clusters);
    std::cout << std::endl;

    // Terminating algorithm.
    std::cout << "Terminating clustering pipeline..." << std::endl;
    fw.print();
}

int main(int argc, char** argv) {
    // Default values
    int threads = 1;
    std::string dataFilepath = "data/dataset_867/binaryExpression.csv";
    SimilarityMetric metric = Cosine_Similarity;
    CommunityDetection cd = Louvain;
    int k = 10;
    float resolution_parameter = 0.2;

    for (int i = 1; i < argc; i++) {
        std::string arg = argv[i];
        if (arg == "-DATA") {
            if (i + 1 == argc) {
                std::cerr << "-DATA is not specified, no more arguments found\n" << std::endl;
            }

            dataFilepath = argv[++i];
        } else if (arg == "-SIM") {
            if (i + 1 == argc) {
                std::cerr << "-SIM is not specified, no more arguments found\n" << std::endl;
            }

            metric = metricFromName(argv[++i]);
        } else if (arg == "-COM") {
            if (i + 1 == argc) {
                std::cerr << "-COM is not specified, no more arguments found\n" << std::endl;
            }

            cd = communityFromName(argv[++i]);
        } else if (arg == "-THREADS") {
            if (i + 1 == argc) {
                std::cerr << "-THREADS is not specified, no more arguments found\n" << std::endl;
            }

            std::stringstream ss(argv[++i]);
            if (!(ss >> threads)) {
                std::cerr << "Invalid number of threads.\n";
                return -1;
            }
        } else if (arg == "-K") {
            if (i + 1 == argc) {
                std::cerr << "-K is not specified, no more arguments found\n" << std::endl;
            }

            std::stringstream ss(argv[++i]);
            if (!(ss >> k)) {
                std::cerr << "Invalid number of K parameter.\n";
                return -1;
            }
        } else {
            std::cerr << "Unknown argument: " << arg << "\n";
            return -1;
        }
    }

    runClustering(threads, dataFilepath, metric, cd, k, resolution_parameter);

    return 0;
}