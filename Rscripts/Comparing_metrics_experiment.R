library(Rcpp)

Sys.setenv(CXX = "/usr/bin/g++")
Sys.setenv(PKG_CXXFLAGS = "-fdiagnostics-color=always -Ofast -march=native",
           PKG_CPPFLAGS = "-I/home/m/Downloads/vcpkg/installed/x64-linux/include -I/usr/local/include -I/home/m/rp/workspace/binary_clustering_final/binary_clustering -I/home/m/R/x86_64-pc-linux-gnu-library/4.3/Rcpp/include -I/home/m/R/x86_64-pc-linux-gnu-library/4.3/Rcpp/include -I/home/m/R/x86_64-pc-linux-gnu-library/4.3/RcppAnnoy/include",
           PKG_LIBS = "-L/home/m/Downloads/vcpkg/installed/x64-linux/lib -L/usr/local/lib /home/m/R/x86_64-pc-linux-gnu-library/4.3/Rcpp/libs/Rcpp.so -ligraph -lboost_system -lboost_thread -lR -pthread")

currentDir <- "/home/m/rp/workspace/binary_clustering_final/binary_clustering/"
dataDir <- "/home/m/rp/workspace/binary_clustering_v2/binaryclustering/"

sourceCpp(file.path(currentDir, "src/Rinterface.cpp"))

dataset <- "1000x2000"

dataPath <- file.path(dataDir, "Binary clustering/data/HLCA_samples/", paste0(dataset, ".csv"))
outputPath <- file.path(currentDir, "Binary clustering/output/comparing_algorithms/HLCA_samples/", dataset, "binary_exact")

metrics <- c("Bray-Curtis", "Ochiai", "Dice", "Euclidean", "Manhattan", "Jaccard", "Kulsinski", "Pearson", "YuleQ", "Zero")

sampleSize <- 10

run <- function(dataPath, sampleSize, outputPath, metrics, write = TRUE) {

    randomMetrics <- sample(metrics)

    for (i in 1:length(metrics)) {
        read(dataPath)

        preprocess()

        startTime <- Sys.time()

        findNeighboursExact(randomMetrics[i], 30)

        endTime <- Sys.time()

        executionTime <- difftime(endTime, startTime, units = "secs")
        if(write) {
            write(executionTime, file = file.path(outputPath, paste0(randomMetrics[i], ".txt")), append = TRUE)
        }
    }
    
}

if (!dir.exists(outputPath)) {
    dir.create(outputPath, recursive = TRUE)
}

for (i in 1:sampleSize) {

    run(dataPath, sampleSize, outputPath, metrics)
}