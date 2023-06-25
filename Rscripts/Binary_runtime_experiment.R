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
outputPath <- file.path(currentDir, "Binary clustering/output/comparing_algorithms/HCLA_samples/", dataset, "binary_exact")

sampleSize <- 1

run_exact <- function(dataPath, sampleSize, outputPath) {


    for (i in 1:sampleSize) {

        read(dataPath)
        
        startTime <- Sys.time()

        preprocess()

        intermediateTime <- Sys.time()

        findNeighboursExact("Euclidean", 30)

        endTime <- Sys.time()

        totalTime <- difftime(endTime, startTime, units = "secs")
        prepTime <- difftime(intermediateTime, startTime, units = "secs")
        kNNTime <- difftime(endTime, intermediateTime, units = "secs")

        write(totalTime, file = file.path(outputPath, "total.txt"), append = TRUE)
        write(prepTime, file = file.path(outputPath, "prep.txt"), append = TRUE)
        write(kNNTime, file = file.path(outputPath, "knn.txt"), append = TRUE)
    }
}

run_approx <- function(dataPath, sampleSize, outputPath) {

    for (i in 1:sampleSize) {

        read(dataPath)
        
        startTime <- Sys.time()

        preprocess()

        intermediateTime <- Sys.time()

        findNeighboursApprox(30)

        endTime <- Sys.time()

        totalTime <- difftime(endTime, startTime, units = "secs")
        prepTime <- difftime(intermediateTime, startTime, units = "secs")
        kNNTime <- difftime(endTime, intermediateTime, units = "secs")
        

        write(totalTime, file = file.path(outputPath, "total.txt"), append = TRUE)
        write(prepTime, file = file.path(outputPath, "prep.txt"), append = TRUE)
        write(kNNTime, file = file.path(outputPath, "knn.txt"), append = TRUE)
    }
}



if (!dir.exists(outputPath)) {
    dir.create(outputPath, recursive = TRUE)
}

run_exact(dataPath, sampleSize, outputPath)
# run_approx(dataPath, sampleSize, outputPath)
