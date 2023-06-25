library(Seurat)
currentDir <- "/home/m/rp/workspace/binary_clustering_final/binary_clustering/Binary clustering/"
dataDir <- "/home/m/rp/workspace/binary_clustering_v2/binaryclustering/Binary clustering/"

dataset <- "1000x500"

dataPath <- file.path(dataDir, "data/HLCA_samples/", paste0(dataset, ".csv"))
outputPath <- file.path(currentDir, "output/comparing_algorithms/HLCA_samples/", dataset, "seurat")

sampleSize <- 1

run <- function(dataPath, sampleSize, outputPath) {

    data = read.csv(dataPath)

    for (i in 1:sampleSize) {

        seurat_object <- CreateSeuratObject(data)
        
        startTime <- Sys.time()

        all.genes <- 1:nrow(data)
        seurat_object <- NormalizeData(seurat_object)
        seurat_object <- ScaleData(seurat_object, features = all.genes)

        intermediateTime1 <- Sys.time()

        seurat_object <- RunPCA(seurat_object, features = all.genes)

        intermediateTime2 <- Sys.time()

        seurat_object <- FindNeighbors(seurat_object, dims = 1:10, k.param = 30)

        endTime <- Sys.time()

        totalTime <- difftime(endTime, startTime, units = "secs")
        prepTime <- difftime(intermediateTime1, startTime, units = "secs")
        PCATime <- difftime(intermediateTime2, intermediateTime1, units = "secs")
        kNNTime <- difftime(endTime, intermediateTime2, units = "secs")

        write(totalTime, file = file.path(outputPath, "total.txt"), append = TRUE)
        write(prepTime, file = file.path(outputPath, "prep.txt"), append = TRUE)
        write(PCATime, file = file.path(outputPath, "pca.txt"), append = TRUE)
        write(kNNTime, file = file.path(outputPath, "knn.txt"), append = TRUE)
    }
}

if (!dir.exists(outputPath)) {
    dir.create(outputPath, recursive = TRUE)
}
run(dataPath, sampleSize, outputPath)

