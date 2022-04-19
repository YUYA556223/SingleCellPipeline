run_single_seupipe <- function(seurat, resolution = 0.8) {
    seurat %>%
        Seurat::RunPCA(
            npcs = 30, verbose = TRUE
        ) %>%
        Seurat::RunUMAP(
            reduction = "pca", dims = 1:30,
            verbose = TRUE
        ) %>%
        Seurat::FindNeighbors(
            reduction = "umap",
            dims = 1:2, verbose = TRUE
        ) %>%
        Seurat::FindClusters(
            reduction.type = "umap",
            algorithm = 1, resolution = resolution, verbose = TRUE
        ) %>%
        return()
}