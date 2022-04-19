read.seurat <- function(filepath) {
    library(SeuratDisk)
    library(dplyr)
    message("Start reading seurat.")
    seurat <-
        SeuratDisk::LoadH5Seurat(
            filepath
        )
    message("Finished reading seurat.")
    return(seurat)
}