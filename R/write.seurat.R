write.seurat <- function(seurat, filepath, overwrite = TRUE) {
    library(SeuratDisk)
    library(dplyr)
    message("Start writing seurat.")
    seurat %>%
        SeuratDisk::SaveH5Seurat(
            filename = filepath,
            overwrite = overwrite
        )
    message("Finished writing seurat.")
    return(seurat)
}