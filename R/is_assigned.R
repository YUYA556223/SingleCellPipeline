is_assigned <- function(seurat) {
    return(
        "celltype" %in% colnames(seurat@meta.data)
    )
}