set_default_assay <- function(seurat, assay) {
    SeuratObject::DefaultAssay(seurat) <- assay
    return(seurat)
}