set_sample.marker <- function(assigned_seurat, assay, sample.name, markers) { # nolint
    assigned_seurat@assays[[assay]]@misc$sample.marker[[sample.name]] <- markers
    return(assigned_seurat)
}