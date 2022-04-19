get_sample.marker <- function(assigned_seurat, sample.name, assay, out_lambda = NULL) { # nolint
    markers <-
        assigned_seurat@assays[[assay]]@misc$sample.marker[[sample.name]]
    if (is.null(out_lambda)) {
        return(markers)
    } else {
        out_lambda(markers)
        return(assigned_seurat)
    }
}