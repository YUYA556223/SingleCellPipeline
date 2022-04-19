get_cell.meta <- function(seurat, out_lambda = NULL) { # nolint
    if (is.null(out_lambda)) {
        return(seurat@meta.data)
    } else {
        out_lambda(seurat@meta.data)
        return(seurat)
    }
}