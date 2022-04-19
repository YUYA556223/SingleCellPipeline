subscribe <- function(seurat, out_lambda) {
    if (!is.null(out_lambda)) {
        out_lambda(seurat)
    }
    return(seurat)
}