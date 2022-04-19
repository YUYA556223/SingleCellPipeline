get_trajectory_linages <- function(seurat, out_lambda = NULL) {
    path <-
        seurat@reductions$umap@misc$trajectory$linages
    pseudotime <-
        seurat@reductions$umap@misc$trajectory$pseudotime
    if (is.null(path)) {
        stop(
            "Could not find trajectory linage.",
            "Please run 'estimate_trajectory()' first."
        )
        return(seurat)
    }
    if (is.null(out_lambda)) {
        return(path)
    } else {
        out_lambda(path)
        return(seurat)
    }
}