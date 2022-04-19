get_trajectory_pseudotime <- function(seurat, out_lambda = NULL) {
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
        return(pseudotime)
    } else {
        out_lambda(pseudotime)
        return(seurat)
    }
}