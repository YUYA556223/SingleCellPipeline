where <- function(seurat, condition) {
    return(subset(
        seurat,
        condition
    ))
}