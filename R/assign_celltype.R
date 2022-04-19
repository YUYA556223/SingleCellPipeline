assign_celltype <- function(seurat, annotation, cell.meta, sample.meta) { # nolint
    library(ggplot2, quietly = TRUE)
    library(Seurat, quietly = TRUE)
    library(SeuratDisk, quietly = TRUE)
    library(dplyr, quietly = TRUE)
    if (SingleCellPipeline::is_assigned(seurat)) {
        warning("Given seurat object is already seurat")
        return(seurat)
    }
    type_assigned <- annotation$celltype
    names(type_assigned) <- levels(seurat)
    seurat <-
        Seurat::RenameIdents(
            seurat, type_assigned
        )
    seurat@meta.data <-
        tibble::rownames_to_column(
            seurat@meta.data
        ) %>%
        dplyr::rename(
            cell = rowname
        ) %>%
        dplyr::mutate(
            seurat_clusters = as.integer(seurat_clusters) - 1
        ) %>%
        dplyr::left_join(
            annotation,
            by = c("seurat_clusters" = "cluster")
        ) %>%
        dplyr::left_join(
            cell.meta %>%
                dplyr::rename(
                    color.celltype = color
                ),
            by = "celltype"
        ) %>%
        dplyr::left_join(
            sample.meta %>%
                dplyr::rename(
                    color.sample = color
                ),
            by = "sample"
        ) %>%
        dplyr::mutate(
            celltype = factor(
                celltype,
                levels = cell.meta$celltype
            ),
            sample = factor(
                sample,
                levels = sample.meta$sample
            )
        )
    rownames(seurat@meta.data) <-
        seurat@meta.data$cell
    seurat@meta.data <-
        seurat@meta.data %>%
        dplyr::select(-cell)
    order_celltype <-
        cell.meta %>%
        dplyr::distinct(celltype)
    diff <-
        setdiff(levels(seurat), order_celltype$celltype)
    cnt <- 0
    diff_str <- ""
    for (diff_el in diff) {
        if (cnt == 0) {
            diff_str <- diff_el
        } else {
            diff_str <- paste(
                diff_str,
                diff_el,
                sep = ", "
            )
        }
        cnt <- cnt + 1
    }
    if (length(diff) != 0) {
        stop(
            paste(
                "Cluster type [",
                diff_str,
                " ] is unknown."
            )
        )
    }
    levels(seurat) <- order_celltype$celltype
    return(seurat)
}