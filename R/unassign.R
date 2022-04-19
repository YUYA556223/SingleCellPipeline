unassign <- function(assigned_seurat) {
    if (SingleCellPipeline::is_assigned(assigned_seurat)) {
        SeuratObject::Idents(assigned_seurat) <-
            assigned_seurat@meta.data$seurat_clusters
        assigned_seurat@meta.data <-
            assigned_seurat@meta.data %>%
            dplyr::select(-celltype, -evidence, color.sample)
        message("Successfly finalized unassign operation")
        return(assigned_seurat)
    }
}