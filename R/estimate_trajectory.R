estimate_trajectory <- function(assigned_seurat, root_celltype, recalculate = FALSE) { # nolint
    message("@ Start executing estimate_trajectory >>>")
    library(ggplot2, quietly = TRUE)
    library(Seurat, quietly = TRUE)
    library(SeuratDisk, quietly = TRUE)
    library(dplyr, quietly = TRUE)
    library(slingshot, quietly = TRUE)
    library(SingleCellExperiment, quietly = TRUE)
    if (is.null(assigned_seurat@reductions$umap)) {
        stop(paste(
            "Could not find dimention reduction: UMAP",
            "Please reduct dimention with UMAP"
        ))
    }
    umap.path <- # nolint
        assigned_seurat@reductions$umap@misc$trajectory$linages # nolint
    if (!recalculate & !is.null(umap.path)) {
        warning(
            "Trajectory linages are already calculated.",
            "If you recalculate linage, designate recalculate = TRUE."
        )
        return(assigned_seurat)
    }
    message(
        "Start estimating trajetory"
    )
    target_clusters <-
        tibble::rownames_to_column(
            assigned_seurat@meta.data
        ) %>%
        dplyr::filter(
            celltype == root_celltype
        ) %>%
        dplyr::distinct(seurat_clusters)
    message(
        "Estimating trajectory pathway with slingshot..."
    )
    sds <-
        slingshot(
            Embeddings(assigned_seurat, "umap"),
            clusterLabels = assigned_seurat$seurat_clusters,
            start.clus = target_clusters$seurat_clusters,
            stretch = 0
        )
    message(
        "<<< Finished."
    )
    linages <- sds@metadata$lineages
    cnt.all <- 0
    for (linage in names(linages)) {
        pre <- -1
        cnt <- 0
        for (cluster in linages[[linage]]) {
            if (pre != -1) {
                df <- data.frame(
                    from = pre,
                    to = cluster
                )
                if (cnt == 0) {
                    all.df <- df
                } else {
                    all.df <- rbind(all.df, df)
                }
                cnt <- cnt + 1
            }
            pre <- cluster
        }
        if (cnt.all == 0) {
            path <- all.df
        } else {
            path <-
                rbind(path, all.df)
        }
        cnt.all <- cnt.all + 1
    }
    message(
        "<<< Finished calculating all lineages."
    )
    message(
        "<<< Extracting pathway data."
    )
    path <-
        path %>%
        dplyr::distinct(from, to) %>%
        dplyr::mutate(
            from = as.integer(from),
            to = as.integer(to)
        )
    umap.summarized <-
        tibble::rownames_to_column(
            as.data.frame(
                assigned_seurat@reductions$umap@cell.embeddings
            )
        ) %>%
        dplyr::rename(
            cell = rowname
        ) %>%
        dplyr::left_join(
            tibble::rownames_to_column(
                assigned_seurat@meta.data
            ) %>%
                dplyr::rename(
                    cell = rowname,
                    ident = seurat_clusters
                ) %>%
                dplyr::select(cell, ident),
            by = "cell"
        ) %>%
        dplyr::group_by(ident) %>%
        dplyr::summarise(
            UMAP_1 = mean(UMAP_1),
            UMAP_2 = mean(UMAP_2)
        )
    umap.path <-
        path %>%
        dplyr::left_join(
            umap.summarized,
            by = c("from" = "ident")
        ) %>%
        dplyr::rename(
            UMAP_x_from = UMAP_1,
            UMAP_y_from = UMAP_2
        ) %>%
        dplyr::left_join(
            umap.summarized,
            by = c("to" = "ident")
        ) %>%
        dplyr::rename(
            UMAP_x_to = UMAP_1,
            UMAP_y_to = UMAP_2
        )
    assigned_seurat@reductions$umap@misc$trajectory$pseudotime <-
        TrajectoryUtils::pathStats(sds)$pseudotime
    assigned_seurat@reductions$umap@misc$trajectory$linages <- umap.path
    message(
        "<<< Saved."
    )
    return(assigned_seurat)
}