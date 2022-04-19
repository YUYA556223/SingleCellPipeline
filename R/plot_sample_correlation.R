plot_sample_correlation <- function(assigned_seurat, assay, sample.name.1, sample.name.2, out_lambda = NULL) { # nolint
    library(viridis)
    marker.1 <- # nolint
        assigned_seurat@assays[[assay]]@misc$sample.marker[[sample.name.1]]
    marker.2 <- # nolint
        assigned_seurat@assays[[assay]]@misc$sample.marker[[sample.name.2]]
    if (is.null(marker.1)) {
        message(
            paste(
                "Could not find marker1 info.",
                "Start finding marker. It may take some time.",
                sep = "\n"
            )
        )
        markers <-
            assigned_seurat %>%
            subset(
                subset = sample == sample.name.1
            ) %>%
            Seurat::FindAllMarkers(
                assay = assay,
                verbose = TRUE,
                densify = TRUE
            )
        assigned_seurat@assays$SCT@misc$sample.marker[[sample.name.1]] <-
            markers
        marker.1 <- markers # nolint
    }
    if (is.null(marker.2)) {
        message(
            paste(
                "Could not find marker2 info.",
                "Start finding marker. It may take some time.",
                sep = "\n"
            )
        )
        markers <-
            assigned_seurat %>%
            subset(
                subset = sample == sample.name.2
            ) %>%
            Seurat::FindAllMarkers(
                assay = assay,
                verbose = TRUE,
                densify = TRUE
            )
        assigned_seurat@assays$SCT@misc$sample.marker[[sample.name.2]] <-
            markers
        marker.2 <- markers # nolint
    }
    if (!is.null(out_lambda)) {
        out_lambda(
            list(
                marker.1 = marker.1,
                marker.2 = marker.2
            )
        )
    }
    message("Calculating correlation...")
    celltypes <-
        assigned_seurat %>%
        SingleCellPipeline::select_cell.meta(out_lambda = NULL) %>%
        dplyr::arrange(celltype) %>%
        dplyr::distinct(celltype)
    celltypes <- as.character(celltypes$celltype)
    combinations <-
        as.data.frame(
            gtools::permutations(
                n = length(celltypes), r = 2, v = celltypes,
                repeats.allowed = TRUE
            )
        ) %>%
        dplyr::rename(
            primary = V1,
            secondary = V2
        ) %>%
        dplyr::mutate(
            num = row_number()
        )
    for (cnt in combinations$num) {
        target <-
            combinations %>%
            dplyr::filter(
                num == cnt
            )
        markers <-
            marker.1 %>%
            dplyr::filter(cluster == target$primary) %>%
            tibble::rownames_to_column() %>%
            dplyr::select(gene, avg_log2FC) %>%
            dplyr::rename(
                primary = avg_log2FC
            ) %>%
            dplyr::full_join(
                marker.2 %>%
                    dplyr::filter(cluster == target$secondary) %>%
                    tibble::rownames_to_column() %>%
                    dplyr::select(gene, avg_log2FC) %>%
                    dplyr::rename(
                        secondary = avg_log2FC
                    ),
                by = "gene"
            ) %>%
            dplyr::mutate_all(
                ~ replace(., is.na(.), 0)
            )
        corelation <-
            target %>%
            dplyr::mutate(
                cor(
                    x = markers %>%
                        dplyr::select(primary) %>%
                        dplyr::pull(),
                    y = markers %>%
                        dplyr::select(secondary) %>%
                        dplyr::pull(),
                    method = "pearson"
                )
            )
        if (cnt == 1) {
            results <- corelation
        } else {
            results <- rbind(results, corelation)
        }
    }
    results <-
        results %>%
        dplyr::select(-num)
    colnames(results) <-
        c("primary", "secondary", "cor")
    cell_level <-
        assigned_seurat %>%
        SingleCellPipeline::select_cell.meta(out_lambda = NULL) %>%
        dplyr::distinct(celltype) %>%
        dplyr::arrange(celltype) %>%
        dplyr::select(celltype) %>%
        dplyr::mutate(
            celltype = as.character(celltype)
        ) %>%
        dplyr::pull()
    results <-
        results %>%
        dplyr::mutate(
            primary = factor(
                primary,
                levels = cell_level
            ),
            secondary = factor(
                secondary,
                levels = cell_level
            )
        )
    plt <-
        results %>%
        ggplot2::ggplot() +
        ggplot2::geom_tile(
            ggplot2::aes(
                x = primary,
                y = secondary,
                fill = cor
            )
        ) +
        scale_fill_viridis_c(option = "C") +
        ggplot2::xlab(sample.name.1) +
        ggplot2::ylab(sample.name.2) +
        ggplot2::scale_fill_gradientn(
            colors = c("#002897", "#F2EA91", "#cf0d24"),
            limits = c(-1, 1)
        )
    return(plt)
}