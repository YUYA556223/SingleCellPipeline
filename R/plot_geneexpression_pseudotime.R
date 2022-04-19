plot_geneexpression_pseudotime <- function(assigned_seurat, feature_gene, assay, target_cells = c("all"), excl_samples = c(), scale.free = TRUE, facet.wrap = TRUE, layout.nrows = NULL, layout.ncols = NULL, out_lambda = NULL) { # nolint
    dp <-
        assigned_seurat %>%
        Seurat::DoHeatmap(
            features = feature_gene,
            assay = assay, slot = "data"
        )
    message("@ Start executing plot_geneexpression_pseudotime >>>")
    if (!SingleCellPipeline::is_assigned(assigned_seurat)) {
        stop(
            paste(
                "No cell type has been assigned_seurat.",
                "Run 'assign_celltype' function first."
            )
        )
        return(NULL)
    }
    message("<<< Extracting")
    if (target_cells[1] == "all") {
        # assign all cell types
        target <-
            assigned_seurat@meta.data %>%
            dplyr::distinct(celltype)
        target_cells <- target$celltype
        rm(target)
    }
    message("<<< Running pseudotime calculation")
    pseudotime <-
        as.data.frame(
            assigned_seurat@reductions$umap@misc$trajectory$pseudotime
        ) %>%
        tibble::rownames_to_column("cell")
    pseudotime.summarized <- # nolint
        as.data.frame(pseudotime) %>%
        dplyr::mutate(
            pseudotime = dplyr::select(
                .,
                dplyr::starts_with("Lineage")
            ) %>%
                rowMeans(na.rm = TRUE)
        ) %>%
        dplyr::select(cell, pseudotime)
    all.samples <-
        assigned_seurat@meta.data %>%
        dplyr::distinct(sample) %>%
        dplyr::filter(!sample %in% excl_samples)
    message("<<< Extracting samples")
    all.samples <- all.samples$sample
    feature <-
        dp$data %>%
        dplyr::filter(
            !is.na(Expression)
        )
    feature <-
        feature %>%
        dplyr::left_join(
            pseudotime.summarized,
            by = c("Cell" = "cell")
        ) %>%
        dplyr::left_join(
            assigned_seurat@meta.data %>%
                tibble::rownames_to_column("cell") %>%
                dplyr::select(
                    cell, sample
                ),
            by = c("Cell" = "cell")
        )
    size <- feature %>%
        dplyr::filter(Expression != 0)
    sample.color <- # nolint
        assigned_seurat@meta.data %>%
        tibble::rownames_to_column("cell") %>%
        dplyr::distinct(sample, color.sample) %>%
        dplyr::mutate(
            fill = "#a3a3a3"
        )
    if (scale.free) {
        scale.param <- "free"
    } else {
        scale.param <- "fixed"
    }
    excl_sample_geneset <-
        feature %>%
        group_by(Feature, sample) %>%
        summarise(
            mean = mean(Expression),
            .groups = "drop"
        ) %>%
        dplyr::filter(
            mean == 0
        )
    plt <-
        feature %>%
        dplyr::filter(
            !(Feature %in% excl_sample_geneset$Feature &
                sample %in% excl_sample_geneset$sample)
        ) %>%
        ggplot2::ggplot(
            aes(
                x = pseudotime,
                y = Expression,
                color = sample,
                fill = sample
            )
        ) +
        scale_color_manual(
            values = sample.color$color.sample
        ) +
        scale_fill_manual(
            values = sample.color$fill
        ) +
        stat_smooth(
            method = "gam",
            formula = y ~ s(x, bs = "cs"),
            na.rm = TRUE
        ) +
        theme_minimal() +
        theme(
            legend.justification = "bottom",
            legend.text.align = 0,
            panel.background = element_rect(
                fill = "white"
            ),
            plot.background = element_rect(
                fill = "white",
                color = "white",
                size = 0
            ),
            text = element_text(
                colour = "black"
            ),
            strip.text.x = element_text(
                face = "italic"
            ),
            panel.border = element_blank()
        ) +
        xlim(0, NA) +
        ylim(0, NA) +
        xlab("Pseudotime") +
        ylab("Expression")
    if (facet.wrap) {
        plt <-
            plt +
            facet_wrap(
                . ~ Feature,
                scales = scale.param,
                nrow = layout.nrows,
                ncol = layout.ncols
            )
    } else {
        plt <-
            plt +
            facet_grid(
                . ~ Feature,
                scales = scale.param
            )
    }
    plt$data <-
        plt$data %>%
        dplyr::mutate(
            Feature = factor(
                Feature,
                levels = feature_gene
            )
        ) %>%
        dplyr::arrange(Feature)
    if (is.null(out_lambda)) {
        return(plt)
    }
    out_lambda(plt)
    return(assigned_seurat)
}