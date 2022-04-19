plot_feature_backSPIN <- function(assigned_seurat, cell_order, assay, target.celltypes, target.features, ctrl, out_lambda = NULL) {
    library(Seurat)
    library(SeuratDisk)
    library(SingleCellPipeline)
    library(tidyverse)
    library(patchwork)
    target_seurat <-
        assigned_seurat %>%
        subset(
            subset = celltype %in% target.celltypes
        )
    target_seurat@meta.data <-
        target_seurat@meta.data %>%
        mutate(
            cellnum = cell_order
        )
    data <- data.frame()
    for (target.feature in target.features) {
        plt <-
            target_seurat %>%
            subset(
                celltype %in% target.celltypes
            ) %>%
            set_default_assay(assay = assay) %>%
            FeaturePlot(
                features = target.feature
            )
        data <-
            data %>%
            rbind(
                plt$data %>%
                    tibble::rownames_to_column("cell") %>%
                    rename(
                        expression = dplyr::starts_with(target.feature)
                    ) %>%
                    mutate(
                        feature = target.feature,
                        expression = as.double(expression)
                    ) %>%
                    select(cell, expression, feature) %>%
                    left_join(
                        target_seurat@meta.data %>%
                            tibble::rownames_to_column("cell") %>%
                            select(cell, cellnum, sample),
                        by = "cell"
                    )
            )
    }
    p1 <-
        data %>%
        mutate(
            feature = factor(
                feature,
                levels = target.features
            )
        ) %>%
        ggplot(
            aes(
                x = cellnum,
                y = expression,
                fill = feature
            )
        ) +
        geom_bar(stat = "identity") +
        facet_grid(. ~ feature) +
        theme(
            axis.text.x.bottom = element_blank(),
            axis.ticks.x.bottom = element_blank(),
            axis.title.x.bottom = element_blank(),
            panel.background = element_rect(
                colour = "#a8a8a8",
                fill = "white"
            ),
            panel.grid = element_line(
                colour = "#e7e7e7"
            ),
            legend.position = "none",
            strip.text.x = element_text(
                face = "italic"
            )
        ) +
        scale_y_continuous(expand = c(0, 0)) +
        ylab("")
    data <- data.frame()
    samples <- target_seurat@meta.data %>%
        distinct(sample) %>%
        pull()
    for (target.feature in target.features) {
        plt <-
            target_seurat %>%
            subset(
                celltype %in% target.celltypes
            ) %>%
            set_default_assay(assay = assay) %>%
            FeaturePlot(
                features = target.feature
            )
        for (target.sample in samples) {
            data <-
                data %>%
                rbind(
                    plt$data %>%
                        tibble::rownames_to_column("cell") %>%
                        rename(
                            expression = dplyr::starts_with(target.feature)
                        ) %>%
                        mutate(
                            feature = target.feature,
                            expression = as.double(expression),
                            target_sample = target.sample
                        ) %>%
                        select(cell, expression, feature, target_sample) %>%
                        left_join(
                            target_seurat@meta.data %>%
                                tibble::rownames_to_column("cell") %>%
                                select(cell, cellnum, sample),
                            by = "cell"
                        ) %>%
                        mutate(
                            color = if_else(
                                target_sample == sample,
                                true = as.character(sample),
                                false = "others"
                            )
                        )
                )
        }
    }
    target_seurat_ctrl <-
        target_seurat %>%
        subset(
            subset = sample == ctrl
        )
    plt <-
        target_seurat_ctrl %>%
        set_default_assay(assay = assay) %>%
        DoHeatmap(
            features = target.features,
            raster = FALSE,
            slot = "data"
        )
    p2 <-
        plt$data %>%
        dplyr::filter(
            !is.na(Expression)
        ) %>%
        dplyr::left_join(
            target_seurat_ctrl@meta.data %>%
                tibble::rownames_to_column("Cell") %>%
                arrange(cellnum) %>%
                mutate(
                    cellnum = row_number()
                ) %>%
                select(Cell, cellnum),
            by = "Cell"
        ) %>%
        mutate(
            sample = ctrl
        ) %>%
        arrange(cellnum) %>%
        ggplot(
            aes(
                x = cellnum,
                y = Expression,
                fill = Feature
            )
        ) +
        geom_bar(stat = "identity") +
        facet_grid(. ~ Feature) +
        theme(
            axis.text.x.bottom = element_blank(),
            axis.ticks.x.bottom = element_blank(),
            legend.position = "none",
            plot.background = element_rect(
                fill = "white",
                colour = "white"
            ),
            panel.background = element_rect(
                fill = "white", color = "#585858"
            ),
            panel.grid = element_line(colour = "#dadada"),
            strip.background.x = element_blank(),
            strip.text.x = element_blank(),
            axis.title.x.bottom = element_blank()
        ) +
        scale_y_continuous(
            expand = c(0, 0)
        ) +
        labs(
            y = "Coltrol"
        )
    p3 <-
        data %>%
        mutate(
            color = factor(
                color,
                levels = c(
                    data %>%
                        arrange(sample) %>%
                        distinct(sample) %>%
                        pull() %>%
                        as.character(),
                    "others"
                )
            ),
            feature = factor(
                feature,
                levels = target.features
            ),
            target_sample = factor(
                target_sample,
                levels = data %>%
                    distinct(sample) %>%
                    pull() %>%
                    levels()
            )
        ) %>%
        ggplot(
            aes(
                x = cellnum,
                y = expression,
                fill = color
            )
        ) +
        geom_bar(stat = "identity") +
        facet_grid(
            target_sample ~ feature
        ) +
        theme(
            legend.justification = "bottom",
            axis.ticks.x.bottom = element_blank(),
            axis.text.x.bottom = element_blank(),
            panel.background = element_rect(
                fill = "white",
                color = "#585858"
            ),
            panel.grid = element_line(
                colour = "#dadada"
            ),
            strip.text.x = element_blank(),
            strip.background.x = element_blank(),
            strip.text.y = element_text(
                size = 7
            )
        ) +
        scale_fill_manual(
            values = c(
                target_seurat@meta.data %>%
                    arrange(sample) %>%
                    distinct(color.sample) %>%
                    pull(),
                "#e7e7e7"
            )
        ) +
        labs(
            x = "BackSPINed Cell Order",
            y = "Feature Expression",
            fill = "Sample"
        )
    rm(data)
    plt <-
        p1 + p2 + p3 +
        plot_layout(
            ncol = 1, heights = c(.05, .05, 1)
        ) +
        theme(
            plot.margin = unit(c(0, 0, 0, 0), "cm"),
            legend.key.size = unit(.3, "cm")
        )
    if (is.null(out_lambda)) {
        return(plt)
    } else {
        out_lambda(plt)
    }
}