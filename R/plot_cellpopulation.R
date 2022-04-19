plot_cellpopulation <- function(assigned_seurat, show_cellcount = TRUE, out_lambda = NULL) {
    library(ggplot2, quietly = TRUE)
    library(Seurat, quietly = TRUE)
    library(SeuratDisk, quietly = TRUE)
    library(dplyr, quietly = TRUE)
    library(patchwork, quietly = TRUE)
    ## p1
    p1 <-
        assigned_seurat@meta.data %>%
        group_by(celltype) %>%
        summarise(
            cellnum = n()
        ) %>%
        left_join(
            assigned_seurat@meta.data %>%
                tibble::rownames_to_column("cell") %>%
                distinct(celltype),
            by = "celltype"
        ) %>%
        ggplot(
            aes(
                x = celltype,
                y = cellnum,
                fill = celltype
            )
        ) +
        geom_bar(stat = "identity") +
        theme_classic() +
        scale_fill_manual(
            values = assigned_seurat@meta.data %>%
                arrange(celltype) %>%
                distinct(color.celltype) %>%
                pull()
        ) +
        theme(
            axis.text.x.bottom = element_blank(),
            axis.line.x.bottom = element_blank(),
            legend.position = "None",
            axis.ticks.x.bottom = element_blank(),
            axis.title.x.bottom = element_blank(),
            axis.line.y.left = element_blank(),
            axis.title.y = element_text(
                size = 14
            )
        ) +
        ylab("All Cell Count")
    ##

    meta <-
        assigned_seurat@meta.data %>%
        tibble::rownames_to_column("cell") %>%
        dplyr::select(
            cell, sample, seurat_clusters, celltype
        )

    all.cell.count <- # nolint
        meta %>%
        dplyr::group_by(sample) %>%
        dplyr::summarise(
            all = n(), .groups = "drop"
        )

    grouped <-
        meta %>%
        dplyr::group_by(
            sample, celltype
        ) %>%
        dplyr::summarise(
            count = n(), .groups = "drop"
        ) %>%
        dplyr::left_join(
            all.cell.count,
            by = "sample"
        ) %>%
        dplyr::mutate(
            percent = count / all * 100
        )
    sample.meta <-
        assigned_seurat@meta.data %>%
        dplyr::arrange(sample) %>%
        dplyr::distinct(color.sample)
    p2 <-
        grouped %>%
        dplyr::left_join(
            grouped %>%
                dplyr::group_by(celltype) %>%
                dplyr::summarise(
                    allcellcount = sum(percent)
                ),
            by = "celltype"
        ) %>%
        dplyr::mutate(
            percent = percent / allcellcount * 100
        ) %>%
        dplyr::select(
            sample, celltype, percent
        ) %>%
        ggplot(aes(celltype, percent)) +
        geom_bar(
            aes(fill = sample),
            position = "stack",
            stat = "identity"
        ) +
        scale_y_continuous(labels = scales::comma, expand = c(0.01, 0)) +
        coord_cartesian(clip = "off") +
        scale_fill_manual(name = "Sample", values = sample.meta$color.sample) +
        theme_bw() +
        theme(
            legend.position = "right",
            plot.title = element_text(hjust = 0.5),
            text = element_text(size = 16),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.title = element_blank(),
            plot.margin = margin(t = 20, r = 0, b = 0, l = 10, unit = "pt"),
            axis.text.x.bottom = element_text(
                angle = 45,
                hjust = 1,
                size = 15,
                colour = "black"
            ),
            axis.title.y.left = element_text(
                margin = margin(t = 0, r = 20, b = 0, l = 10, unit = "pt")
            ),
            axis.text.y.left = element_text(
                colour = "black"
            ),
            legend.justification = "bottom"
        ) +
        ylab("Normalized Sample Percent [%]")
    if (show_cellcount) {
        plt <-
            p1 + p2 +
            plot_layout(ncol = 1, heights = c(.2, .8)) +
            theme(
                plot.margin = unit(
                    c(-20, 0, 0, 0), "cm"
                )
            )
    } else{
        plt <- p2
    }
    if (is.null(out_lambda)) {
        return(plt)
    } else {
        out_lambda(plt)
        return(assigned_seurat)
    }
}