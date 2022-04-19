plot_samples_with_celltype <- function(assigned_seurat, disable.color = "gray", out_lambda = NULL, raster = FALSE) { # nolint
    library(Seurat)
    library(SeuratDisk)
    library(ggplot2)
    library(dplyr)
    library(tibble)
    library(crayon)
    if (!SingleCellPipeline::is_assigned(assigned_seurat)) {
        stop(
            paste(
                "No cell type has been assigned.",
                "Run 'assign_celltype' function first."
            )
        )
        return(NULL)
    }
    samples_all <-
        assigned_seurat@meta.data %>%
        dplyr::arrange(sample) %>%
        dplyr::distinct(sample)
    cnt <- 0
    for (samplename in assigned_seurat@meta.data %>%
        distinct(sample) %>%
        pull() %>%
        levels()) {
        cat(green(
            "Extracting " %+% samplename %+% " \n"
        ))
        target.assigned_seurat <- assigned_seurat
        target.assigned_seurat@meta.data <-
            target.assigned_seurat@meta.data %>%
            dplyr::select(sample, celltype, color.sample, color.celltype) %>%
            dplyr::mutate(
                celltype = if_else(
                    sample == samplename,
                    true = as.character(celltype),
                    false = "Disabled"
                ),
                sample = factor(
                    sample,
                    levels = samples_all %>%
                        arrange(sample) %>%
                        mutate(
                            factor.data = if_else(
                                sample == samplename,
                                true = 1,
                                false = 0
                            ),
                            sample = as.character(sample)
                        ) %>%
                        arrange(factor.data) %>%
                        distinct(sample) %>%
                        pull()
                ),
                color.celltype = if_else(
                    sample == samplename,
                    true = as.character(color.celltype),
                    false = disable.color
                )
            )
        plt <-
            target.assigned_seurat %>%
            Seurat::DimPlot(raster = raster) +
            ggplot2::ggtitle(samplename)

        plt$data <-
            plt$data %>%
            tibble::rownames_to_column("cell") %>%
            left_join(
                target.assigned_seurat@meta.data %>%
                    tibble::rownames_to_column("cell") %>%
                    select(cell, celltype, color.celltype),
                by = "cell"
            ) %>%
            mutate(
                ident = factor(
                    celltype,
                    levels = c(assigned_seurat@meta.data %>%
                        distinct(celltype) %>%
                        pull() %>%
                        levels(), "Disabled")
                )
            ) %>%
            select(-celltype) %>%
            tibble::column_to_rownames("cell") %>%
            arrange(ident)
        plt <-
            plt +
            scale_color_manual(
                values = plt$data %>%
                    arrange(ident) %>%
                    distinct(color.celltype) %>%
                    pull()
            ) + Seurat::NoLegend()
        plt$data <-
            plt$data %>%
            arrange(desc(ident))
        if (cnt == 0) {
            plt.all <- plt
        } else {
            plt.all <-
                plt.all + plt
        }
        rm(target.assigned_seurat)
        cnt <- cnt + 1
    }
    if (is.null(out_lambda)) {
        return(plt.all)
    } else {
        out_lambda(plt.all)
    }
    return(assigned_seurat)
}