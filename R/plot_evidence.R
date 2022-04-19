plot_evidence <- function(assigned_seurat, as.seurat_cluster = TRUE, assay, display_summary = FALSE, out_lambda = NULL) { # nolint
    library(ggplot2, quietly = TRUE)
    library(Seurat, quietly = TRUE)
    library(SeuratDisk, quietly = TRUE)
    library(dplyr, quietly = TRUE)
    if (!is_assigned(assigned_seurat)) {
        stop(
            paste(
                "No cell type has been assigned_seurat.",
                "Run 'assign_celltype' function first."
            )
        )
        return(NULL)
    }
    features <-
        assigned_seurat@meta.data %>%
        dplyr::arrange(celltype) %>%
        dplyr::distinct(evidence)
    cnt <- 0
    for (gene in features$evidence) {
        for (splitted in strsplit(gene, ", ")) {
            target.colorset <- # nolint
                as.data.frame(splitted) %>%
                dplyr::mutate(
                    nonsplitted = gene
                ) %>%
                dplyr::rename(
                    gene = splitted
                ) %>%
                dplyr::left_join(
                    assigned_seurat@meta.data %>%
                        dplyr::select(evidence, celltype),
                    by = c("nonsplitted" = "evidence")
                )
            if (cnt == 0) {
                all.features <- target.colorset
            } else {
                all.features <-
                    rbind(all.features, target.colorset)
            }
        }
        cnt <- cnt + 1
    }
    all.features <-
        all.features %>%
        dplyr::left_join(
            assigned_seurat@meta.data %>%
                dplyr::select(celltype, color.celltype),
            by = "celltype"
        ) %>%
        dplyr::distinct(gene, celltype, color.celltype) %>%
        dplyr::rename(
            color = color.celltype
        )

    nondup.gene <-
        all.features %>%
        group_by(gene) %>%
        dplyr::filter(
            n() <= 1
        ) %>%
        dplyr::select(gene, color)
    gene.dup <-
        all.features %>%
        group_by(gene) %>%
        dplyr::filter(
            n() > 1
        ) %>%
        dplyr::distinct(gene) %>%
        dplyr::mutate(
            color = "#666666"
        )
    merged <-
        rbind(gene.dup, nondup.gene)
    message("Constructing vln plot")
    warning("These genes are dupulicated")
    message(gene.dup$gene)
    if (as.seurat_cluster) {
        idents <-
            as.data.frame(SeuratObject::Idents(assigned_seurat))
        colnames(idents) <- c("idents")
        cluster_idents <-
            tibble::rownames_to_column(
                idents
            ) %>%
            dplyr::left_join(
                tibble::rownames_to_column(assigned_seurat@meta.data) %>%
                    dplyr::select(rowname, seurat_clusters, celltype),
                by = "rowname"
            ) %>%
            dplyr::rename(
                cell = rowname
            )
        ordered_seurat_clusters <-
            cluster_idents %>%
            dplyr::arrange(celltype) %>%
            dplyr::distinct(seurat_clusters)
        cluster <-
            cluster_idents %>%
            dplyr::select(cell, seurat_clusters) %>%
            dplyr::mutate(
                idents = factor(
                    seurat_clusters,
                    levels = ordered_seurat_clusters$seurat_clusters
                )
            ) %>%
            dplyr::select(-seurat_clusters)
        rownames(cluster) <- cluster$cell
        cluster <-
            cluster %>%
            dplyr::select(-cell)
        SeuratObject::Idents(assigned_seurat) <-
            as.list(cluster)$idents
    }
    vlnplot <-
        vlnplot <-
        VlnPlot(
            assigned_seurat,
            merged$gene,
            stack = TRUE, flip = TRUE,
            assay = assay
        ) +
        theme(
            legend.position = "none",
            panel.background = element_rect(
                fill = "white",
                color = "white"
            ),
            plot.background = element_rect(
                fill = "white",
                color = "white"
            ),
            strip.text = element_text(
                face = "italic",
                hjust = 0,
                vjust = 0.5
            )
        ) +
        scale_fill_manual(
            values = merged$color
        )

    if (display_summary) {
        vlnplot <-
            vlnplot +
            stat_summary(
                fun = mean, geom = "line", aes(group = 1)
            )
    }
    if (as.seurat_cluster) {
        cluster.info <- # nolint
            assigned_seurat@meta.data %>%
            dplyr::arrange(celltype) %>%
            dplyr::distinct(seurat_clusters)
        vlnplot$data <-
            vlnplot$data %>%
            dplyr::mutate(
                ident = factor(
                    ident,
                    cluster.info$seurat_clusters
                )
            )
    } else {
        cluster.info <- # nolint
            assigned_seurat@meta.data %>%
            dplyr::arrange(celltype) %>%
            dplyr::distinct(celltype)
        vlnplot$data <-
            vlnplot$data %>%
            dplyr::mutate(
                ident = factor(
                    ident,
                    cluster.info$celltype
                )
            )
    }
    message("Finalizing...")
    if (is.null(out_lambda)) {
        stop(
            paste(
                "out_lambda function is null.",
                "designate out_lambda function like",
                "SingleCellPipeline::lambda(plot: ",
                "ggplot2::ggsave('filename', plot))",
                ")"
            )
        )
        return(assigned_seurat)
    } else {
        out_lambda(vlnplot)
    }
    return(assigned_seurat)
}