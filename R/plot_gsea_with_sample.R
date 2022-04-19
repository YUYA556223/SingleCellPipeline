plot_gsea_with_sample <- function(assigned_seurat, target.sample_name.1, target.sample_name.2, assay, out_lambda = NULL, target.celltype = c("all")) { # nolint
    library(clusterProfiler, quietly = TRUE)
    library(DOSE, quietly = TRUE)
    library(org.Hs.eg.db, quietly = TRUE)
    library(Seurat, quietly = TRUE)
    library(SeuratDisk, quietly = TRUE)
    library(dplyr, quietly = TRUE)
    library(tibble, quietly = TRUE)
    library(ggplot2, quietly = TRUE)
    if (target.celltype[1] == "all") {
        tmp <- assigned_seurat@meta.data %>%
            dplyr::distinct(celltype)
        target.celltype <- as.character(tmp$celltype)
    }
    cells.all <-
        assigned_seurat@meta.data %>%
        tibble::rownames_to_column("cell") %>%
        dplyr::select(
            cell,
            sample, seurat_clusters, celltype
        )
    cells.1 <- cells.all %>% dplyr::filter(sample == target.sample_name.1 &
        celltype %in% target.celltype)
    cells.2 <- cells.all %>% dplyr::filter(sample == target.sample_name.2 &
        celltype %in% target.celltype)
    SeuratObject::DefaultAssay(assigned_seurat) <- assay
    fc <- rownames_to_column(Seurat::FoldChange(
        assigned_seurat,
        cells.1$cell, cells.2$cell
    )) %>%
        dplyr::rename(gene = rowname) %>%
        dplyr::arrange(desc(avg_log2FC))
    entrez <- as.data.frame(mapIds(
        org.Hs.eg.db, fc$gene, "ENTREZID",
        "SYMBOL"
    ))
    colnames(entrez) <- c("entrez")
    gseadata <- cbind(entrez, fc) %>% dplyr::filter(isNA(entrez) ==
        FALSE)
    rownames(gseadata) <- gseadata$entrez
    gseadata <- gseadata %>%
        dplyr::select(entrez, avg_log2FC)
    genelist <- gseadata$avg_log2FC
    names(genelist) <- gseadata$entrez
    gse_result <- gseGO(
        geneList = genelist,
        OrgDb = org.Hs.eg.db,
        ont = "ALL", # "BP","CC","MF","ALL"から選択
        verbose = TRUE
    )
    rp <-
        ridgeplot(
            gse_result
        ) +
        ggplot2::theme(
            legend.justification = "bottom",
            plot.title = element_text(
                hjust = 0.5,
                face = "bold",
                size = 20
            )
        ) +
        ggplot2::scale_colour_viridis_c(
            alpha = 1,
            begin = 0,
            end = 1,
            direction = 1,
            option = "C",
            aesthetics = "fill"
        ) +
        ggplot2::ggtitle(
            paste(
                "GSEA on ",
                target.sample_name.1,
                " (vs.",
                target.sample_name.2,
                ") ",
                sep = ""
            )
        ) +
        xlab(
            paste(
                "log2FoldChange\n\n",
                target.sample_name.2,
                "  <--                            -->  ",
                target.sample_name.1,
                sep = ""
            )
        ) +
        ylab(
            "GO Terms"
        )
    if (is.null(out_lambda)) {
        return(rp)
    } else {
        out_lambda(rp)
        return(assigned_seurat)
    }
}