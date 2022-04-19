# Integrate multiple sample data into seuratobject.
# [WARNING!] This process requires plenty of time and memory.
# @
# > project <- "vascular"
# > basepath <-
# >      "/home/yuyasato/work2/_organoid_scRNA/raw/matrix/__filtered/__fetal"
# > samples <- list(
# >     fetalPCW14 = list(
# >         path = paste(basepath, "fetalPCW14_0rsydy7.tsv", sep = "/")
# >     ),
# >     fetalPCW15 = list(
# >         path = paste(basepath, "fetalPCW15_0rsydy7.tsv", sep = "/")
# >     ),
# >     fetalPCW16 = list(
# >         path = paste(basepath, "fetalPCW16_GSE162170.tsv", sep = "/")
# >     ),
# >     fetalPCW19 = list(
# >         path = paste(basepath, "fetalPCW19_0rsydy7.tsv", sep = "/")
# >     ),
# >     fetalPCW20 = list(
# >         path = paste(basepath, "fetalPCW20_GSE162170.tsv", sep = "/")
# >     ),
# >     fetalPCW21 = list(
# >         path = paste(basepath, "fetalPCW21_GSE162170.tsv", sep = "/")
# >     ),
# >     fetalPCW24 = list(
# >         path = paste(basepath, "fetalPCW24_GSE162170.tsv", sep = "/")
# >     )
# > )
# > resolution <- 2.0
# > markers <- "/home/yuyasato/work2/_organoid_scRNA/ref/__markers_refined.tsv"


integrate <-
    function(samples, project_name, low_feature.cutoff, outputdir, resolution, marker_path, useSCT) { # nolint
        library(Seurat)
        library(SeuratDisk)
        library(SeuratData)
        library(ggplot2)
        library(dplyr)
        output <- list(
            clustered = paste(outputdir, "seurat", sep = "/"),
            markers = paste(outputdir, "markers.tsv", sep = "/"),
            clusteringGraph = paste(outputdir, "cluster.png", sep = "/"),
            cluster_each = paste(outputdir, "sample.png", sep = "/"),
            graphfeature = paste(outputdir, "feature.png", sep = "/"),
            graphfeatureclusters = paste(
                outputdir, "cluster_feature.png",
                sep = "/"
            ),
            merged = paste(outputdir, "merged.png", sep = "/"),
            feature_dot = paste(outputdir, "dotplot.png", sep = "/")
        )
        seurats <- list()
        gc()
        if (!dir.exists(outputdir)) {
            dir.create(outputdir)
        }

        for (sample in names(samples)) {
            message(paste(
                "start loading ",
                sample,
                sep = " "
            ))
            seurats[[sample]] <- read.table(
                samples[[sample]]$path,
                sep = "\t", row.names = 1, header = TRUE
            )
            message(paste(
                "finished loading ",
                sample,
                sep = " "
            ))
            seurats[[sample]] <-
                Seurat::CreateSeuratObject(
                    counts = seurats[[sample]],
                    project = project_name, min.cells = 5
                )
            seurats[[sample]]@meta.data$sample <- sample
            seurats[[sample]] <-
                Seurat::PercentageFeatureSet(
                    seurats[[sample]],
                    pattern = "^MT-",
                    col.name = "percent.mt"
                )
            seurats[[sample]] <- subset(
                seurats[[sample]],
                subset = nFeature_RNA > low_feature.cutoff
            )
            seurats[[sample]] <-
                Seurat::SCTransform(
                    seurats[[sample]],
                    vars.to.regress = "percent.mt", verbose = TRUE
                )
        }

        features <- SelectIntegrationFeatures(
            object.list = seurats,
            verbose = TRUE
        )
        seurats <-
            PrepSCTIntegration(
                object.list = seurats,
                anchor.features = features,
                verbose = TRUE
            )
        # seurats <- lapply(
        #     X = seurats,
        #     FUN = RunPCA,
        #     features = features
        # )
        # Find Anchors
        anchors <- Seurat::FindIntegrationAnchors(
            object.list = seurats,
            anchor.features = features,
            normalization.method = "SCT",
            verbose = TRUE
        )
        rm(seurats)
        rm(features)
        merged <- Seurat::IntegrateData(
            anchorset = anchors,
            normalization.method = "SCT",
            new.assay.name = "integrated",
            verbose = TRUE
        )
        rm(anchors)

        SeuratObject::DefaultAssay(merged) <- "integrated"
        merged <- Seurat::RunPCA(
            merged,
            npcs = 30,
            verbose = TRUE
        )
        merged <- Seurat::RunUMAP(
            merged,
            reduction = "pca",
            dims = 1:30,
            verbose = TRUE
        )

        merged <- Seurat::FindNeighbors(
            merged,
            reduction = "umap",
            dims = 1:2,
            verbose = TRUE
        )
        merged <- Seurat::FindClusters(
            merged,
            reduction.type = "umap",
            algorithm = 1,
            resolution = resolution,
            verbose = TRUE
        )

        p1 <-
            Seurat::DimPlot(
                merged,
                reduction = "umap",
                group.by = "sample"
            )
        p2 <- Seurat::DimPlot(
            merged,
            reduction = "umap",
            label = TRUE,
            repel = TRUE
        ) + NoLegend()

        ggplot2::ggsave(output$clusteringGraph, p1 + p2,
            width = 40, height = 20, units = "cm", dpi = 400
        )

        p4 <- Seurat::DimPlot(merged, reduction = "umap", split.by = "sample")
        ggplot2::ggsave(output$cluster_each, p4,
            width = 80, height = 20, units = "cm", dpi = 400
        )

        plt <-
            Seurat::DimPlot(
                merged,
                reduction = "umap", group.by = "sample"
            ) +
            ggplot2::theme(
                legend.justification = "bottom"
            ) +
            ggplot2::ggtitle("")
        ggplot2::ggsave(
            output$merged,
            plt,
            width = 25, height = 20, units = "cm", dpi = 400
        )



        p5 <- VlnPlot(
            merged,
            features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
            ncol = 3, group.by = "sample"
        )
        ggplot2::ggsave(output$graphfeature, p5,
            width = 40, height = 40, units = "cm", dpi = 400
        )
        p6 <- VlnPlot(
            merged,
            features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
            ncol = 3, group.by = "seurat_clusters"
        )

        ggplot2::ggsave(output$graphfeatureclusters, p6,
            width = 80, height = 40, units = "cm", dpi = 400
        )

        marker <-
            read.table(
                marker_path,
                sep = "\t",
                header = TRUE
            ) %>%
            dplyr::select(gene)
        color <- "white"
        dp <- Seurat::DotPlot(
            merged,
            features = marker$gene,
            scale.min = 0,
            dot.scale = 13,
            cols = c("#ffffff00", "#d60000"),
            col.min = 0,
            assay = "SCT"
        ) +
            theme(
                rect = element_rect(colour = color, fill = color),
                text = element_text(color = "#000000"),
                plot.background = element_rect(colour = color, fill = color),
                panel.background = element_rect(colour = color, fill = color),
                panel.border =
                    element_rect(fill = NA, colour = "#1b1b1b", size = 1),
                panel.grid.major = element_line(colour = "grey60"),
                panel.grid.minor = element_line(colour = "grey30"),
                legend.key = element_rect(colour = color, fill = color),
                axis.text = element_text(colour = "#000000"),
                axis.text.x.bottom = element_text(
                    size = 18,
                    hjust = .5, vjust = .5
                ),
                axis.text.y.left = element_text(
                    size = 18
                ),
                axis.title.y.left = element_text(
                    size = 20,
                    margin = unit(c(0, 1, 0, 1), "lines")
                ),
                axis.title.x.bottom = element_text(
                    size = 20,
                    margin = unit(c(1, 0, 1, 0), "lines")
                )
            ) +
            xlab("Gene Features") +
            ylab("Cluster") + coord_flip()

        ggplot2::ggsave(
            output$feature_dot,
            dp,
            width = 90, height = 50, units = "cm", dpi = 200
        )

        SeuratDisk::SaveH5Seurat(
            object = merged,
            filename = output$clustered,
            overwrite = TRUE
        )


        markers <- Seurat::FindAllMarkers(
            merged,
            only.pos = TRUE,
            verbose = TRUE
        )
        write.table(
            markers,
            file = output$markers,
            sep = "\t"
        )
    }