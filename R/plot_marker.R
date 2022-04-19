plot_marker <- function(seurat, marker.genes, assay, out_lambda = NULL) { # nolint
    library(dplyr, quietly = TRUE)
    library(ggplot2, quietly = TRUE)
    plt <-
        seurat %>%
        Seurat::DotPlot(
            features = marker.genes, scale.min = 0,
            dot.scale = 13, cols = c("#ffffff00", "#d60000"), col.min = 0,
            assay = assay
        ) + theme(
            rect = element_rect(
                colour = "white",
                fill = "white"
            ), text = element_text(color = "#000000"),
            plot.background = element_rect(colour = "white", fill = "white"),
            panel.background = element_rect(colour = "white", fill = "white"),
            panel.border = element_rect(
                fill = NA, colour = "#1b1b1b",
                size = 1
            ), panel.grid.major = element_line(colour = "grey60"),
            panel.grid.minor = element_line(
                colour = "grey30"
            ), legend.key = element_rect(
                colour = "white",
                fill = "white"
            ), axis.text = element_text(colour = "#000000"),
            axis.text.x.bottom = element_text(
                size = 18, hjust = 0.5,
                vjust = 0.5
            ), axis.text.y.left = element_text(size = 18),
            axis.title.y.left = element_text(size = 20, margin = unit(c(
                0,
                1, 0, 1
            ), "lines")), axis.title.x.bottom = element_text(
                size = 20,
                margin = unit(c(1, 0, 1, 0), "lines")
            )
        ) + xlab("Gene Features") +
        ylab("Cluster") + coord_flip()
    if (is.null(out_lambda)) {
        return(plt)
    }
    out_lambda(plt)
    return(seurat)
}