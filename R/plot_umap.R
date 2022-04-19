plot_umap <- function(assigned_seurat, out_lambda = NULL, plot_trajectory_path = FALSE, raster = TRUE) { # nolint
  library(ggplot2)
  library(Seurat)
  library(SeuratDisk)
  library(dplyr)
  if (!is_assigned(assigned_seurat)) {
    stop(
      paste(
        "No cell type has been assigned.",
        "Run 'assign_celltype' function first."
      )
    )
    return(NULL)
  }

  celltype_color <-
    tibble::rownames_to_column(
      assigned_seurat@meta.data
    ) %>%
    dplyr::arrange(celltype) %>%
    dplyr::distinct(color.celltype)
  celltype <-
    tibble::rownames_to_column(
      assigned_seurat@meta.data
    ) %>%
    dplyr::arrange(celltype) %>%
    dplyr::distinct(celltype)
  p <-
    Seurat::DimPlot(
      assigned_seurat,
      reduction = "umap", label = FALSE, repel = TRUE,
      label.color = "#000000", raster = raster
    ) +
    scale_color_manual("Cell type", values = celltype_color$color.celltype) +
    ggplot2::theme(
      legend.justification = "bottom"
    )
  p$data$ident <-
    factor(
      p$data$ident,
      levels = levels(
        celltype$celltype
      )
    )
  if (plot_trajectory_path) {
    arrow_style <- arrow(
      angle = 20,
      ends = "last",
      type = "closed",
      length = grid::unit(5, "pt")
    )
    p <-
      p +
      geom_segment(
        data = SingleCellPipeline::get_trajectory_linages(
          assigned_seurat,
          out_lambda = NULL
        ),
        aes(
          x = UMAP_x_from, y = UMAP_y_from,
          xend = UMAP_x_to, yend = UMAP_y_to
        ),
        colour = "#161616", arrow = arrow_style, size = .5
      )
  }
  if (is.null(out_lambda)) {
    return(p)
  }
  out_lambda(p)
  return(assigned_seurat)
}