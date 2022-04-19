plot_features <- function(assigned_seurat, target.features, cell_order, assay, title = "", out_lambda = NULL) {
    cnt <- 0
    all.samples <-
        assigned_seurat@meta.data %>%
        arrange(sample) %>%
        distinct(sample) %>%
        pull()
    for (target.sample in all.samples) {
        orders <-
            assigned_seurat@meta.data %>%
            tibble::rownames_to_column("cell") %>%
            mutate(
                order = cell_order
            ) %>%
            filter(
                sample == target.sample
            ) %>%
            arrange(order) %>%
            mutate(
                order = row_number(),
                min = min(order),
                max = max(order),
                projected = 1 / (max - min) * (order - min)
            ) %>%
            select(cell, sample, projected, color.sample)
        if (cnt == 0) {
            all.orders <- orders
        } else {
            all.orders <- all.orders %>%
                rbind(orders)
        }
        cnt <- cnt + 1
    }

    cnt <- 0
    for (target.feature in target.features) {
        plt <-
            assigned_seurat %>%
            set_default_assay(assay) %>%
            FeaturePlot(
                features = target.feature,
                slot = "data",
                raster = FALSE
            )
        data <-
            plt$data %>%
            tibble::rownames_to_column("cell") %>%
            rename(
                expression = dplyr::starts_with(target.feature)
            ) %>%
            mutate(
                expression = as.double(expression),
                feature = target.feature
            ) %>%
            select(cell, expression, feature) %>%
            left_join(
                all.orders,
                by = "cell"
            )
        if (cnt == 0) {
            all.data <- data
        } else {
            all.data <- all.data %>%
                rbind(data)
        }
        cnt <- cnt + 1
    }
    cnt <- 0
    for (target.sample in all.samples) {
        plt <-
            all.data %>%
            filter(
                sample == target.sample
            ) %>%
            mutate(
                feature = factor(
                    feature,
                    levels = target.features
                )
            ) %>%
            ggplot(
                aes(
                    x = projected,
                    y = expression,
                    color = sample
                )
            ) +
            geom_smooth(
                method = "loess",
                level = 0.2,
                formula = y ~ x
            ) +
            facet_grid(
                feature ~ .,
                scales = "free"
            ) +
            theme_minimal() +
            theme(
                axis.ticks.y.left = element_blank(),
                axis.title.y.left = element_blank(),
                axis.text.x.top = element_blank(),
                axis.ticks.x.top = element_blank(),
                legend.position = "none",
                plot.background = element_rect(
                    fill = "white",
                    colour = "white"
                ),
                panel.background = element_rect(
                    colour = "#424242"
                ),
                plot.title = element_text(
                    hjust = 0.5,
                    size = 15
                ),
                strip.text.y = element_text(
                    face = "italic"
                ),
                panel.grid = element_blank()
            ) +
            xlab(
                target.sample
            ) +
            scale_color_manual(
                values = all.data %>%
                    arrange(sample) %>%
                    filter(
                        sample == target.sample
                    ) %>%
                    distinct(color.sample) %>%
                    pull()
            ) +
            theme(
                plot.margin = unit(c(0, 0, 0, 0), "cm")
            ) +
            scale_x_continuous(position = "top")
        if (cnt != length(all.samples) - 1) {
            plt <-
                plt +
                theme(
                    strip.background.y = element_blank(),
                    strip.text.y = element_blank()
                )
        }
        if (cnt == 0) {
            all.plt <- plt
        } else {
            all.plt <- all.plt + plt
        }
        cnt <- cnt + 1
    }
    all.plt <-
        all.plt +
        plot_layout(nrow = 1) +
        plot_annotation(
            title = title,
            theme = theme(
                plot.title = element_text(
                    hjust = 0.5,
                    size = 20
                )
            )
        )
    if (is.null(out_lambda)) {
        return(all.plt)
    } else {
        out_lambda(all.plt)
    }
}