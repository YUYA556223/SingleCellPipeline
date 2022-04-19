write.tsv <- function(df, path, quote = FALSE, row.names = FALSE) { # nolint
    library(dplyr)
    df %>%
        write.table(
            path,
            sep = "\t",
            quote = quote,
            row.names = row.names
        )
    return()
}