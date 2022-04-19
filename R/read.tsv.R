read.tsv <- function(path) { # nolint
    return(
        read.table(path, sep = "\t", header = TRUE)
    )
}