lambda <- function(..., envir = parent.frame()) {
    if (!require(lazyeval)) stop("Please install.packages('lazyeval')")
    library(lazyeval)
    args <- lazyeval::lazy_dots(...)
    if (length(args) == 0) {
        return(function() {})
    }
    args <- Map(function(x) x$expr, args)
    vars <- unlist(Map(function(x) deparse(x), args[-length(args)]))
    expr <- as.character(args[length(args)])
    expr <- strsplit(expr, ":")[[1]]
    var <- expr[1]
    expr <- paste0(expr[-1], collapse = ":")
    vars <- ifelse(is.null(vars),
        var,
        paste(paste(vars, collapse = ","), var, sep = ",")
    )
    fun_str <- sprintf("function(%s) %s", vars, expr)
    eval(parse(text = fun_str), envir = envir)
}