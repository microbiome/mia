#' Subset features or samples by conditions on the data
#'
#' \code{subsetBy} keeps the features (rows) and samples (columns) that meet
#' conditions written as R expressions. The conditions are evaluated on the
#' object that is passed to the function, so when several \code{subsetBy}
#' calls are piped, each condition is computed from the data that remains
#' after the previous steps.
#'
#' @param x a \code{\link[SummarizedExperiment]{SummarizedExperiment}} object.
#'
#' @param rows An expression that returns a \code{logical} vector with one
#' value per feature, or \code{NULL} to keep all features. (Default:
#' \code{NULL})
#'
#' @param cols An expression that returns a \code{logical} vector with one
#' value per sample, or \code{NULL} to keep all samples. (Default:
#' \code{NULL})
#'
#' @param ... additional arguments. Currently not used.
#'
#' @details
#' The expressions can refer to the assays by their names (e.g.
#' \code{counts}), which are available as matrices. \code{rows} can also
#' refer to the columns of \code{rowData(x)}, and \code{cols} to the columns of
#' \code{colData(x)}. If an assay and a \code{rowData} or \code{colData}
#' column have the same name, the assay is used. Other variables, such as
#' thresholds, are taken from the environment where \code{subsetBy} is called.
#'
#' \code{rows} is evaluated before \code{cols}, and \code{cols} is evaluated
#' on the data that remains after subsetting the rows. \code{NA} values in
#' the result of a condition are treated as \code{FALSE}.
#'
#' @return
#' \code{x} with the features and samples that meet the conditions.
#'
#' @seealso
#' \code{\link[=getPrevalence]{getPrevalent}},
#' \code{\link[=addQC]{addQC}}
#'
#' @name subsetBy
#'
#' @examples
#' data(GlobalPatterns)
#' tse <- GlobalPatterns
#'
#' # Keep features with more than 10 counts in total
#' tse_sub <- subsetBy(tse, rows = rowSums(counts) > 10)
#'
#' # Filter in steps. Each step is computed from the data left by the
#' # previous step, so library sizes are calculated after the rare features
#' # have been removed.
#' tse_sub <- tse |>
#'     subsetBy(rows = rowSums(counts) > 10) |>
#'     subsetBy(rows = rowMeans(counts > 0) > 0.25) |>
#'     subsetBy(cols = colSums(counts) > 100000)
#'
#' # Conditions can also use rowData and colData variables
#' tse_sub <- subsetBy(
#'     tse, rows = Kingdom == "Bacteria", cols = SampleType == "Feces")
#'
NULL

#' @rdname subsetBy
#' @export
setMethod("subsetBy", signature = c(x = "SummarizedExperiment"),
    function(x, rows = NULL, cols = NULL, ...){
        rows <- rlang::enquo(rows)
        cols <- rlang::enquo(cols)
        if( !rlang::quo_is_null(rows) ){
            i <- .eval_subset_condition(
                rows, x, rowData(x), n = nrow(x), arg = "rows")
            x <- x[i, ]
        }
        if( !rlang::quo_is_null(cols) ){
            j <- .eval_subset_condition(
                cols, x, colData(x), n = ncol(x), arg = "cols")
            x <- x[, j]
        }
        return(x)
    }
)

# Evaluate a condition with the assays and rowData/colData variables as a
# data mask, and return the indices to keep.
.eval_subset_condition <- function(quo, x, df, n, arg){
    data <- c(as.list(df), as.list(assays(x)))
    res <- rlang::eval_tidy(quo, data = data)
    if( !is.logical(res) || length(res) != n ){
        stop("'", arg, "' must return a logical vector with one value per ",
            if( arg == "rows" ) "feature" else "sample", " (", n, ").",
            call. = FALSE)
    }
    return(which(res))
}
