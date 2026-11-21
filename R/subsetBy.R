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
#' @param collapse \code{Logical scalar}. Should the features that do not
#' meet the \code{rows} condition be merged into one feature instead of being
#' removed? (Default: \code{FALSE})
#'
#' @param collapse.name \code{Character scalar}. Name of the merged feature
#' when \code{collapse = TRUE}. (Default: \code{"Other"})
#'
#' @param ... additional arguments passed to
#' \code{\link[=agglomerateByVariable]{agglomerateByVariable}} when
#' \code{collapse = TRUE}, e.g. \code{update.tree}.
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
#' With \code{collapse = TRUE}, the features that do not meet the \code{rows}
#' condition are merged into a single feature named \code{collapse.name},
#' placed after the kept features. Its values are the sums of the merged
#' features, so the sample totals stay the same. This is useful, for example,
#' for plotting the most abundant features together with the rest. Summing is
#' meaningful for counts and relative abundances but not for transformed
#' values such as clr. The merging is done with
#' \code{\link[=agglomerateByVariable]{agglomerateByVariable}}, which drops the
#' \code{rowTree} unless \code{update.tree = TRUE}. \code{collapse} affects
#' only features; samples that do not meet the \code{cols} condition are
#' always removed.
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
#' # Keep the features with a mean relative abundance above 1% and merge the
#' # rest into one feature called "Other"
#' tse <- transformAssay(tse, method = "relabundance")
#' tse_sub <- subsetBy(
#'     tse, rows = rowMeans(relabundance) > 0.01, collapse = TRUE)
#' tail(rownames(tse_sub))
#'
#' # Combine different subsettings in a pipe: first remove samples with
#' # fewer than 500,000 reads, then keep the 5 most abundant features and
#' # merge the rest into "Other". The result can be plotted, e.g., with
#' # miaViz::plotBoxplot() or miaViz::plotAbundance().
#' tse_top <- GlobalPatterns |>
#'     subsetBy(cols = colSums(counts) >= 500000) |>
#'     transformAssay(method = "relabundance") |>
#'     subsetBy(rows = rank(-rowMeans(relabundance)) <= 5, collapse = TRUE)
#' dim(tse_top)
#' rownames(tse_top)
#'
#' # Plot the top features and "Other" by sample type
#' if( requireNamespace("miaViz", quietly = TRUE) ){
#'     miaViz::plotBoxplot(
#'         tse_top, assay.type = "relabundance",
#'         x = "SampleType", fill.by = "SampleType",
#'         facet.by = "rownames", scales = "free_y")
#' }
#'
#' # Conditions can also use rowData and colData variables
#' tse_sub <- subsetBy(
#'     tse, rows = Kingdom == "Bacteria", cols = SampleType == "Feces")
#'
#' # This includes quality control metrics stored with addQC(). Note that
#' # stored metrics are not updated when the data is subsetted, so in stepwise
#' # filtering, compute the metric in the condition or call addQC() again.
#' tse <- addQC(tse, by = "rows")
#' tse <- addQC(tse, by = "cols")
#' tse_sub <- subsetBy(tse, rows = sd > 0 & detected >= 3, cols = total > 1e6)
#'
#' # Keep the top 10 features, i.e., the features with the highest mean
#' # abundance
#' tse_top <- subsetBy(tse, rows = rank(-mean) <= 10)
#'
#' # Keep features that are abundant (over 1%) in at least one sample, or
#' # that have a median count above zero
#' tse_sub <- subsetBy(tse, rows = rowMaxs(relabundance) > 0.01)
#' tse_sub <- subsetBy(tse, rows = rowMedians(counts) > 0)
#'
#' # Remove constant features, i.e., features with zero variance
#' tse_sub <- subsetBy(tse, rows = rowVars(counts) > 0)
#'
#' # Standard deviation is on the scale of the data, e.g., clr
#' tse <- transformAssay(tse, method = "clr", pseudocount = 1)
#' tse_sub <- subsetBy(tse, rows = rowSds(clr) > 1)
#'
#' # Keep features and samples with at most 10% missing values
#' tse_sub <- subsetBy(
#'     tse,
#'     rows = rowMeans(is.na(counts)) <= 0.1,
#'     cols = colMeans(is.na(counts)) <= 0.1)
#'
#' # Keep rare features, i.e., features present in at most 50% of samples
#' tse_sub <- subsetBy(tse, rows = rowMeans(counts > 0) <= 0.5)
#'
#' # Use a detection limit: keep features whose relative abundance is above
#' # 0.1% in more than 20% of samples
#' tse_sub <- subsetBy(tse, rows = rowMeans(relabundance > 0.001) > 0.2)
#'
#' # Calculate prevalence in each sample type. Keep features that are
#' # prevalent in at least one sample type, so that group-specific features
#' # are not removed, e.g., before differential abundance analysis.
#' prev <- sapply(splitOn(tse, group = "SampleType"), getPrevalence)
#' tse_sub <- subsetBy(tse, rows = rowSums(prev > 0.5) > 0)
#'
#' # Keep features that are rare in all sample types
#' tse_sub <- subsetBy(tse, rows = rowSums(prev > 0.5) == 0)
#'
NULL

#' @rdname subsetBy
#' @export
setMethod("subsetBy", signature = c(x = "SummarizedExperiment"),
    function(
        x, rows = NULL, cols = NULL, collapse = FALSE,
        collapse.name = "Other", ...){
        if( !.is_a_bool(collapse) ){
            stop("'collapse' must be TRUE or FALSE.", call. = FALSE)
        }
        if( !.is_a_string(collapse.name) ){
            stop("'collapse.name' must be a single character value.",
                call. = FALSE)
        }
        rows <- rlang::enquo(rows)
        cols <- rlang::enquo(cols)
        if( !rlang::quo_is_null(rows) ){
            i <- .eval_subset_condition(
                rows, x, rowData(x), n = nrow(x), arg = "rows")
            if( collapse ){
                x <- .collapse_rows(x, i, collapse.name, ...)
            } else{
                x <- x[i, ]
            }
        }
        if( !rlang::quo_is_null(cols) ){
            j <- .eval_subset_condition(
                cols, x, colData(x), n = ncol(x), arg = "cols")
            x <- x[, j]
        }
        return(x)
    }
)

# Merge the features that are not in 'keep' into one feature, placed after
# the kept features.
.collapse_rows <- function(x, keep, name, ...){
    # Nothing to merge
    if( length(keep) == nrow(x) ){
        return(x)
    }
    if( is.null(rownames(x)) || anyDuplicated(rownames(x)) ){
        stop("'collapse = TRUE' requires unique rownames.", call. = FALSE)
    }
    if( name %in% rownames(x)[keep] ){
        stop("A kept feature is already named '", name, "'. Choose ",
            "another 'collapse.name'.", call. = FALSE)
    }
    group <- rep(name, nrow(x))
    group[keep] <- rownames(x)[keep]
    # If only one feature is collapsed, every group has a single member and
    # agglomerateByVariable() returns the data unchanged, without renaming the
    # collapsed feature. Rename it directly in that case.
    if( sum(group == name) == 1L ){
        rownames(x) <- group
    } else{
        x <- agglomerateByVariable(x, by = "rows", group = group, ...)
    }
    x <- x[c(group[keep], name), ]
    return(x)
}

# Evaluate a condition with the assays and rowData/colData variables as a
# data mask, and return the indices to keep.
.eval_subset_condition <- function(quo, x, df, n, arg){
    # Assays take precedence over rowData/colData variables with the same name
    data <- as.list(df)
    # rowData/colData can have duplicated or empty column names, which cannot
    # be used in a data mask. The first column with each name is used.
    data <- data[ !duplicated(names(data)) & names(data) != "" ]
    assays <- as.list(assays(x))
    data[names(assays)] <- assays
    res <- rlang::eval_tidy(quo, data = data)
    if( !is.logical(res) || length(res) != n ){
        stop("'", arg, "' must return a logical vector with one value per ",
            if( arg == "rows" ) "feature" else "sample", " (", n, ").",
            call. = FALSE)
    }
    return(which(res))
}
