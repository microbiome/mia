#' @name
#' subsetBy
#'
#' @title
#' Subset features and samples
#'
#' @description
#' These functions filter features (rows) and samples (columns) of a
#' \code{SummarizedExperiment} object. They all return the filtered object.
#'
#' \itemize{
#'   \item \code{subsetByAbundance()}: removes features (or samples) with low
#'   abundance, e.g., low total, mean or maximum abundance.
#'   \item \code{subsetByVariance()}: removes features (or samples) with low
#'   variance. By default, removes features with zero variance.
#'   \item \code{subsetBySd()}: removes features (or samples) with low
#'   standard deviation. By default, removes features with zero standard
#'   deviation.
#'   \item \code{subsetByMissing()}: removes features (or samples) with
#'   missing values.
#'   \item \code{subsetByTop()}: keeps the top features based on, e.g., mean
#'   abundance or variance. See \code{\link[=getTop]{getTop()}}.
#' }
#'
#' @details
#' mia provides the following functions to filter data:
#'
#' \tabular{ll}{
#'   \strong{Function} \tab \strong{Keeps} \cr
#'   \code{subsetByAbundance()} \tab Features or samples with total, mean,
#'   median or maximum abundance above a threshold \cr
#'   \code{subsetByVariance()} \tab Features or samples with variance above a
#'   threshold \cr
#'   \code{subsetBySd()} \tab Features or samples with standard deviation
#'   above a threshold \cr
#'   \code{subsetByMissing()} \tab Features or samples with a proportion of
#'   missing values at or below a threshold \cr
#'   \code{subsetByTop()} \tab The \code{top} features based on a statistic \cr
#'   \code{\link[=getPrevalence]{subsetByPrevalent()}} \tab Prevalent
#'   features, optionally in at least one group \cr
#'   \code{\link[=getPrevalence]{subsetByRare()}} \tab Rare features,
#'   optionally in all groups \cr
#' }
#'
#' Functions with the prefix \code{get} return the names of features instead
#' of the filtered object, e.g., \code{\link[=getTop]{getTop()}},
#' \code{\link[=getPrevalence]{getPrevalent()}} and
#' \code{\link[=getPrevalence]{getRare()}}. To merge non-prevalent features
#' instead of removing them, see
#' \code{\link[=agglomerateByPrevalence]{agglomerateByPrevalence()}}.
#'
#' \strong{subsetByAbundance}
#'
#' Features (or samples) whose abundance is strictly greater than
#' \code{threshold} are kept. The abundance is summarized over samples (or
#' features) with \code{method}: \code{"sum"} (total abundance),
#' \code{"mean"}, \code{"median"} or \code{"max"}. For example, use
#' \code{method = "mean"} with relative abundances to keep features whose mean
#' relative abundance exceeds a threshold, or \code{method = "max"} to keep
#' features that are abundant in at least one sample. With the
#' default \code{threshold = 0} and \code{method = "sum"}, features that are
#' absent from all samples are removed. To filter both samples and features,
#' call the function twice with different \code{MARGIN}. To filter features by
#' prevalence, use \code{\link[=getPrevalence]{subsetByPrevalent()}}.
#'
#' \strong{subsetByVariance}
#'
#' Features (or samples) whose variance is strictly greater than
#' \code{threshold} are kept. With the default \code{threshold = 0}, features
#' with zero variance, i.e., constant features, are removed. Features whose
#' variance cannot be calculated (e.g., all values are missing) are removed.
#'
#' \strong{subsetBySd}
#'
#' Same as \code{subsetByVariance()} but uses standard deviation, so the
#' \code{threshold} is on the same scale as the abundances. With the default
#' \code{threshold = 0}, both functions give the same result.
#'
#' \strong{subsetByMissing}
#'
#' Features (or samples) whose proportion of missing values (\code{NA}) is
#' less than or equal to \code{threshold} are kept. With the default
#' \code{threshold = 0}, all features with missing values are removed.
#'
#' \strong{subsetByTop}
#'
#' Features are selected with \code{\link[=getTop]{getTop()}}. The original
#' order of features is preserved.
#'
#' @return
#' Filtered \code{SummarizedExperiment} object.
#'
#' @inheritParams getDissimilarity
#' @inheritParams getPrevalence
#' @inheritParams summary
#'
#' @param threshold \code{Numeric scalar}. Features (or samples) with
#' abundance (\code{subsetByAbundance}), variance (\code{subsetByVariance})
#' or standard deviation (\code{subsetBySd}) less than or equal to this value
#' are removed. For \code{subsetByMissing}, the maximum allowed proportion of
#' missing values (between 0 and 1). (Default: \code{0})
#'
#' @param method \code{Character scalar}. Specifies how abundances are
#' summarized in \code{subsetByAbundance}. Must be one of \code{"sum"},
#' \code{"mean"}, \code{"median"} or \code{"max"}. (Default: \code{"sum"})

#'
#' @param MARGIN \code{Character scalar}. Determines whether features
#' (\code{"features"}) or samples (\code{"samples"}) are filtered.
#' (Default: \code{"features"})
#'
#' @param na.rm \code{Logical scalar}. Should missing values be ignored when
#' calculating the statistics? (Default: \code{FALSE})
#'
#' @param ... additional arguments.
#' \itemize{
#'   \item \code{subsetByTop()}: passed to \code{\link[=getTop]{getTop()}},
#'   e.g., \code{detection} with \code{method = "prevalence"}.
#' }
#'
#' @examples
#' data(GlobalPatterns)
#' tse <- GlobalPatterns
#'
#' # Remove features that are absent from all samples
#' tse_sub <- subsetByAbundance(tse)
#' dim(tse_sub)
#'
#' # Remove samples with total counts of 100,000 or less, and features with
#' # total counts of 10 or less
#' tse_sub <- tse |>
#'     subsetByAbundance(threshold = 1e5, MARGIN = "samples") |>
#'     subsetByAbundance(threshold = 10)
#' dim(tse_sub)
#'
#' # Keep features whose mean relative abundance exceeds 0.1%
#' tse_sub <- tse |>
#'     transformAssay(method = "relabundance") |>
#'     subsetByAbundance(
#'         assay.type = "relabundance", method = "mean", threshold = 0.001)
#' dim(tse_sub)
#'
#' # Keep features that exceed 1% relative abundance in at least one sample
#' tse_sub <- tse |>
#'     transformAssay(method = "relabundance") |>
#'     subsetByAbundance(
#'         assay.type = "relabundance", method = "max", threshold = 0.01)
#' dim(tse_sub)
#'
#' # Additionally, keep features that are present in over 30% of samples
#' tse_sub <- subsetByPrevalent(tse_sub, prevalence = 0.3, detection = 0)
#' dim(tse_sub)
#'
#' # Remove features with zero variance
#' tse_sub <- subsetByVariance(tse)
#' dim(tse_sub)
#'
#' # Keep features whose variance of relative abundances exceeds a threshold
#' tse <- transformAssay(tse, method = "relabundance")
#' tse_sub <- subsetByVariance(
#'     tse, assay.type = "relabundance", threshold = 1e-4)
#' dim(tse_sub)
#'
#' # Keep features whose standard deviation of relative abundances exceeds 1%
#' tse_sub <- subsetBySd(tse, assay.type = "relabundance", threshold = 0.01)
#' dim(tse_sub)
#'
#' # Add missing values to the first 5 features in 2 samples
#' assay(tse, "na") <- assay(tse, "counts")
#' assay(tse, "na")[1:5, 1:2] <- NA
#'
#' # Remove features with missing values
#' tse_sub <- subsetByMissing(tse, assay.type = "na")
#' dim(tse_sub)
#'
#' # Allow missing values in up to 10% of samples
#' tse_sub <- subsetByMissing(tse, assay.type = "na", threshold = 0.1)
#' dim(tse_sub)

#'
#' # Keep the 10 most variable features
#' tse_sub <- subsetByTop(tse, top = 10, method = "sd")
#' rownames(tse_sub)
#'
#' # Keep the 5 most abundant genera
#' tse_sub <- subsetByTop(tse, top = 5, rank = "Genus", na.rm = TRUE)
#' rownames(tse_sub)
#'
#' @seealso
#' \code{\link[=getPrevalence]{subsetByPrevalent}},
#' \code{\link[=getTop]{getTop}},
#' \code{\link[=agglomerateByPrevalence]{agglomerateByPrevalence}}
#'
NULL

############################### subsetByAbundance ##############################

#' @rdname subsetBy
#' @importFrom DelayedMatrixStats rowSums2 colSums2 rowMeans2 colMeans2
#'   rowMedians colMedians rowMaxs colMaxs
#' @export
setMethod("subsetByAbundance", signature = c(x = "SummarizedExperiment"),
    function(
        x,
        assay.type = "counts",
        threshold = 0,
        method = "sum",
        MARGIN = "features",
        na.rm = FALSE,
        ...
    ) {
        FUNS <- list(
            sum = list(rowSums2, colSums2),
            mean = list(rowMeans2, colMeans2),
            median = list(rowMedians, colMedians),
            max = list(rowMaxs, colMaxs)
        )
        if( !(.is_a_string(method) && method %in% names(FUNS)) ){
            stop("'method' must be one of: '",
                paste0(names(FUNS), collapse = "', '"), "'.", call. = FALSE)
        }
        FUNS <- FUNS[[method]]
        x <- .subset_by_statistic(
            x, assay.type, threshold, MARGIN, na.rm, FUNS[[1L]], FUNS[[2L]])
        return(x)
    }
)

############################### subsetByVariance ###############################

#' @rdname subsetBy
#' @importFrom DelayedMatrixStats rowVars colVars
#' @export
setMethod("subsetByVariance", signature = c(x = "SummarizedExperiment"),
    function(
        x,
        assay.type = "counts",
        threshold = 0,
        MARGIN = "features",
        na.rm = FALSE,
        ...
    ) {
        x <- .subset_by_statistic(
            x, assay.type, threshold, MARGIN, na.rm, rowVars, colVars)
        return(x)
    }
)

################################## subsetBySd ##################################

#' @rdname subsetBy
#' @importFrom DelayedMatrixStats rowSds colSds
#' @export
setMethod("subsetBySd", signature = c(x = "SummarizedExperiment"),
    function(
        x,
        assay.type = "counts",
        threshold = 0,
        MARGIN = "features",
        na.rm = FALSE,
        ...
    ) {
        x <- .subset_by_statistic(
            x, assay.type, threshold, MARGIN, na.rm, rowSds, colSds)
        return(x)
    }
)

############################### subsetByMissing ################################

#' @rdname subsetBy
#' @importFrom DelayedMatrixStats rowMeans2 colMeans2
#' @export
setMethod("subsetByMissing", signature = c(x = "SummarizedExperiment"),
    function(
        x,
        assay.type = "counts",
        threshold = 0,
        MARGIN = "features",
        ...
    ) {
        .check_assay_present(assay.type, x)
        if( !(.is_a_numeric(threshold) && !is.na(threshold) &&
                threshold >= 0 && threshold <= 1) ){
            stop("'threshold' must be a single numeric value between 0 ",
                "and 1.", call. = FALSE)
        }
        MARGIN <- .check_MARGIN(MARGIN)
        #
        missing <- is.na(assay(x, assay.type))
        FUN <- if( MARGIN == 1L ) rowMeans2 else colMeans2
        keep <- FUN(missing) <= threshold
        x <- if( MARGIN == 1L ) x[keep, ] else x[, keep]
        return(x)
    }
)

################################# subsetByTop ##################################

#' @rdname subsetBy
#' @export
setMethod("subsetByTop", signature = c(x = "SummarizedExperiment"),
    function(x, top = 5L, rank = NULL, ...){
        x <- .merge_features(x, rank = rank, ...)
        # Rownames must be set and unique to map the names back to rows
        x_named <- .norm_rownames(x)
        features <- getTop(x_named, top = top, ...)
        # Preserve the original order
        index <- sort(match(features, rownames(x_named)))
        x <- x[index, ]
        return(x)
    }
)

################################ HELP FUNCTIONS ################################

# Keeps features (or samples) whose statistic is strictly greater than the
# threshold. Features whose statistic cannot be calculated (NA) are removed.
.subset_by_statistic <- function(
        x, assay.type, threshold, MARGIN, na.rm, row_FUN, col_FUN){
    .check_assay_present(assay.type, x)
    if( !(.is_a_numeric(threshold) && !is.na(threshold) && threshold >= 0) ){
        stop("'threshold' must be a single non-negative numeric value.",
            call. = FALSE)
    }
    MARGIN <- .check_MARGIN(MARGIN)
    if( !.is_a_bool(na.rm) ){
        stop("'na.rm' must be TRUE or FALSE.", call. = FALSE)
    }
    #
    mat <- assay(x, assay.type)
    FUN <- if( MARGIN == 1L ) row_FUN else col_FUN
    values <- FUN(mat, na.rm = na.rm)
    keep <- !is.na(values) & values > threshold
    x <- if( MARGIN == 1L ) x[keep, ] else x[, keep]
    return(x)
}

# Gets sample groups from colData as a factor
.get_sample_groups <- function(x, col.var){
    if( !(.is_a_string(col.var) && col.var %in% colnames(colData(x))) ){
        stop("'col.var' must specify a column from colData(x).", call. = FALSE)
    }
    groups <- colData(x)[[col.var]]
    if( anyNA(groups) ){
        stop("'col.var' must not contain missing values.", call. = FALSE)
    }
    groups <- factor(groups)
    return(groups)
}
