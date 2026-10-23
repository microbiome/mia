#' Calculate quality control metrics
#'
#' These functions calculate simple quality control metrics, such as library
#' sizes, for samples or features.
#'
#' @param x a \code{\link[SummarizedExperiment]{SummarizedExperiment}} object.
#'
#' @param assay.type \code{Character scalar}. Specifies the name of assay
#' used in calculation. (Default: \code{"counts"})
#'
#' @param by \code{Character scalar}. Determines if metrics are calculated
#' for samples (\code{'cols'}) or features (\code{'rows'}).
#' (Default: \code{"cols"})
#'
#' @param name \code{Character vector}. Specifies the names of the metrics in
#' the result. Must have the same length as the number of metrics.
#' (Default: \code{c("total", "detected", "mean", "median", "sd")})
#'
#' @param ... additional arguments passed to \code{getQCMetrics}.
#'
#' @details
#' The following metrics are calculated for each sample (or feature):
#' \itemize{
#'   \item \code{total}: Sum of values, e.g., library size.
#'   \item \code{detected}: Number of values greater than zero.
#'   \item \code{mean}: Mean of values.
#'   \item \code{median}: Median of values.
#'   \item \code{sd}: Standard deviation of values.
#' }
#'
#' @return
#' \code{getQCMetrics} returns a \code{DataFrame} with one row per sample (or
#' feature). \code{addQCMetrics} returns \code{x} with the metrics added to
#' \code{colData(x)} (or \code{rowData(x)}).
#'
#' @name addQCMetrics
#'
#' @examples
#' data(GlobalPatterns)
#' tse <- GlobalPatterns
#'
#' # Calculate metrics for samples
#' getQCMetrics(tse)
#'
#' # Add metrics for features to rowData
#' tse <- addQCMetrics(tse, by = "rows")
#' rowData(tse)
#'
NULL

#' @rdname addQCMetrics
#' @export
setMethod("addQCMetrics", signature = c(x = "SummarizedExperiment"),
    function(x, by = "cols", ...){
        by <- .check_MARGIN(by)
        res <- getQCMetrics(x, by = by, ...)
        x <- .add_values_to_colData(
            x, unname(as.list(res)), colnames(res),
            altexp = list(...)[["altexp"]], MARGIN = by)
        return(x)
    }
)

#' @rdname addQCMetrics
#' @importFrom DelayedMatrixStats colSums2 colMeans2 colMedians colSds
#' @export
setMethod("getQCMetrics", signature = c(x = "SummarizedExperiment"),
    function(
        x, assay.type = "counts", by = "cols",
        name = c("total", "detected", "mean", "median", "sd"), ...){
        ############################## Input check #############################
        # Support altExp hiddenly
        x <- .check_and_get_altExp(x, ...)
        .check_assay_present(assay.type, x)
        by <- .check_MARGIN(by)
        if( !(.is_non_empty_character(name) && length(name) == 5L &&
                !anyDuplicated(name)) ){
            stop("'name' must be a unique character vector of length 5.",
                call. = FALSE)
        }
        ############################ Input check end ###########################
        mat <- assay(x, assay.type)
        # Calculate metrics by columns
        if( by == 1L ){
            mat <- t(mat)
        }
        res <- list(
            total = colSums2(mat),
            detected = colSums2(mat > 0),
            mean = colMeans2(mat),
            median = colMedians(mat),
            sd = colSds(mat))
        res <- DataFrame(lapply(res, unname), row.names = colnames(mat))
        colnames(res) <- name
        return(res)
    }
)
