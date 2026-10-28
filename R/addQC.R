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
#' @param ... additional arguments passed to \code{getQC}:
#'
#' \itemize{
#'   \item \code{name}: Specifies the names of the metrics in
#'   the result. Must have the same length as the number of metrics.
#'   (Default: \code{c("total", "detected", "mean", "median", "sd")})
#'
#'   \item \code{na.rm}: \code{Logical scalar}. Should NA values be omitted?
#'   (Default: \code{TRUE})
#' }
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
#' \code{getQC} returns a \code{DataFrame} with one row per sample (or
#' feature). \code{addQC} returns \code{x} with the metrics added to
#' \code{colData(x)} (or \code{rowData(x)}).
#'
#' @name addQC
#'
#' @examples
#' data(GlobalPatterns)
#' tse <- GlobalPatterns
#'
#' # Calculate metrics for samples
#' getQC(tse)
#'
#' # Add metrics for features to rowData
#' tse <- addQC(tse, by = "rows")
#' rowData(tse)
#'
NULL

#' @rdname addQC
#' @export
setMethod("addQC", signature = c(x = "SummarizedExperiment"),
    function(x, by = "cols", ...){
        by <- .check_MARGIN(by)
        res <- getQC(x, by = by, ...)
        x <- .add_values_to_colData(
            x, unname(as.list(res)), colnames(res),
            altexp = list(...)[["altexp"]], MARGIN = by)
        return(x)
    }
)

#' @rdname addQC
#' @export
setMethod("getQC", signature = c(x = "SummarizedExperiment"),
    function(
        x, assay.type = "counts", by = "cols", ...){
        ############################## Input check #############################
        # Support altExp hiddenly
        x <- .check_and_get_altExp(x, ...)
        .check_assay_present(assay.type, x)
        by <- .check_MARGIN(by)
        ############################ Input check end ###########################
        mat <- assay(x, assay.type)
        # Calculate metrics by columns
        if( by == 1L ){
            mat <- t(mat)
        }
        res <- .calculate_qc_metrics(mat, ...)
        return(res)
    }
)

#' @importFrom DelayedMatrixStats colSums2 colMeans2 colSds colQuantiles
.calculate_qc_metrics <- function(
        mat,
        name = c("total", "detected", "mean", "sd", paste0("q", probs*100)),
        probs = seq(from = 0, to = 1, by = 0.25),
        na.rm = FALSE,
        ...
){
    if( !(.is_non_empty_character(name) && length(name) == 4+length(probs) &&
          !anyDuplicated(name)) ){
        stop("'name' must be a unique character vector of length ",
             4+length(probs), ".",
             call. = FALSE)
    }
    if( !is.numeric(probs) ){
        stop("'probs' must be numeric vector.", call. = FALSE)
    }
    if( !.is_a_bool(na.rm) ){
        stop("'na.rm' must be TRUE or FALSE.", call. = FALSE)
    }
    #
    # Calculate metrics
    res <- list(
        total = colSums2(mat, na.rm = na.rm),
        detected = colSums2(mat > 0, na.rm = na.rm),
        mean = colMeans2(mat, na.rm = na.rm),
        sd = colSds(mat, na.rm = na.rm)
    )
    res <- DataFrame(lapply(res, unname), row.names = colnames(mat))
    # Calculate quantiles
    quantiles <- colQuantiles(mat, probs = probs, na.rm = na.rm, ...)
    res <- cbind(res, quantiles)

    # Adjust names
    colnames(res) <- name

    return(res)
}
