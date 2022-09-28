#' Calculate sample distances with \code{vegan}
#' 
#' Will be removed by Bioc 3.15
#'
#' \code{calculateDistance} calculates a distance matrix between samples. The
#' type of distance calculated can be modified by setting \code{FUN}, which
#' expects a function with a matrix input as its first argument.
#'
#' @param x a
#'   \code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#'   object containing a tree.
#'
#' @param FUN a \code{function} for distance calculation. The function must
#'   expect the input matrix as its first argument. With rows as samples 
#'   and columns as features.
#'
#' @param assay_name a single \code{character} value for specifying which
#'   assay to use for calculation.
#'
#' @param exprs_values a single \code{character} value for specifying which
#'   assay to use for calculation. 
#'   (Please use \code{assay_name} instead.)
#'   
#' @param abund_values a single \code{character} value for specifying which
#'   assay to use for calculation.
#'   (Please use \code{assay_name} instead. At some point \code{abund_values}
#'   will be disabled.)
#'   
#' @param transposed Logical scalar, is x transposed with cells in rows?
#'
#' @param ... other arguments passed onto \code{FUN}
#'
#' @return a sample-by-sample distance matrix, suitable for NMDS, etc.
#'
#' @name calculateDistance
#'
#' @export
#'
#' @examples
#' # generate some example data
#' mat <- matrix(1:60, nrow = 6)
#' df <- DataFrame(n = c(1:6))
#' tse <- TreeSummarizedExperiment(assays = list(counts = mat),
#'                                 rowData = df)
#' \dontrun{
#' calculateDistance(tse)
#' }
#' 
NULL


#' @rdname calculateDistance
#' @export
setGeneric("calculateDistance", signature = c("x"),
           function(x, FUN = stats::dist, ...)
             standardGeneric("calculateDistance"))

#' @rdname calculateDistance
#' @export
setMethod("calculateDistance", signature = c(x = "SummarizedExperiment"),
    function(x, FUN = stats::dist, assay_name = abund_values, abund_values = exprs_values, 
             exprs_values = "counts", transposed = FALSE,
             ...){
        mat <- assay(x, assay_name)
        if(!transposed){
            mat <- t(mat)
        }
        calculateDistance(mat, FUN, ...)
    }
)


.calculate_distance <- function(mat, FUN = stats::dist, ...){
    # Distance between all samples against all samples
    # Distance between all samples against all samples
    do.call(FUN, c(list(mat),list(...)))
}


