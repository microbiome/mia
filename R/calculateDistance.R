#' Calculate sample distances with \code{vegan}
#'
#' \code{calculateDistance} calculates a distance matrix between samples. The
#' type of distance calculate can be modifier by setting \code{FUN}, which
#' expects a function with a matrix input as its first argument.
#'
#' @param x a
#'   \code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#'   object containing a tree.
#'
#' @param FUN a \code{function} for distance calculation. The function must
#'   expect the input matrix as its first argument. Samples have to be rows and
#'   feature columns.
#'
#' @param exprs_values a single \code{character} value for specifying which
#'   assay to use for calculation.
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
#' se <- SummarizedExperiment(assays = list(counts = mat),
#'                            rowData = df)
#' #
#' calculateDistance(se)
NULL


#' @rdname calculateDistance
#' @export
setGeneric("calculateDistance", signature = c("x"),
           function(x, FUN = stats::dist, ...)
             standardGeneric("calculateDistance"))

#' @rdname calculateDistance
#' @export
setMethod("calculateDistance", signature = c(x = "ANY"),
    function(x, FUN = stats::dist, ...){
        do.call(FUN, c(list(x),list(...)))
    }
)

#' @rdname calculateDistance
#' @export
setMethod("calculateDistance", signature = c(x = "SummarizedExperiment"),
    function(x, FUN = stats::dist, exprs_values = "counts", transposed = FALSE,
             ...){
        mat <- assay(x, exprs_values)
        if(!transposed){
            mat <- t(mat)
        }
        calculateDistance(mat, FUN, ...)
    }
)
