#' Perform MDS on sample-level data
#'
#' Perform multi-dimensional scaling (MDS) on samples, based on the data in a
#' SingleCellExperiment object.
#'
#' @param x a \linkS4class{SingleCellExperiment} containing a numeric matrix of
#'   log-expression values where rows are features and columns are cells.
#'
#' @param exprs_values a single \code{character} value for specifying which
#'   assay to use for calculation.
#'
#' @param FUN a \code{function} or \code{character} value with a function name
#'   returning a \code{\link[stats:dist]{dist}} object
#'
#' @param ... additional arguments to pass to
#'   \code{\link[SEtup:calculateMDS2]{calculateMDS2}} and \code{FUN}.
#'
#' @param dimred String or integer scalar specifying the existing dimensionality
#'   reduction results to use.
#'
#' @param n_dimred Integer scalar or vector specifying the dimensions to use if
#'   dimred is specified.
#'
#' @param altexp String or integer scalar specifying an alternative experiment
#'   containing the input data.
#'
#' @param name String specifying the name to be used to store the result in the
#'   reducedDims of the output.
#'
#' @return For \code{calculateMDS}, a matrix is returned containing the MDS
#'   coordinates for each cell (row) and dimension (column).
#'
#' @details The function \code{\link{cmdscale}} is used internally to compute
#'   the MDS components.
#'
#' @name runMDS2
#' @seealso
#' \code{\link{cmdscale}}, to perform the underlying calculations.
#'
#' \code{\link[scater:plotReducedDim]{plotMDS}}, to quickly visualize the
#' results.
#'
#' @author Aaron Lun, based on code by Davis McCarthy, modified for flexible
#'   distance function input by Felix G.M. Ernst
#'
#' @importFrom SEtup calculateMDS2
#'
#' @examples
#' data(esophagus, package="MicrobiomeExperiment")
#' esophagus <- runMDS2(esophagus, FUN = calculateUniFrac, name = "UniFrac",
#'                      tree = rowTree(esophagus))
#' reducedDim(esophagus)
NULL

#' @rdname runMDS2
#' @export
setMethod("calculateMDS2", "SingleCellExperiment",
    function(x, ..., exprs_values = "counts", dimred = NULL, n_dimred = NULL,
           FUN = calculateDistance){
        mat <- .get_mat_from_sce(x, exprs_values = exprs_values,
                                 dimred = dimred, n_dimred = n_dimred)
        calculateMDS2(mat, transposed = !is.null(dimred), FUN = FUN,...)
    }
)

#' @rdname runMDS2
#' @export
#' @importFrom SingleCellExperiment reducedDim<- altExp
runMDS2 <- function(x, ..., altexp = NULL, name = "MDS2") {
    if (!is.null(altexp)) {
        y <- altExp(x, altexp)
    } else {
        y <- x
    }
    reducedDim(x, name) <- calculateMDS2(y, ...)
    x
}
