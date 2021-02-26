#' Perform MDS on sample-level data
#'
#' Perform multi-dimensional scaling (MDS) on samples, based on the data in a
#' SingleCellExperiment object.
#'
#' @param x a \linkS4class{SingleCellExperiment} containing a numeric matrix of
#'   expression values where rows are features and columns are cells.
#'
#' @param ncomponents Numeric scalar indicating the number of MDS dimensions
#'   to obtain.
#'
#' @param ntop Numeric scalar specifying the number of features with the highest
#'   variances to use for dimensionality reduction.
#'
#' @param subset_row Vector specifying the subset of features to use for
#'   dimensionality reduction. This can be a character vector of row names, an
#'   integer vector of row indices or a logical vector.
#'
#' @param scale Logical scalar, should the expression values be standardized?
#'
#' @param transposed Logical scalar, is x transposed with cells in rows?
#'
#' @param exprs_values a single \code{character} value for specifying which
#'   assay to use for calculation.
#'
#' @param FUN a \code{function} or \code{character} value with a function name
#'   returning a \code{\link[stats:dist]{dist}} object
#'
#' @param ... additional arguments to pass to \code{FUN}.
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
#' @param keep_dist \code{TRUE} or \code{FALSE}: Should the dist object be
#'   returned as attribute of the MDS result? (default: code{keep_dist = FALSE})
#'
#' @return For \code{calculateMDS}, a matrix is returned containing the MDS
#'   coordinates for each cell (row) and dimension (column).
#'
#' @details The function \code{\link{cmdscale}} is used internally to compute
#'   the MDS components.
#'
#' @name runMDS2
#'
#' @seealso
#' \code{\link{cmdscale}}, to perform the underlying calculations.
#'
#' \code{\link[scater:plotReducedDim]{plotMDS}}, to quickly visualize the
#' results.
#'
#' @author Aaron Lun, based on code by Davis McCarthy, modified for flexible
#'   distance function input by Felix G.M. Ernst
#'
#' @examples
#' # generate some example data
#' mat <- matrix(1:60, nrow = 6)
#' df <- DataFrame(n = c(1:6))
#' se <- SummarizedExperiment(assays = list(counts = mat),
#'                            rowData = df)
#' #
#' calculateMDS2(se)
#'
#' #
#' data(esophagus)
#' esophagus <- runMDS2(esophagus, FUN = calculateUniFrac, name = "UniFrac",
#'                      tree = rowTree(esophagus))
#' reducedDim(esophagus)
NULL


#' @rdname runMDS2
#' @export
setGeneric("calculateMDS2", function(x, ...) standardGeneric("calculateMDS2"))

#' @importFrom stats cmdscale dist
.calculate_mds2 <- function(x, FUN = calculateDistance,
                            ncomponents = 2, ntop = 500, subset_row = NULL,
                            scale = FALSE, transposed = FALSE,
                            keep_dist = FALSE, ...){
    if(!transposed) {
        x <- .get_mat_for_reddim(x, subset_row = subset_row, ntop = ntop,
                                 scale = scale)
    }

    x <- as.matrix(x)
    cell_dist <- do.call(FUN, c(list(x),list(...)))
    mds <- cmdscale(cell_dist, k = ncomponents, eig = TRUE)
    ans <- mds$points
    attr(ans,"eig") <- mds$eig
    attr(ans,"GOF") <- mds$GOF
    if (keep_dist) {
        attr(ans,"dist") <- cell_dist
    }
    ans
}

#' @rdname runMDS2
#' @export
setMethod("calculateMDS2", "ANY", .calculate_mds2)

#' @rdname runMDS2
#' @importFrom SummarizedExperiment assay
#' @export
setMethod("calculateMDS2", "SummarizedExperiment",
    function(x, ..., exprs_values = "counts", FUN = calculateDistance) {
        .calculate_mds2(assay(x, exprs_values), FUN = FUN, ...)
    }
)

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
