#' Perform non-metric MDS on sample-level data
#'
#' Perform non-metric multi-dimensional scaling (nMDS) on samples, based on the
#' data in a \code{SingleCellExperiment} object.
#'
#' @param x For \code{calculateNMDS}, a numeric matrix of expression values
#'   where rows are features and columns are cells.
#'   Alternatively, a \code{TreeSummarizedExperiment} containing such a matrix.
#'
#'   For \code{runNMDS} a \linkS4class{SingleCellExperiment}
#'
#' @param ncomponents Numeric scalar indicating the number of NMDS dimensions
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
#' @param keep_dist Logical scalar indicating whether the \code{dist} object
#'   calculated by \code{FUN} should be stored as \sQuote{dist} attribute of
#'   the matrix returned/stored by \code{calculateNMDS}/ \code{runNMDS}.
#'
#' @param transposed Logical scalar, is x transposed with cells in rows?
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
#' @param FUN a \code{function} or \code{character} value with a function
#'   name returning a \code{\link[stats:dist]{dist}} object
#'
#' @param nmdsFUN a \code{character} value to choose the scaling
#'   implementation, either \dQuote{isoMDS} for
#'   \code{\link[MASS:isoMDS]{MASS::isoMDS}} or \dQuote{monoMDS} for
#'   \code{\link[vegan:monoMDS]{vegan::monoMDS}}
#'
#' @param ... additional arguments to pass to \code{FUN} and
#'   \code{nmdsFUN}.
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
#' @return For \code{calculateNMDS}, a matrix is returned containing the MDS
#'   coordinates for each sample (row) and dimension (column).
#'
#' @details
#' Either \code{\link[MASS:isoMDS]{MASS::isoMDS}} or
#' \code{\link[vegan:monoMDS]{vegan::monoMDS}} are used internally to compute
#' the NMDS components. If you supply a custom \code{FUN}, make sure that
#' the arguments of \code{FUN} and \code{nmdsFUN} do not collide.
#'
#' @name runNMDS
#'
#' @seealso
#' \code{\link[MASS:isoMDS]{MASS::isoMDS}},
#' \code{\link[vegan:monoMDS]{vegan::monoMDS}}
#' for NMDS component calculation.
#'
#' \code{\link[scater:plotReducedDim]{plotMDS}}, to quickly visualize the
#' results.
#'
#' @author
#' Felix Ernst
#'
#' @examples
#' # generate some example data
#' mat <- matrix(1:60, nrow = 6)
#' df <- DataFrame(n = c(1:6))
#' tse <- TreeSummarizedExperiment(assays = list(counts = mat),
#'                                 rowData = df)
#' #
#' calculateNMDS(tse)
#'
#' #
#' data(esophagus)
#' esophagus <- runNMDS(esophagus, FUN = vegan::vegdist, name = "BC")
#' esophagus <- runNMDS(esophagus, FUN = vegan::vegdist, name = "euclidean",
#'                      method = "euclidean")
#' reducedDims(esophagus)
NULL


#' @rdname runNMDS
#' @export
setGeneric("calculateNMDS", function(x, ...) standardGeneric("calculateNMDS"))


.format_nmds_isoMDS <- function(nmds){
    ans <- nmds$points
    attr(ans, "Stress") <- nmds$stress
    ans
}

#' @importFrom vegan scores
.format_nmds_monoMDS <- function(nmds){
    ans <- scores(nmds)
    colnames(ans) <- NULL
    attr(ans, "Stress") <- nmds$stress
    ans
}

.format_nmds <- function(nmds, nmdsFUN, sample_names){
    ans <- switch(nmdsFUN,
                  "isoMDS" = .format_nmds_isoMDS(nmds),
                  "monoMDS" = .format_nmds_monoMDS(nmds))
    rownames(ans) <- sample_names
    ans
}

.get_nmds_args_isoMDS <- function(args){
    args[c("maxit","trace","tol","p")]
}

.get_nmds_args_monoMDS <- function(args){
    args[c("model","threshold","maxit","weakties","stress","scaling","pc",
           "smin","sfgrmin","sratmax")]
}

.get_nmds_args <- function(nmdsFUN, ...){
    args <- list(...)
    args <- switch(nmdsFUN,
                   "isoMDS" = .get_nmds_args_isoMDS(args),
                   "monoMDS" = .get_nmds_args_monoMDS(args))
    args <- args[!vapply(args,is.null,logical(1))]
    args
}

#' @importFrom MASS isoMDS
#' @importFrom stats cmdscale
#' @importFrom vegan vegdist monoMDS
.calculate_nmds <- function(x, FUN = vegdist,
                            nmdsFUN = c("isoMDS","monoMDS"),
                            ncomponents = 2, ntop = 500, subset_row = NULL,
                            scale = FALSE, transposed = FALSE,
                            keep_dist = FALSE, ...){
    nmdsFUN <- match.arg(nmdsFUN)
    nmdsArgs <- .get_nmds_args(nmdsFUN, ...)
    if(!transposed) {
        x <- .get_mat_for_reddim(x, subset_row = subset_row, ntop = ntop,
                                 scale = scale)
    }
    x <- as.matrix(x)
    sample_names <- rownames(x)
    sample_dist <- do.call(FUN,
                           c(list(x),
                             list(...)))
    attributes(sample_dist) <- attributes(sample_dist)[c("class","Size")]
    y <- cmdscale(sample_dist, k = ncomponents)
    ans <- do.call(nmdsFUN,
                   c(list(sample_dist, y = y, k = ncomponents),
                     nmdsArgs))
    ans <- .format_nmds(ans, nmdsFUN, sample_names)
    if (keep_dist) {
        attr(ans,"dist") <- sample_dist
    }
    ans
}

#' @rdname runNMDS
#' @export
setMethod("calculateNMDS", "ANY", .calculate_nmds)

#' @rdname runNMDS
#' @importFrom SummarizedExperiment assay
#' @export
setMethod("calculateNMDS", "SummarizedExperiment",
    function(x, ..., assay_name = abund_values, abund_values = exprs_values, 
             exprs_values = "counts", FUN = vegdist) {
        .calculate_nmds(assay(x, assay_name), FUN = FUN, ...)
    }
)

#' @rdname runNMDS
#' @export
setMethod("calculateNMDS", "SingleCellExperiment",
    function(x, ..., assay_name = abund_values, abund_values = exprs_values, 
             exprs_values = "counts", dimred = NULL, n_dimred = NULL,
             FUN = vegdist){
        mat <- .get_mat_from_sce(x, exprs_values = assay_name,
                                 dimred = dimred, n_dimred = n_dimred)
        calculateNMDS(mat, transposed = !is.null(dimred), FUN = FUN,...)
    }
)

#' @rdname runNMDS
#' @export
#' @importFrom SingleCellExperiment reducedDim<- altExp
runNMDS <- function(x, ..., altexp = NULL, name = "NMDS") {
    if (!is.null(altexp)) {
        y <- altExp(x, altexp)
    } else {
        y <- x
    }
    reducedDim(x, name) <- calculateNMDS(y, ...)
    x
}

#' @rdname runNMDS
#' @export
plotNMDS <- function(x, ..., ncomponents = 2){
    plotReducedDim(x, ncomponents = ncomponents, dimred = "NMDS",
                   ...)
}
