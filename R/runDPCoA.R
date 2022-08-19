#' Calculation of Double Principal Correspondance analysis
#'
#' Double Principal Correspondance analysis is made available via the
#' \code{ade4} package in typical fashion. Results are stored in the
#' \code{reducedDims} and are available for all the expected functions.
#'
#' @param x For \code{calculateDPCoA}, a numeric matrix of expression values
#'   where rows are features and columns are cells.
#'   Alternatively, a \code{TreeSummarizedExperiment} containing such a matrix.
#'
#'   For \code{runDPCoA} a \linkS4class{TreeSummarizedExperiment} containing the
#'   expression values as well as a \code{rowTree} to calculate \code{y} using
#'   \code{\link[ape:cophenetic.phylo]{cophenetic.phylo}}.
#'
#' @param y a \code{dist} or a symmetric \code{matrix} compatible with
#'   \code{ade4:dpcoa}
#'
#' @param ncomponents Numeric scalar indicating the number of DPCoA dimensions
#'   to obtain.
#'
#' @param ntop Numeric scalar specifying the number of features with the highest
#'   variances to use for dimensionality reduction. Alternatively \code{NULL},
#'   if all features should be used. (default: \code{ntop = NULL})
#'
#' @param subset_row Vector specifying the subset of features to use for
#'   dimensionality reduction. This can be a character vector of row names, an
#'   integer vector of row indices or a logical vector.
#'
#' @param scale Logical scalar, should the expression values be standardized?
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
#' @param tree_name a single \code{character} value for specifying which
#'   rowTree will be used in calculation. 
#'   (By default: \code{tree_name = "phylo"})
#'
#' @param altexp String or integer scalar specifying an alternative experiment
#'   containing the input data.
#'
#' @param name String specifying the name to be used to store the result in the
#'   reducedDims of the output.
#'
#' @param ... Currently not used.
#'
#' @details
#' In addition to the reduced dimension on the features, the reduced dimension
#' for samples are returned as well as \code{sample_red} attribute.
#' \code{eig}, \code{feature_weights} and \code{sample_weights} are
#' returned as attributes as well.
#'
#' @returns
#' For \code{calculateDPCoA} a matrix with samples as rows and CCA dimensions as
#' columns
#'
#' For \code{runDPCoA} a modified \code{x} with the results stored in
#' \code{reducedDim} as the given \code{name}
#'
#'
#' @name runDPCoA
#' @seealso
#' \code{\link[scater:plotReducedDim]{plotReducedDim}}
#' \code{\link[SingleCellExperiment:reducedDims]{reducedDims}}
#'
#' @examples
#' data(esophagus)
#' dpcoa <- calculateDPCoA(esophagus)
#' head(dpcoa)
#'
#' esophagus <- runDPCoA(esophagus)
#' reducedDims(esophagus)
#'
#' library(scater)
#' plotReducedDim(esophagus, "DPCoA")
NULL

#' @export
#' @rdname runDPCoA
setGeneric("calculateDPCoA", signature = c("x", "y"),
           function(x, y, ...)
               standardGeneric("calculateDPCoA"))

.calculate_dpcoa <- function(x, y, ncomponents = 2, ntop = NULL,
                             subset_row = NULL, scale = FALSE,
                             transposed = FALSE, ...)
{
    .require_package("ade4")
    # input check
    # Check ncomponents
    if( !.is_an_integer(ncomponents)  ){
        stop("'ncomponents' must be a single integer value specifying the number ",
             "of DPCoA dimensions.", call. = FALSE)
    }
    # Check ntop
    if( !(is.null(ntop) || .is_an_integer(ntop))  ){
        stop("'ntop' must be NULL or a single integer value specifying the number ",
             "of features with the highest variance.", call. = FALSE)
    }
    # Check subset_row
    y <- as.matrix(y)
    if(length(unique(dim(y))) != 1L){
        stop("'y' must be symmetric.", call. = FALSE)
    }
    #
    # Get NAs. ade4:dpcoa lead to an error if there are any NAs
    if( any( is.na(x) ) ){
        stop("'x' includes NAs. Please try to convert them into numeric values.",
             call. = FALSE)
    }
    if(!transposed) {
        if(is.null(ntop)){
            ntop <- nrow(x)
        }
        x <- .get_mat_for_reddim(x, subset_row = subset_row, ntop = ntop,
                                 scale = scale)
    }
    y <- y[rownames(y) %in% colnames(x),
           colnames(y) %in% colnames(x),
           drop = FALSE]
    if(nrow(y) != ncol(x)){
        stop("x and y must have corresponding dimensions.", call. = FALSE)
    }
    y <- sqrt(y)
    y <- as.dist(y)
    #
    dpcoa <- ade4::dpcoa(data.frame(x), y, scannf = FALSE, nf = ncomponents)
    ans <- as.matrix(dpcoa$li)
    rownames(ans) <- rownames(x)
    colnames(ans) <- NULL
    attr(ans,"eig") <- dpcoa$eig
    tmp <- as.matrix(dpcoa$dls)
    rownames(tmp) <- colnames(x)
    colnames(tmp) <- NULL
    attr(ans,"sample_red") <- tmp
    attr(ans,"feature_weights") <- unname(dpcoa$dw)
    attr(ans,"sample_weights") <- unname(dpcoa$lw)
    ans
}

#' @export
#' @rdname runDPCoA
setMethod("calculateDPCoA", c("ANY","ANY"), .calculate_dpcoa)

#' @export
#' @importFrom ape cophenetic.phylo
#' @rdname runDPCoA
setMethod("calculateDPCoA", signature = c("TreeSummarizedExperiment","missing"),
    function(x, ..., assay_name = abund_values, abund_values = exprs_values, 
             exprs_values = "counts", tree_name = "phylo")
    {
        .require_package("ade4")
        # Check assay_name
        .check_assay_present(assay_name, x)
        # Check tree_name
        if( !(.is_a_string(tree_name) && tree_name %in% names(x@rowTree)) ){
            stop("'tree_name' must specify a rowTree from 'x'.",
                 call. = FALSE)
        }
        #
        # Get tree
        tree <- rowTree(x, tree_name)
        # Select only those features that are in the rowTree
        whichTree <- rowLinks(x)[ , "whichTree"] == tree_name
        if( any(!whichTree) ){
            warning("Not all rows were present in the rowTree specified by 'tree_name'.",
                    "'x' is subsetted.", call. = FALSE)
            # Subset the data
            x <- x[ whichTree, ]
        }
        dist <- cophenetic.phylo(tree)
        # Get assay
        mat <- assay(x, assay_name)
        calculateDPCoA(mat, dist, ...)
    }
)

#' @export
#' @rdname runDPCoA
#' @importFrom SingleCellExperiment reducedDim<-
runDPCoA <- function(x, ..., altexp = NULL, name = "DPCoA"){
    if (!is.null(altexp)) {
        y <- altExp(x, altexp)
    } else {
        y <- x
    }
    reducedDim(x, name) <- calculateDPCoA(y, ...)
    x
}
