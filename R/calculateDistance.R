#' Calculate sample distances with \code{vegan}
#' 
#' Will be removed by Bioc 3.15
#'
#' \code{getDissimilarity} calculates a distance matrix between samples. The
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
#' @name getDissimilarity
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
#' getDissimilarity(tse)
#' }
#' 
NULL

#' @rdname getDissimilarity
#' @export
setGeneric(
    "getDissimilarity", signature = c("x"), function(x, method, ...)
        standardGeneric("getDissimilarity"))

#' @rdname getDissimilarity
#' @export
setMethod(
    "getDissimilarity", signature = c(x = "TreeSummarizedExperiment"),
    function(
        x, method, assay.type = "counts", tree.name = "phylo", transposed = FALSE, ...){
    mat <- assay(x, assay.type)
    if(!transposed){
        mat <- t(mat)
    }
    tree <- rowTree(tse, tree.name)
    res <- getDissimilarity(mat, method = method, tree = tree, ...)
    return(res)
    }
)

#' @rdname getDissimilarity
#' @export
setMethod(
    "getDissimilarity", signature = c(x = "SummarizedExperiment"),
    function(
        x, method, assay.type = "counts", transposed = FALSE, ...){
    mat <- assay(x, assay.type)
    if(!transposed){
        mat <- t(mat)
    }
    res <- getDissimilarity(mat, method = method, ...)
    return(res)
    }
)

#' @rdname getDissimilarity
#' @export
setMethod(
    "getDissimilarity", signature = c(x = "ANY"),
    function(
        x, method, ...){
    # Input check
    if( !.is_a_string(method) ){
        stop("'method' must be a single character value.", call. = FALSE)
    }
    #
    # Calculate dissimilarity
    mat <- .calculate_dissimilarity(mat = x, method = method, ...)
    return(mat)
    }
)

.calculate_dissimilarity <- function(
        mat, method, diss.fun = NULL, tree = NULL, ...){
    # inout check
    if( !(is.null(diss.fun) || is.function(diss.fun)) ){
        stop("'diss.fun' must be NULL or a function.", call. = FALSE)
    }
    #
    args <- c(list(mat, method = method), list(...))
    # If the dissimilarity functon is not specified, get default choice
    if( is.null(diss.fun) ){
        if( method %in% c("overlap") ){
            diss.fun <- calculateOverlap
            message("'diss.fun' defaults to calculateOverlap.")
        } else if( method %in% c("unifrac")  ){
            args[["tree"]] <- tree
            diss.fun <- calculateUnifrac
            message("'diss.fun' defaults to calculateUnifrac.")
        } else if( method %in% c("jsd")  ){
            diss.fun <- runJSD
            message("'diss.fun' defaults to runJSD.")
        } else if( require("vegan") ){
            diss.fun <- vegan::vegdist
            message("'diss.fun' defaults to vegan::vegdist.")
        } else{
            diss.fun <- stats::dist
            message("'diss.fun' defaults to stats::dist.")
        }
    }
    # Calcuate dissimilarity with specified function
    res <- do.call(diss.fun, args)
    return(res)
}
