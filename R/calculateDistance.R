#' Calculate dissimilarities
#'
#' These functions calculate dissimilarities on data stored in a 
#'  \code{\link[TreeSummarizedExperiment:TreeSummarizedExperiment-class]{TreeSummarizedExperiment}}
#'  object.
#'
#' @param x a \code{\link[TreeSummarizedExperiment:TreeSummarizedExperiment-class]{TreeSummarizedExperiment}}
#'   object.
#'
#' @param method Character scalar. Specifies which distance to calculate.
#'
#' @param assay.type Character scalar. Specifies which assay to use for 
#'   calculation. (Default: \code{"counts"})
#'
#' @param assay_name Character scalar. Specifies which assay to use for 
#'   calculation. (Please use \code{assay.type} instead. At 
#'   some point \code{assay_name} will be disabled.)
#'
#' @param tree.name Character scalar. Specifies which tree will be used in 
#'   calculation. (Default: \code{"phylo"})
#'   
#' @param transposed Logical scalar. Specifies if x is transposed with cells in
#'   rows. (Default: \code{FALSE})
#'
#' @param ... other arguments passed onto \code{\link[vegan:vegdist]{vegdist}}
#'
#' @return 
#' \code{getDissimilarity} returns a distance matrix.
#' 
#' \code{addDissimilarity} returns a
#'   \code{\link[TreeSummarizedExperiment:TreeSummarizedExperiment-class]{TreeSummarizedExperiment}}
#'   with distance matrix added to reducedDim slot.
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
    tree <- rowTree(x, tree.name)
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

#' @rdname getDissimilarity
#' @export
setGeneric(
  "addDissimilarity", signature = c("x"), function(x, method, ...)
    standardGeneric("addDissimilarity"))

#' @rdname getDissimilarity
#' @export
setMethod(
  "addDissimilarity", signature = c(x = "SummarizedExperiment"),
  function(
    x, method, assay.type = "counts", transposed = FALSE, name = method, ...){
    res <- getDissimilarity(x, method = method, assay.type = assay.type, 
                            transposed = transposed, ...)
    if ( !identical(rownames(as.matrix(res)), colnames(assay(x, assay.type))) ){
      warning("Samples of the dissimilarity matrix should be the same as the
            samples in columns of the assay specified with 'assay.type'. The 
            result is not added to reducedDim.")
      return(res)
    }
    else{
      .add_values_to_reducedDims(x, res, name = name)
    }
      
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
            diss.fun <- getOverlap
            message("'diss.fun' defaults to getOverlap.")
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
