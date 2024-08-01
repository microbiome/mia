#' Calculate dissimilarities
#'
#' These functions are designed to calculate dissimilarities on data stored 
#' within a 
#' \code{\link[TreeSummarizedExperiment:TreeSummarizedExperiment-class]{TreeSummarizedExperiment}}
#'  object. For overlap, Unifrac, and Jensen-Shannon Divergence (JSD) 
#'  dissimilarities, the functions use mia internal algorithms, while for other 
#'  types of dissimilarities, they rely on \code{\link[vegan:vegdist]{vegdist}}.
#'
#' @param x a \code{\link[TreeSummarizedExperiment:TreeSummarizedExperiment-class]{TreeSummarizedExperiment}}
#'   object.
#'
#' @param method \code{Character scalar}. Specifies which distance to calculate.
#' 
#' @param name \code{Character scalar}. The name to be used to store the result 
#'  in the reducedDims of the output. (Default: \code{method})
#'
#' @param assay.type \code{Character scalar}. Specifies which assay to use for 
#'   calculation. (Default: \code{"counts"})
#'
#' @param assay_name Deprecated. Use \code{assay.type} instead.
#'   
#' @param exprs_values Deprecated. Use \code{assay.type} instead.
#'   
#' @param transposed \code{Logical scalar}. Specifies if x is transposed with cells in
#'   rows. (Default: \code{FALSE})
#'   
#' @param chunkSize \code{Integer scalar}. Defines the size of data send
#'   to the individual worker. Only has an effect, if \code{BPPARAM} defines
#'   more than one worker. (Default: \code{nrow(x)})
#' 
#' @param BPPARAM A
#'   \code{\link[BiocParallel:BiocParallelParam-class]{BiocParallelParam}}
#'   object specifying whether the calculation should be parallelized.
#'   
#' @param detection \code{Numeric scalar}. Specific to overlap dissimilarity.
#'   Defines detection threshold for absence/presence of features. Feature that 
#'   has abundance under threshold in either of samples, will be discarded when 
#'   evaluating overlap between samples. (Default: \code{0}) 
#'   
#' @param tree if \code{x} is a matrix, a
#'   \code{\link[TreeSummarizedExperiment:phylo]{phylo}} object matching the
#'   matrix. This means that the phylo object and the columns should relate
#'   to the same type of features (aka. microorganisms).
#'
#' @param tree.name \code{Character scalar}. Specific to unifrac dissimilarity. 
#' Specifies the name of the tree used in calculation. (Default: \code{"phylo"})
#' 
#' @param tree_name Deprecated. Use \code{tree.name} instead.
#'   
#' @param weighted \code{Logical scalar}. Should use weighted-Unifrac
#'   calculation? Weighted-Unifrac takes into account the relative abundance of
#'   species/taxa shared between samples, whereas unweighted-Unifrac only
#'   considers presence/absence. Default is \code{FALSE}, meaning the
#'   unweighted-Unifrac distance is calculated for all pairs of samples.
#'   (Default: \code{FALSE})
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
#' data(GlobalPatterns)
#' tse <- GlobalPatterns
#' tse <- addDissimilarity(tse, method = "overlap", detection = 0.25)
#' reducedDim(tse, "overlap")
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
        x, method, exprs_values = "counts", assay_name = exprs_values, 
        assay.type = assay_name, tree.name = "phylo", transposed = FALSE, ...){
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
        x, method, exprs_values = "counts", assay_name = exprs_values, 
        assay.type = assay_name, transposed = FALSE, ...){
    # Input checks
    if(!.is_non_empty_string(assay.type)){
        stop("'assay.type' must be a non-empty single character value",
            call. = FALSE)
    }
    if(!.is_non_empty_string(method)){
        stop("'method' must be a non-empty single character value",
            call. = FALSE)
    }
    if(!.is_a_bool(transposed)){
        stop("'na.rm' must be TRUE or FALSE.", call. = FALSE)
    }
    mat <- assay(x, assay.type)
    if(!transposed){
        mat <- t(mat)
    }
    if( method == "unifrac"){
      node.label <- .get_rowlinks(x, assay.type)
      mat <- .calculate_dissimilarity(mat = x, method = method, 
                                      node.label = node.label, ...)
    }
    else{
      mat <- .calculate_dissimilarity(mat = x, method = method, ...)
    }
    return(mat)
    }
)

#' @rdname getDissimilarity
#' @export
setMethod(
    "getDissimilarity", signature = c(x = "ANY"),
    function(
        x, method, exprs_values = "counts", assay_name = exprs_values, 
        assay.type = assay_name, ...){
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
    x, method, exprs_values = "counts", assay_name = exprs_values, 
    assay.type = assay_name, name = method, transposed = FALSE, ...){
    res <- getDissimilarity(x, method = method, assay.type = assay.type, 
                            transposed = transposed, ...)
    # Input checks
    if(!.is_non_empty_string(assay.type)){
      stop("'assay.type' must be a non-empty single character value",
           call. = FALSE)
    }
    if(!.is_non_empty_string(method)){
      stop("'method' must be a non-empty single character value",
           call. = FALSE)
    }
    if(!.is_non_empty_string(name)){
      stop("'name' must be a non-empty single character value",
           call. = FALSE)
    }
    if(!.is_a_bool(transposed)){
      stop("'na.rm' must be TRUE or FALSE.", call. = FALSE)
    }
    if( !identical(rownames(as.matrix(res)), colnames(assay(x, assay.type))) ){
      warning("Samples of the dissimilarity matrix should be the same as the
            samples in columns of the assay specified with 'assay.type'. The 
            result is not added to reducedDim.")
      return(res)
    }
    else{
      .add_values_to_reducedDims(x, as.matrix(res), name)
    }
      
  }
)

.calculate_dissimilarity <- function(
        mat, method, node.label = NULL, diss.fun = NULL, tree = NULL, ...){
    # input check
    if( !(is.null(diss.fun) || is.function(diss.fun)) ){
        stop("'diss.fun' must be NULL or a function.", call. = FALSE)
    }
    #
    args <- c(list(mat, method = method, node.label = node.label, list(...)))
    # If the dissimilarity functon is not specified, get default choice
    if( is.null(diss.fun) ){
        if( method %in% c("overlap") ){
            diss.fun <- .get_overlap
            message("'diss.fun' defaults to .get_overlap.")
        } else if( method %in% c("unifrac")  ){
            args[["tree"]] <- tree
            diss.fun <- .get_unifrac
            message("'diss.fun' defaults to .get_unifrac.")
        } else if( method %in% c("jsd")  ){
            diss.fun <- .get_jsd
            message("'diss.fun' defaults to .get_jsd.")
        } else if( requireNamespace("vegan") ){
            diss.fun <- vegan::vegdist
            message("'diss.fun' defaults to vegan::vegdist.")
        } else{
            diss.fun <- stats::dist
            message("'diss.fun' defaults to stats::dist.")
        }
    }
    # Calculate dissimilarity with specified function
    res <- do.call(diss.fun, args)
    return(res)
}
