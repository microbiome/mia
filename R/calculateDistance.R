#' Calculate dissimilarities
#'
#' These functions are designed to calculate dissimilarities on data stored 
#' within a 
#' \code{\link[TreeSummarizedExperiment:TreeSummarizedExperiment-class]{TreeSummarizedExperiment}}
#'  object. For overlap, Unifrac, and Jensen-Shannon Divergence (JSD) 
#'  dissimilarities, the functions use mia internal functions, while for other 
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
#' @param tree.name \code{Character scalar}. Specific to unifrac dissimilarity. 
#' Specifies the name of the tree used in calculation. (Default: \code{"phylo"})
#' 
#' @param tree_name Deprecated. Use \code{tree.name} instead.
#'
#' @param ... other arguments passed onto \code{\link[vegan:vegdist]{vegdist}}, 
#'  or the following arguments passed onto mia internal functions for overlap,
#'  unifrac and JSD dissimilarities:
#' \itemize{
#'   \item \code{weighted}: Specific to unifrac dissimilarity.
#'   \code{Logical scalar}. Should use weighted-Unifrac calculation? 
#'   Weighted-Unifrac takes into account the relative abundance of
#'   species/taxa shared between samples, whereas unweighted-Unifrac only
#'   considers presence/absence. Default is \code{FALSE}, meaning the
#'   unweighted-Unifrac distance is calculated for all pairs of samples.
#'   (Default: \code{FALSE})
#'   
#'   \item \code{chunkSize}: \code{Integer scalar}. Defines the size of data 
#'   send to the individual worker. Only has an effect, if \code{BPPARAM} 
#'   defines more than one worker. (Default: \code{nrow(x)})
#'   
#'   \item \code{BPPARAM}: 
#'   \code{\link[BiocParallel:BiocParallelParam-class]{BiocParallelParam} object}.
#'   Specifies whether the calculation should be parallelized.
#'   
#'   \item \code{detection}: \code{Numeric scalar}. Specific to overlap dissimilarity.
#'   Defines detection threshold for absence/presence of features. Feature that 
#'   has abundance under threshold in either of samples, will be discarded when 
#'   evaluating overlap between samples. (Default: \code{0}) 
#' }
#'
#' @return 
#' \code{getDissimilarity} returns a sample-by-sample distance matrix, suitable 
#'   for NMDS, etc.
#' 
#' \code{addDissimilarity} returns \code{x} that includes distance matrix in its 
#'   reducedDim. 
#'   
#' @details 
#'   Overlap reflects similarity between sample-pairs. When overlap is 
#'   calculated using relative abundances, the higher the value the higher the 
#'   similarity is. 
#'   When using relative abundances, overlap value 1 means that 
#'   all the abundances of features are equal between two samples, and 0 means 
#'   that samples have completely different relative abundances. 
#'   
#'   Mia unifrac internal function utilizes 
#'   \code{\link[rbiom:unifrac]{rbiom:unifrac()}}.
#'   
#' @name getDissimilarity
#' 
#' @author
#' For overlap implementation: 
#' Leo Lahti and Tuomas Borman. Contact: \url{microbiome.github.io}
#' 
#' For JSD implementation:
#' Susan Holmes \email{susan@@stat.stanford.edu}.
#' Adapted for phyloseq by Paul J. McMurdie.
#' Adapted for mia by Felix G.M. Ernst
#' 
#' @seealso
#' \url{http://en.wikipedia.org/wiki/Jensen-Shannon_divergence}
#' 
#' @references
#' For unifrac dissimilarity: \url{http://bmf.colorado.edu/unifrac/}
#' 
#' See also additional descriptions of Unifrac in the following articles:
#'
#' Lozupone, Hamady and Knight, ``Unifrac - An Online Tool for Comparing
#' Microbial Community Diversity in a Phylogenetic Context.'', BMC
#' Bioinformatics 2006, 7:371
#'
#' Lozupone, Hamady, Kelley and Knight, ``Quantitative and qualitative (beta)
#' diversity measures lead to different insights into factors that structure
#' microbial communities.'' Appl Environ Microbiol. 2007
#'
#' Lozupone C, Knight R. ``Unifrac: a new phylogenetic method for comparing
#' microbial communities.'' Appl Environ Microbiol. 2005 71 (12):8228-35.
#' 
#' For JSD dissimilarity: 
#' Jensen-Shannon Divergence and Hilbert space embedding.
#' Bent Fuglede and Flemming Topsoe University of Copenhagen,
#' Department of Mathematics
#' \url{http://www.math.ku.dk/~topsoe/ISIT2004JSD.pdf}
#'
#' @export
#'
#' @examples
#' # load data
#' data(GlobalPatterns)
#' tse <- GlobalPatterns
#' 
#' ### Overlap dissimilarity
#' 
#' tse <- addDissimilarity(tse, method = "overlap", detection = 0.25)
#' reducedDim(tse, "overlap")[1:6, 1:6]
#' 
#' ### JSD dissimilarity
#' 
#' tse <- addDissimilarity(tse, method = "jsd")
#' reducedDim(tse, "jsd")[1:6, 1:6]
#' 
#' # Multi Dimensional Scaling applied to JSD distance matrix
#' tse <- runMDS(tse, FUN = getDissimilarity, method = "overlap", 
#'               assay.type = "counts")
#' reducedDim(tse, "MDS")[1:6, ]
#'               
#' ### Unifrac dissimilarity
#' 
#' res <- getDissimilarity(tse, method = "unifrac", weighted = FALSE)
#' dim(as.matrix((res)))
#' 
#' tse <- addDissimilarity(tse, method = "unifrac", weighted = TRUE)
#' reducedDim(tse, "unifrac")[1:6, 1:6]
#' 
NULL

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
    x, method, assay_name = "counts", assay.type = assay_name, name = method,
    transposed = FALSE, tree_name = "phylo", tree.name = tree_name, ...){
    #
    res <- getDissimilarity(
      x, method = method, assay.type = assay.type, transposed = transposed,
      tree.name = tree.name, ...)
    
    .add_values_to_reducedDims(x, as.matrix(res), name)
  }
)

#' @rdname getDissimilarity
#' @export
setGeneric(
    "getDissimilarity", signature = c("x"), function(x, method, ...)
        standardGeneric("getDissimilarity"))

#' @rdname getDissimilarity
#' @export
setMethod(
    "getDissimilarity", signature = c(x = "SummarizedExperiment"),
    function(
        x, method, exprs_values = "counts", assay_name = exprs_values, 
        assay.type = assay_name, transposed = FALSE, tree_name = "phylo",
        tree.name = tree_name, ...){
    # Input checks
    .check_assay_present(assay.type, x)
    if( !.is_non_empty_string(method) ){
        stop("'method' must be a non-empty single character value",
            call. = FALSE)
    }
    if( !.is_a_bool(transposed) ){
        stop("'na.rm' must be TRUE or FALSE.", call. = FALSE)
    }
    #
    # If mrthod is unifrac, the object is TreeSE and tree was not provided by
    # user, get tree arguments from TreeSE in addition to matrix and method.
    if( method %in% c("unifrac") && !"tree" %in% names(list(...)) &&
            is(x, "TreeSummarizedExperiment") ){
        args <- .get_tree_args(
            x,  method = method, assay.type = assay.type,
            transposed = transposed, tree.name = tree.name, ...)
    } else{
      # For other methods, get only matrix and method for arguments.
        mat <- assay(x, assay.type)
        if( !transposed ){
            mat <- t(mat)
        }
        args <- c(list(x = mat, method = method), list(...))
    }
    # Calculate dissimilarity based on matrix
    mat <- do.call(getDissimilarity, args)
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

.calculate_dissimilarity <- function(
        mat, method, node.label = NULL, diss.fun = NULL, tree = NULL, ...){
    # input check
    if( !(is.null(diss.fun) || is.function(diss.fun)) ){
        stop("'diss.fun' must be NULL or a function.", call. = FALSE)
    }
    #
    args <- c(list(mat, method = method), list(...))
    # If the dissimilarity functon is not specified, get default choice
    if( is.null(diss.fun) ){
        if( method %in% c("overlap") ){
            diss.fun <- .get_overlap
            message("'diss.fun' defaults to .get_overlap.")
        } else if( method %in% c("unifrac")  ){
            args <- c(args, list(tree = tree, node.label = node.label))
            diss.fun <- .get_unifrac
            message("'diss.fun' defaults to .get_unifrac.")
        } else if( method %in% c("jsd")  ){
            diss.fun <- .get_jsd
            message("'diss.fun' defaults to mia:::.get_jsd.")
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

.get_tree_args <- function(
        x, method, assay.type = assay_name, assay_name = exprs_values, 
        exprs_values = "counts", tree.name = tree_name, 
        tree_name = "phylo", transposed = FALSE, ...){
    # Get functions and parameters based on direction
    tree_present_FUN <- if (transposed) .check_colTree_present
    else .check_rowTree_present
    tree_FUN <- if (transposed) colTree else rowTree
    links_FUN <- if (transposed) colLinks else rowLinks
    margin_name <- if (transposed) "col" else "row"
    # Check assay.type
    .check_assay_present(assay.type, x)
    # Check tree.name
    tree_present_FUN(tree.name, x)
    #
    # Select only those features/samples that are in the tree
    links <- links_FUN(x)
    present_in_tree <- links[, "whichTree"] == tree.name
    if( any(!present_in_tree) ){
      warning(
        "Not all ", margin_name, "s were present in the ", margin_name,
        "Tree specified by 'tree.name'. 'x' is subsetted.",
        call. = FALSE)
      # Subset the data
      if( transposed ){
        x <- x[, present_in_tree]
      } else{
        x <- x[present_in_tree, ]
      }
    }
    # Get assay. By default, dissimilarity between samples is calculated. In
    # dissimilarity functions, features must be in columns and samples in rows
    # in this case.
    mat <- assay(x, assay.type)
    if( !transposed ){
        mat <- t(mat)
    }
    # Get tree
    tree <- tree_FUN(x, tree.name)
    # Get links and take only nodeLabs
    links <- links_FUN(x)
    links <- links[ , "nodeLab" ]
    node.label <- links
    
    # Create an arument list that includes matrix, and tree-related parameters.
    args <- list(x = mat, method = method, tree = tree, node.label = node.label)
    args <- c(args, list(...))
    return(args)
}
