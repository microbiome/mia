#' @name
#' addMDS
#'
#' @title
#' Perform multi-dimensional scaling (MDS)
#'
#' @description
#' Perform multi-dimensional scaling (MDS) also know as Principal Coordinate
#' Analysis (PCoA). These functions are wrappers for
#' \code{\link[scater:runMDS]{scater::calculateMDS}}.
#'
#' @details
#' These functions are wrappers for
#' \code{\link[scater:runMDS]{scater::calculateMDS}} and
#' \code{\link[scater:runMDS]{scater::runMDS}}. While \code{getMDS}
#' returns the results, \code{addMDS} adds them to \code{reducedDim(x)}. The
#' difference is that these functions apply microbiome-specific options such
#' as \code{getDissimilarity} function by default.
#'
#' See \code{\link[scater:runMDS]{scater::calculateMDS}} for details.
#'
#' @return
#' \code{getMDS} returns a MDS results.
#' \code{addMDS} returns a \code{x} with MDS results added to its
#' \code{reducedDim(x, name)}.
#'
#' @inheritParams addAlpha
#'
#' @param assay.type \code{Character scalar}. Specifies the name of assay
#' used in calculation. (Default: \code{"counts"})
#'
#' @param name \code{Character scalar}. A name for the \code{reducedDim()}
#' where results will be stored. (Default: \code{"MDS"})
#'
#' @param ... additional arguments.
#' \itemize{
#'     \item \code{FUN}: \code{Function}. A function that is applied to
#'     calculate dissimilarity. (Default: \code{getDissimilarity})
#'
#'     \item \code{subset.result}: \code{Logical result}. Specifies whether to
#'     subset \code{x} to match the result if some samples were removed during
#'     calculation. (Default: \code{TRUE})
#' }
#'
#' @examples
#' library(mia)
#' library(scater)
#' library(patchwork)
#'
#' data(GlobalPatterns)
#' tse <- GlobalPatterns
#'
#' # Calculate PCoA with Bray-Curtis dissimilarity
#' tse <- transformAssay(tse, method = "relabundance")
#' tse <- addMDS(tse, assay.type = "relabundance", method = "bray")
#'
#' # Calculate PCoA with Unifrac and rarefaction. (Note: increase iterations)
#' tse <- addMDS(tse, method = "unifrac", name = "unifrac")
#'
#' # Calculate PCoA with Unifrac and rarefaction. (Note: increase iterations)
#' tse <- addMDS(tse, method = "unifrac", name = "unifrac_rare", niter = 2L)
#'
#' # Visualize results
#' p1 <- plotReducedDim(tse, "unifrac", colour_by = "SampleType") +
#'     labs(title = "Not rarefied")
#' p2 <- plotReducedDim(tse, "unifrac_rare", colour_by = "SampleType") +
#'     labs(title = "Rarefied")
#' p1 + p2
#'
#' @seealso
#' \code{\link[scater:runMDS]{scater::calculateMDS}} and
#' \code{\link[=getDissimilarity]{getDissimilarity}}
#'
#'
NULL

#' @rdname addMDS
#' @export
setMethod("addMDS", signature = c(x = "SingleCellExperiment"),
    function(x, name = "MDS", ...){
        if( !.is_a_string(name) ){
            stop("'name' must be a single character value.", call. = FALSE)
        }
        # Hiddenly support altExp
        x <- .check_and_get_altExp(x, ...)
        # Calculate indices
        args <- c(list(x = x), list(...))
        args <- args[ !names(args) %in% c("altexp") ]
        res <- do.call(getMDS, args)
        # Add object to reducedDim
        x <- .add_object_to_reduceddim(x, res, name = name, ...)
        return(x)
    }
)

#' @rdname addMDS
#' @export
#' @importFrom scater calculateMDS
setMethod("getMDS", signature = c(x = "SingleCellExperiment"),
    function(x, assay.type = "counts", ...){
        .check_assay_present(assay.type, x)
        args <- .get_mds_args(x, assay.type = assay.type, ...)
        res <- do.call(calculateMDS, args)
        return(res)
    }
)

#' @rdname addMDS
#' @export
#' @importFrom scater calculateMDS
setMethod("getMDS", signature = c(x = "TreeSummarizedExperiment"),
    function(x, assay.type = "counts", ...){
        .check_assay_present(assay.type, x)
        args <- .get_mds_args_treese(x, assay.type = assay.type, ...)
        res <- do.call(calculateMDS, args)
        return(res)
    }
)

################################ HELP FUNCTIONS ################################

# This function is used to set default options for SCE
.get_mds_args <- function(x, assay.type, FUN = getDissimilarity, ...){
    args <- c(list(x = x, assay.type = assay.type, FUN = FUN), list(...))
    return(args)
}

# For TreeSE, we also feed rowTree and node.labels as default
.get_mds_args_treese <- function(
        x, tree.name = "phylo", tree = NULL, node.label = NULL, ...){
    # Get tree and corresponding node.labels
    if( is.null(tree) ){
        tree <- rowTree(x, tree.name)
    }
    if( is.null(node.label) ){
        node.labels <- rowLinks(x)
        node.labels <- node.labels[
            node.labels[["whichTree"]] %in% tree.name, "nodeLab"]
    }
    #
    args <- c(.get_mds_args(x, ...), list(tree = tree, node.label = node.label))
    return(args)
}
