#' @name
#' addPCA
#'
#' @title
#' Perform principal component analysis (PCA)
#'
#' @description
#' Perform principal component analysis (PCA). These functions are wrappers
#' for \code{\link[scater:runPCA]{scater::calculatePCA}}.
#'
#' @details
#' These functions are wrappers for
#' \code{\link[scater:runPCA]{scater::calculatePCA}} and
#' \code{\link[scater:runPCA]{scater::runPCA}}. While \code{getPCA}
#' returns the results, \code{addPCA} adds them to \code{reducedDim(x)}.
#'
#' Unlike \code{scater::calculatePCA}, which uses the 500 most variable
#' features by default, these functions use all features (\code{ntop = Inf}).
#'
#' See \code{\link[scater:runPCA]{scater::calculatePCA}} for details.
#'
#' @return
#' \code{getPCA} returns PCA results.
#' \code{addPCA} returns \code{x} with PCA results added to its
#' \code{reducedDim(x, name)}.
#'
#' @inheritParams addMDS
#'
#' @param name \code{Character scalar}. A name for the \code{reducedDim()}
#' where results will be stored. (Default: \code{"PCA"})
#'
#' @param ... additional arguments passed to
#' \code{\link[scater:runPCA]{scater::calculatePCA}}.
#' \itemize{
#'   \item \code{ntop}: \code{Numeric scalar}. Number of the most variable
#'   features used. (Default: \code{Inf})
#' }
#'
#' @examples
#' library(mia)
#' library(scater)
#'
#' data(GlobalPatterns)
#' tse <- GlobalPatterns
#'
#' # Calculate PCA from CLR-transformed abundances
#' tse <- transformAssay(tse, method = "clr", pseudocount = TRUE)
#' tse <- addPCA(tse, assay.type = "clr")
#'
#' # Use only the 500 most variable features
#' tse <- addPCA(tse, assay.type = "clr", ntop = 500, name = "PCA_top500")
#'
#' plotReducedDim(tse, "PCA", colour_by = "SampleType")
#'
#' @seealso
#' \code{\link[scater:runPCA]{scater::calculatePCA}}
#'
NULL

#' @rdname addPCA
#' @export
setMethod("addPCA", signature = c(x = "SingleCellExperiment"),
    function(x, name = "PCA", ...){
        if( !.is_a_string(name) ){
            stop("'name' must be a single character value.", call. = FALSE)
        }
        # Hiddenly support altExp
        x <- .check_and_get_altExp(x, ...)
        # Calculate PCA
        args <- c(list(x = x), list(...))
        args <- args[ !names(args) %in% c("altexp") ]
        res <- do.call(getPCA, args)
        # Add object to reducedDim
        x <- .add_object_to_reduceddim(x, res, name = name, ...)
        return(x)
    }
)

#' @rdname addPCA
#' @export
#' @importFrom scater calculatePCA
setMethod("getPCA", signature = c(x = "SingleCellExperiment"),
    function(x, assay.type = "counts", ...){
        .check_assay_present(assay.type, x)
        args <- c(list(x = x, assay.type = assay.type), list(...))
        # Use all features by default
        if( !"ntop" %in% names(args) ){
            args$ntop <- Inf
        }
        res <- do.call(calculatePCA, args)
        return(res)
    }
)
