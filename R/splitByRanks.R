#' Split a \code{SingleCellExperiment} be taxonomic ranks
#'
#' \code{splitByRanks}
#'
#' @param x \code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#'
#' @param ranks a character vector defining taxonomic ranks. Must all be values
#'   of \code{taxonomicRanks()} function.
#'
#' @param ... arguments passed to \code{agglomerateByRank} function for
#'   \code{SummarizedExperiment} objects and other functions.
#'   See \code{\link[=agglomerate-methods]{agglomerateByRank}} for more details.
#'
#' @return \code{x}, with objects of \code{x} agglomerated for selected ranks
#'   as \code{altExps}.
#'
#' @seealso
#' \code{\link[=merge-methods]{mergeRows}},
#' \code{\link[scuttle:sumCountsAcrossFeatures]{sumCountsAcrossFeatures}},
#' \code{\link[=agglomerate-methods]{agglomerateByRank}},
#' \code{\link[SingleCellExperiment:altExps]{altExps}}
#'
#' @name splitByRanks
#'
#' @examples
#' data(GlobalPatterns)
#' # print the available taxonomic ranks
#' taxonomyRanks(GlobalPatterns)
#'
#' altExps(GlobalPatterns) <- splitByRanks(GlobalPatterns)
#' altExps(GlobalPatterns)
#' altExp(GlobalPatterns,"Kingdom")
#' altExp(GlobalPatterns,"Species")
NULL

setGeneric("splitByRanks",
           signature = "x",
           function(x, ...)
               standardGeneric("splitByRanks"))

#' @rdname splitByRanks
#' @export
setMethod("splitByRanks", signature = c(x = "SingleCellExperiment"),
    function(x, ranks = taxonomyRanks(x), ...){
        # input check
        if(!.is_non_empty_character(ranks)){
            stop("'ranks' must be character vector.",
                 call. = FALSE)
        }
        .check_taxonomic_ranks(ranks,x)
        #
        args <- list(...)
        args[["strip_altexp"]] <- TRUE
        FUN <- function(rank){
            do.call(agglomerateByRank,
                    c(list(x = x, rank = rank), args))
        }
        ans <- lapply(ranks,FUN)
        names(ans) <- ranks
        SimpleList(ans)
    }
)
