#' Getter / setter for relative abundance data
#'
#' This function is being deprecated.
#' Please use \code{assay(x, "relabundance")} instead, which provides a more
#' flexible and robust way to access and modify relative abundance data stored
#' in the assay slot of a \code{\link[TreeSummarizedExperiment:TreeSummarizedExperiment-class]{TreeSummarizedExperiment}} object.
#'
#' @param x a \code{\link[TreeSummarizedExperiment:TreeSummarizedExperiment-class]{TreeSummarizedExperiment}} object
#' @param value a matrix to store as the \sQuote{relabundance} assay
#' @param ... optional arguments not used currently.
#'
#' @return
#' For \code{relabundance}, the matrix stored with the name \dQuote{relabundance}.
#'
#' @name relabundance
#'
#' @export
#'
#' @examples
#' data(GlobalPatterns)
#' # Calculates relative abundances
#' GlobalPatterns <- transformCounts(GlobalPatterns, method="relabundance")
#' # Fetches calculated relative abundances
#' # head(assay(GlobalPatterns, "relabundance"))
NULL

#' @rdname relabundance
setGeneric("relabundance", signature = c("x"),
           function(x, value) standardGeneric("relabundance"))

#' @rdname relabundance
setGeneric("relabundance<-", signature = c("x"),
           function(x, value) standardGeneric("relabundance<-"))

#' @rdname relabundance
#' @importFrom SummarizedExperiment assays
#' @export
setMethod("relabundance", signature = c(x = "SummarizedExperiment"),
    function(x){
        .Deprecated(msg = paste0("'relabundance' is deprecated\n",
                                 "Use 'assay(x, 'relabundance')' instead."))
        assays(x)[["relabundance"]]
    }
)

#' @rdname relabundance
#' @importFrom SummarizedExperiment assays<-
#' @export
setReplaceMethod("relabundance", signature = c(x = "SummarizedExperiment"),
    function(x, value){
        .Deprecated(msg = paste0("'relabundance' is deprecated\n",
                                 "Use 'assay(x, 'relabundance')' instead."))
        assays(x)[["relabundance"]] <- value
        x
    }
)
