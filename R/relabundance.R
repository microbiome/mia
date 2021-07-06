#' Getter / setter for relative abundance data
#'
#' \code{relabundance} is a getter/setter for relative abundance stored in the
#' assay slot \sQuote{relabundance} of a
#' \code{\link[TreeSummarizedExperiment:TreeSummarizedExperiment-class]{TreeSummarizedExperiment}}
#' object. This is a shortcut function for \code{assay(x,"relabundance")}.
#'
#' @param x a
#'   \code{\link[TreeSummarizedExperiment:TreeSummarizedExperiment-class]{TreeSummarizedExperiment}}
#'   object
#' @param value a \code{matrix} to store as the the \sQuote{relabundance} assay
#' @param ... optional arguments not used currently.
#'
#' @return
#' For \code{relabundance} the matrix stored with the name
#' \dQuote{relabundance}.
#'
#'
#' @name relabundance
#'
#' @export
#'
#' @examples
#' data(GlobalPatterns)
#' # Calculates relative abundances
#' GlobalPatterns <- relAbundanceCounts(GlobalPatterns)
#' # Fetches calculated relative abundances
#' head(relabundance(GlobalPatterns))
NULL

#' @rdname relabundance
setGeneric("relabundance", signature = c("x"),
           function(x, ...) standardGeneric("relabundance"))
#' @rdname relabundance
setGeneric("relabundance<-", signature = c("x"),
           function(x, value) standardGeneric("relabundance<-"))

#' @rdname relabundance
#' @importFrom SummarizedExperiment assays
#' @export
setMethod("relabundance",signature = c(x = "SummarizedExperiment"),
    function(x){
        assays(x)[["relabundance"]]
    }
)

#' @rdname relabundance
#' @importFrom SummarizedExperiment assays<-
#' @export
setReplaceMethod("relabundance", signature = c(x = "SummarizedExperiment"),
    function(x, value){
        assays(x)[["relabundance"]] <- value
        x
    }
)

