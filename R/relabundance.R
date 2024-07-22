#' Getter / setter for relative abundance data
#'
#' This function is being deprecated and will be removed in future releases.
#' Please use \code{assay(x, "relabundance")} instead, which provides a more
#' flexible and robust way to access and modify relative abundance data stored
#' in the assay slot of a \code{\link[TreeSummarizedExperiment:TreeSummarizedExperiment-class]{TreeSummarizedExperiment}} object.
#'
#' @inheritParams calculateJSD
#' 
#' @param value \code{Character vector}. A matrix to store as the \sQuote{relabundance} assay
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
#' GlobalPatterns <- transformAssay(GlobalPatterns, method="relabundance")
#' # Fetches calculated relative abundances
#' # head(assay(GlobalPatterns, "relabundance"))
NULL

#' Define a generic function for relabundance
#' The generic function will dispatch to different methods depending on the
#' class of its argument
#' @rdname relabundance
setGeneric("relabundance", signature = c("x"),
           function(x, ...) standardGeneric("relabundance"))
           
#' Define a generic replacement function for relabundance
#' The generic function will dispatch to different methods depending on the
#' class of its argument

#' @rdname relabundance
setGeneric("relabundance<-", signature = c("x"),
           function(x, value) standardGeneric("relabundance<-"))
           

#' Define a method for relabundance for SummarizedExperiment objects
#' This method retrieves the relabundance data from the assay slot of the object
#' and issues a deprecation warning
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
#' Define a replacement method for relabundance for SummarizedExperiment objects
#' This method sets the relabundance data in the assay slot of the object to the
#' provided value and issues a deprecation warning
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
