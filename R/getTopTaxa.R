#' Returns top N \dQuote{FeatureID}
#'
#' This function is will extract top user specified number of \dQuote{FeatureID}
#' in a \code{\link[=SummarizedExperiment-class]{SummarizedExperiment}} object.
#'
#' @param x A
#'  \code{\link[=SummarizedExperiment-class]{SummarizedExperiment}} object.
#'
#' @param top Numeric value, how many top taxa to return. Default return top five taxa.
#'
#' @param method Specify the method to determine top taxa. Either sum, mean or median.
#'               Default is 'mean'
#'
#' @param abund_values a \code{character} value to select an
#'   \code{\link[SummarizedExperiment:SummarizedExperiment-class]{assayNames}}
#'
#' @return A vector of the top N \dQuote{FeatureID}
#' @name getTopTaxa
#'
#' @author
#' Sudarshan A. Shetty
#'
#' @examples
#' data(GlobalPatterns, package = "MicrobiomeExperiment")
#'
#' top_taxa <- getTopTaxa(GlobalPatterns,
#'                        method="mean",
#'                        top=5,
#'                        abund_values="counts")
#' top_taxa
#'
NULL

#' @rdname getTopTaxa
#'
#' @export
setGeneric("getTopTaxa", signature = "x",
           function(x, top=5L, method=c("mean","sum","median"), abund_values = "counts")
               standardGeneric("getTopTaxa"))

.check_max_taxa <- function(x, top, abund_values){
    if(!is.numeric(top) || as.integer(top) != top){
        top("'top' must be integer value", call. = FALSE)
    }
    if(top > nrow(assay(x,abund_values))){
        stop("'top' must be <= nrow(x)", call. = FALSE)
    }
}

#' @rdname getTopTaxa
#'
#' @importFrom DelayedMatrixStats rowSums2 rowMeans2 rowMedians
#'
#' @export
#'

setMethod("getTopTaxa", signature = c(x = "SummarizedExperiment"),
          function(x, top = 5L,
                   method=c("mean","sum","median"),
                   abund_values = "counts"){
              method <- match.arg(method, c("mean","sum","median"))
              # check max taxa
              .check_max_taxa(x, top, abund_values)
              # check assay
              .check_abund_values(abund_values, x)
              #assay(x, abund_values)[,sample_id]
              taxs <- switch(method,
                             mean = rowMeans2(assay(x, abund_values)),
                             sum = rowSums2(assay(x, abund_values)),
                             median = rowMedians(assay(x, abund_values)))
              names(taxs) <- rownames(assay(x))
              taxs <- sort(taxs,decreasing = TRUE)[1:top]
              return(names(taxs))
          })

