#' @export
setGeneric("convert", signature = c("x"),
           function(x,...)
               standardGeneric("convert"))

#' @export
setMethod("convert", signature = c(x = "SummarizedExperiment"),
          function(x,...){
              .makePhyloseqFromSE(x,...)
          }
)

#' @export
setMethod("convert", signature = c(x = "TreeSummarizedExperiment"),
          function(x,...){
              .makePhyloseqFromTreeSE(x,...)
          }
)

#' @export
setMethod("convert", signature = c(x = "dada"),
          function(x,...){
              .makeTreeSEFromDADA2(x,...)
          }
)

#' @export
setMethod("convert", signature = c(x = "phyloseq"),
          function(x,...){
              .makeTreeSEFromPhyloseq(x,...)
          }
)

#' @export
setMethod("convert", signature = c(x = "biom"),
          function(x,...){
              .makeTreeSEFromBiom(x,...)
          }
)
