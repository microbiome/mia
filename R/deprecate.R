#' These functions will be deprecated. Please use other functions instead.
#' 
#' @param x a \code{\link{SummarizedExperiment}} object -
#' 
#' @param obj an object
#' 
#' @param ... -
#' 
#' @name deprecate
NULL

#' @rdname deprecate
#' @export
setGeneric("cluster", signature = c("x"),
            function(x,...)
                standardGeneric("cluster"))

#' @rdname deprecate
#' @export
#' @importFrom bluster clusterRows
setMethod("cluster", signature = c(x = "SummarizedExperiment"),
            function(x,...){
                .Deprecated(msg = paste0("'cluster' is deprecated. ",
                                        "Use 'addCluster' instead."))
                addCluster(x,...)
            }
)

#' @rdname deprecate
#' @export
setGeneric("addTaxonomyTree",
            signature = "x",
            function(x, ...)
                standardGeneric("addTaxonomyTree"))

#' @rdname deprecate
#' @export
setMethod("addTaxonomyTree", signature = c(x = "SummarizedExperiment"),
            function(x,...){
                .Deprecated(msg = paste0("'addTaxonomyTree' is deprecated. ",
                                        "Use 'addHierarchyTree' instead."))
                addHierarchyTree(x,...)
            }
)

#' @rdname deprecate
#' @export
setGeneric("taxonomyTree",
            signature = "x",
            function(x, ...)
                standardGeneric("taxonomyTree"))

#' @rdname deprecate
#' @export
setMethod("taxonomyTree", signature = c(x = "SummarizedExperiment"),
            function(x,...){
                .Deprecated(msg = paste0("'taxonomyTree' is deprecated. ",
                                        "Use 'getHierarchyTree' instead."))
                getHierarchyTree(x,...)
            }
)

#' @rdname deprecate
#' @export

setGeneric("makeTreeSummarizedExperimentFromPhyloseq", signature = c("obj"),
            function(obj)
                standardGeneric("makeTreeSummarizedExperimentFromPhyloseq"))

#' @rdname deprecate
#' @export
setMethod("makeTreeSummarizedExperimentFromPhyloseq", 
            signature = c(obj = "ANY"), function(obj){
                .Deprecated(msg = paste0(
                    "'makeTreeSummarizedExperimentFromPhyloseq' is deprecated.",
                    " Use 'convert' instead."))
                convert(obj)
            }
)

#' @rdname deprecate
#' @export
setGeneric("makePhyloseqFromTreeSummarizedExperiment", signature = c("obj"),
            function(obj, ...)
                standardGeneric("makePhyloseqFromTreeSummarizedExperiment"))

#' @rdname deprecate
#' @export
setMethod("makePhyloseqFromTreeSummarizedExperiment", 
            signature = c(obj = "ANY"), function(obj, ...){
                .Deprecated(msg = paste0(
                    "'makePhyloseqFromTreeSummarizedExperiment' is deprecated.",
                    " Use 'convert' instead."))
                convert(obj, ...)
            }
)

#' @rdname deprecate
#' @export
setGeneric("makePhyloseqFromTreeSE", signature = c("obj"),
           function(obj, ...)
               standardGeneric("makePhyloseqFromTreeSE"))

#' @rdname deprecate
#' @export
setMethod("makePhyloseqFromTreeSE", 
          signature = c(obj = "ANY"), function(obj, ...){
              .Deprecated(msg = paste0(
                  "'makePhyloseqFromTreeSE' is deprecated.",
                  " Use 'convert' instead."))
              convert(obj, ...)
          }
)

loadFromBiom <- function(...) {
    .Deprecated(msg = paste0("'loadFromBiom' is deprecated.",
                            " Use 'importBIOM' instead."))
    importBiom(...)
}

#' @rdname deprecate
#' @export
loadFromQIIME2 <- function(...) {
    .Deprecated(msg = paste0("'loadFromQIIME2' is deprecated.",
                            " Use 'importQIIME2' instead."))
    importQIIME2(...)
}

#' @rdname deprecate
#' @export
readQZA <- function(...) {
    .Deprecated(msg = paste0("'readQZA' is deprecated.",
                            " Use 'importQZA' instead."))
    importQZA(...)
}

#' @rdname deprecate
#' @export
loadFromMothur <- function(...) {
    .Deprecated(msg = paste0("'loadFromMothur' is deprecated.",
                            " Use 'importMothur' instead."))
    importMothur(...)
}

#' @rdname deprecate
#' @export   
loadFromMetaphlan <- function(...) {
    .Deprecated(msg = paste0("'loadFromMetaphlan' is deprecated.",
                            " Use 'importMetaPhlAn' instead."))
    importMetaPhlAn(...)
}

#' @rdname deprecate
#' @export    
loadFromHumann <- function(...) {
    .Deprecated(msg = paste0("'loadFromHumann' is deprecated.",
                            " Use 'importHUMAnN' instead."))
    importHUMAnN(...)
}
