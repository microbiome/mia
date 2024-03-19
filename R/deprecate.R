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
setGeneric("addTaxonomyTree",
            signature = "x",
            function(x, ...)
                standardGeneric("addTaxonomyTree"))

#' @rdname deprecate
#' @export
setMethod("addTaxonomyTree", signature = c(x = "SummarizedExperiment"),
            function(x){
                .Deprecated(msg = paste0("'addTaxonomyTree' is deprecated.",
                                        "Use 'addHierarchyTree' instead."))
                addHierarchyTree(x)
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
            function(x){
                .Deprecated(msg = paste0("'taxonomyTree' is deprecated.",
                                        " Use 'getHierarchyTree' instead."))
                getHierarchyTree(x)
            }
)

#' @rdname deprecate
#' @export
setGeneric("makeTreeSEFromPhyloseq", signature = c("obj"),
           function(obj)
               standardGeneric("makeTreeSEFromPhyloseq"))

#' @rdname deprecate
#' @export
setMethod("makeTreeSEFromPhyloseq", signature = c(obj = "ANY"),
            function(obj){
                .Deprecated(msg = paste0(
                    "'makeTreeSEFromPhyloseq' is deprecated.",
                    " Use 'convert' instead."))
                convert(obj)
            }
)

#' @rdname deprecate
#' @export
setGeneric("makeTreeSEFromDADA2", signature = c("obj"),
            function(obj,...)
                standardGeneric("makeTreeSEFromDADA2"))

#' @rdname deprecate
#' @export
setMethod("makeTreeSEFromDADA2", signature = c(obj = "ANY"),
            function(obj,...){
                .Deprecated(msg = paste0(
                    "'makeTreeSEFromDADA2' is deprecated.",
                    " Use 'convert' instead."))
                convert(obj,...)
            }
)

#' @rdname deprecate
#' @export
setGeneric("makePhyloseqFromTreeSE", signature = c("obj"),
           function(obj, ...)
               standardGeneric("makePhyloseqFromTreeSE"))

#' @rdname deprecate
#' @export
setMethod("makePhyloseqFromTreeSE", signature = c(obj = "ANY"),
          function(obj, ...){
              .Deprecated(msg = paste0(
                  "'makePhyloseqFromTreeSE' is deprecated.",
                  " Use 'convert' instead."))
              convert(obj, ...)
          })

#' @rdname deprecate
#' @export
setGeneric("makeTreeSEFromBiom", signature = c("obj"),
           function(obj, ...)
               standardGeneric("makeTreeSEFromBiom"))

#' @rdname deprecate
#' @export
setMethod("makeTreeSEFromBiom", signature = c(obj = "ANY"),
          function(obj, ...){
              .Deprecated(msg = paste0(
                  "'makeTreeSEFromBiom' is deprecated.",
                  " Use 'convert' instead."))
              convert(obj, ...)
          })
