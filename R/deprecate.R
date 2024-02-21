#' These functions will be deprecated. Please use other functions instead.
#' 
#' @param x a \code{\link{SummarizedExperiment}} object -
#' 
#' @param ... -
#' 
#' @name deprecate
NULL

#' @rdname deprecate
setGeneric("addTaxonomyTree",
           signature = "x",
           function(x, ...)
               standardGeneric("addTaxonomyTree"))

#' @rdname deprecate
#' @export
setMethod("addTaxonomyTree", signature = c(x = "SummarizedExperiment"),
          function(x){
              .Deprecated(msg = paste0("'addTaxonomyTree' is deprecated.\n",
                                       "Use 'addTree' instead."))
              addTree(x)
          }
)

#' @rdname deprecate
setGeneric("taxonomyTree",
           signature = "x",
           function(x, ...)
               standardGeneric("taxonomyTree"))

#' @rdname deprecate
#' @export
setMethod("taxonomyTree", signature = c(x = "SummarizedExperiment"),
          function(x){
              .Deprecated(msg = paste0("'taxonomyTree' is deprecated.\n",
                                       "Use 'getHierarchyTree' instead."))
              getHierarchyTree(x)
          }
)
