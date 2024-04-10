#' These functions will be deprecated. Please use other functions instead.
#' 
#' @param x A
#'   \code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#'   object.
#'   
#' @param ... Additional parameters. See dedicated function.
#' 
#' @name deprecate
NULL

#' @rdname deprecate
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
setGeneric("mergeFeaturesByPrevalence", signature = "x",
           function(x, ...)
               standardGeneric("mergeFeaturesByPrevalence"))

#' @rdname deprecate
#' @export
setMethod("mergeFeaturesByPrevalence", signature = c(x = "SummarizedExperiment"),
          function(x, ...){
              .Deprecated(msg = paste0(
                  "'mergeFeaturesByPrevalence' is deprecated. ",
                  "Use agglomerateByPrevalence instead."))
              x <- agglomerateByPrevalence(x, ...)
              x 
          }
)
