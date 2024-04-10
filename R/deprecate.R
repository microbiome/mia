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
#' @aliases mergeFeatures
#' @export
setGeneric("mergeRows",
            signature = "x",
            function(x, ...)
                standardGeneric("mergeRows"))

#' @rdname deprecate
#' @aliases mergeFeatures
#' @export
setMethod("mergeRows", signature = c(x = "SummarizedExperiment"),
            function(x, ...){
                .Deprecated(msg = paste0("'mergeRows' is deprecated. ",
                                        "Use 'agglomerateByVariable' with ",
                                        "parameter MARGIN = 'rows' instead."))
                agglomerateByVariable(x, MARGIN = "rows", ...)
            }
)

#' @rdname deprecate
#' @export
setMethod("mergeRows", signature = c(x = "TreeSummarizedExperiment"),
          function(x, ...){
              .Deprecated(msg = paste0("'mergeRows' is deprecated. ",
                                       "Use 'agglomerateByVariable' with ",
                                        "parameter MARGIN = 'rows' instead."))
              agglomerateByVariable(x, MARGIN = "rows", ...)
          }
)

#' @rdname deprecate
#' @aliases mergeSamples
#' @export
setGeneric("mergeCols",
            signature = "x",
            function(x, ...)
                standardGeneric("mergeCols"))

#' @rdname deprecate
#' @aliases mergeSamples
#' @export
setMethod("mergeCols", signature = c(x = "SummarizedExperiment"),
            function(x, ...){
                .Deprecated(msg = paste0("'mergeCols' is deprecated. ",
                                        "Use 'agglomerateByVariable' with ", 
                                        "parameter MARGIN = 'cols' instead."))
                agglomerateByVariable(x, MARGIN = "cols", ...)
            }
)

#' @rdname deprecate
#' @aliases mergeSamples
#' @export
setMethod("mergeCols", signature = c(x = "TreeSummarizedExperiment"),
          function(x, ...){
              .Deprecated(msg = paste0("'mergeCols' is deprecated. ",
                                        "Use 'agglomerateByVariable' with ",
                                        "parameter MARGIN = 'cols' instead."))
              agglomerateByVariable(x, MARGIN = "cols", ...)
          }
)

#' @rdname deprecate
#' @aliases mergeRows
#' @export
setGeneric("mergeFeatures",
            signature = "x",
            function(x, ...)
                standardGeneric("mergeFeatures"))

#' @rdname deprecate
#' @export
setMethod("mergeFeatures", signature = c(x = "SummarizedExperiment"),
            function(x, ...){
                .Deprecated(msg = paste0("'mergeFeatures' is deprecated. ",
                                        "Use 'agglomerateByVariable' with ", 
                                        "parameter MARGIN = 'rows' instead."))
                agglomerateByVariable(x, MARGIN = "rows", ...)
            }
)

#' @rdname deprecate
#' @export
setMethod("mergeFeatures", signature = c(x = "TreeSummarizedExperiment"),
          function(x, ...){
              .Deprecated(msg = paste0("'mergeFeatures' is deprecated. ",
                                       "Use 'agglomerateByVariable' with ",
                                        "parameter MARGIN = 'rows' instead."))
              agglomerateByVariable(x, MARGIN = "rows", ...)
          }
)

#' @rdname deprecate
#' @aliases mergeCols
#' @export
setGeneric("mergeSamples",
            signature = "x",
            function(x, ...)
                standardGeneric("mergeSamples"))

#' @rdname deprecate
#' @aliases mergeSamples
#' @export
setMethod("mergeSamples", signature = c(x = "SummarizedExperiment"),
            function(x, ...){
                .Deprecated(msg = paste0("'mergeSamples' is deprecated. ",
                                        "Use 'agglomerateByVariable' with ",
                                        "parameter MARGIN = 'cols' instead."))
                agglomerateByVariable(x, MARGIN = "cols", ...)
            }
)

#' @rdname deprecate
#' @aliases mergeSamples
#' @export
setMethod("mergeSamples", signature = c(x = "TreeSummarizedExperiment"),
          function(x, ...){
              .Deprecated(msg = paste0("'mergeSamples' is deprecated. ",
                                       "Use 'agglomerateByVariable' with ", 
                                        "parameter MARGIN = 'cols' instead."))
              agglomerateByVariable(x, MARGIN = "cols", ...)
          }
)

#' @rdname deprecate
#' @export
setGeneric("mergeFeaturesByRank",
           signature = "x",
           function(x, ...)
               standardGeneric("mergeFeaturesByRank"))

#' @rdname deprecate
#' @importFrom SummarizedExperiment rowData rowData<-
#' @export
setMethod("mergeFeaturesByRank", signature = c(x = "SummarizedExperiment"),
          function(x, ...){
              .Deprecated(msg = paste0("'mergeFeaturesByRank' is deprecated. ",
                                        "Use 'agglomerateByRank' instead."))
              x <- agglomerateByRank(x, ...)
              x
          }
)

#' @rdname deprecate
#' @importFrom SingleCellExperiment altExp altExp<- altExps<-
#' @export
setMethod("mergeFeaturesByRank", signature = c(x = "SingleCellExperiment"),
          function(x, ...){
              .Deprecated(msg = paste0("'mergeFeaturesByRank' is deprecated. ",
                                        "Use 'agglomerateByRank' instead."))
              x <- agglomerateByRank(x, ...)
              x
          }
)

#' @rdname deprecate
#' @export
setMethod("mergeFeaturesByRank", signature = c(x = "TreeSummarizedExperiment"),
          function(x, ...){
              .Deprecated(msg = paste0("'mergeFeaturesByRank' is deprecated. ",
                                        "Use 'agglomerateByRank' instead."))
              x <- agglomerateByRank(x, ...)
              x
          }
)
