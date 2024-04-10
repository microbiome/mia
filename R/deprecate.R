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
                .Deprecated(msg = paste0("'addTaxonomyTree' is deprecated.",
                                        "Use 'addHierarchyTree' instead."))
                addHierarchyTree(x)
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
                .Deprecated(msg = paste0("'taxonomyTree' is deprecated.",
                                        "Use 'getHierarchyTree' instead."))
                getHierarchyTree(x)
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
                                        "Use agglomerateByVariable with ",
                                        "parameter MARGIN = 'rows' instead."))
                agglomerateByVariable(x, MARGIN = "rows", ...)
            }
)

#' @rdname deprecate
#' @export
setMethod("mergeRows", signature = c(x = "TreeSummarizedExperiment"),
          function(x, ...){
              .Deprecated(msg = paste0("'mergeRows' is deprecated. ",
                                       "Use agglomerateByVariable with ",
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
                                        "Use agglomerateByVariable with ", 
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
                                        "Use agglomerateByVariable with ",
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
                                        "Use agglomerateByVariable with ", 
                                        "parameter MARGIN = 'rows' instead."))
                agglomerateByVariable(x, MARGIN = "rows", ...)
            }
)

#' @rdname deprecate
#' @export
setMethod("mergeFeatures", signature = c(x = "TreeSummarizedExperiment"),
          function(x, ...){
              .Deprecated(msg = paste0("'mergeFeatures' is deprecated. ",
                                       "Use agglomerateByVariable with ",
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
                                        "Use agglomerateByVariable with ",
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
                                       "Use agglomerateByVariable with ", 
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
          function(x, rank = taxonomyRanks(x)[1], onRankOnly = FALSE, 
                   na.rm = FALSE, 
                   empty.fields = c(NA, "", " ", "\t", "-", "_"), ...){
              .Deprecated(msg = paste0("mergeFeaturesByRank is deprecated. ",
                                        "Use agglomerateByRank instead."))
              x <- agglomerateByRank(x, rank = rank, onRankOnly = onRankOnly,
                                     na.rm = na.rm,
                                     empty.fields = empty.fields, ...)
              x
          }
)

#' @rdname deprecate
#' @importFrom SingleCellExperiment altExp altExp<- altExps<-
#' @export
setMethod("mergeFeaturesByRank", signature = c(x = "SingleCellExperiment"),
          function(x, ..., altexp = NULL, strip_altexp = TRUE){
              .Deprecated(msg = paste0("mergeFeaturesByRank is deprecated. ",
                                        "Use agglomerateByRank instead."))
              x <- agglomerateByRank(x, ..., altexp = altexp, 
                                     strip_altexp = strip_altexp)
              x
          }
)

#' @rdname deprecate
#' @export
setMethod("mergeFeaturesByRank", signature = c(x = "TreeSummarizedExperiment"),
          function(x, ..., agglomerate.tree = FALSE){
              .Deprecated(msg = paste0("mergeFeaturesByRank is deprecated. ",
                                        "Use agglomerateByRank instead."))
              x <- agglomerateByRank(x, ..., 
                                     agglomerate.tree = agglomerate.tree)
              x
          }
)
