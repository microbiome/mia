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
            function(x, f, archetype = 1L, ...)
                standardGeneric("mergeRows"))

#' @rdname deprecate
#' @aliases mergeFeatures
#' @export
setMethod("mergeRows", signature = c(x = "SummarizedExperiment"),
            function(x, f, archetype = 1L, ...){
                .Deprecated(msg = paste0("'mergeRows' is deprecated.",
                                        "Use AgglomerateByVariable with 
                                        parameter MARGIN = 'rows' instead"))
                AgglomerateByVariable(x, MARGIN = 'rows', f, 
                                    archetype = archetype, ...)
            }
)

#' @rdname deprecate
#' @export
setMethod("mergeRows", signature = c(x = "TreeSummarizedExperiment"),
          function(x, f, archetype = 1L, mergeTree = FALSE, 
                   mergeRefSeq = FALSE, ...){
              .Deprecated(msg = paste0("'mergeRows' is deprecated.",
                                       "Use AgglomerateByVariable with 
                                        parameter MARGIN = 'rows' instead"))
              AgglomerateByVariable(x, MARGIN = 'rows', f, 
                                    archetype = archetype,
                                    mergeTree = mergeTree,
                                    mergeRefSeq = mergeRefSeq,
                                    ...)
          }
)

#' @rdname deprecate
#' @aliases mergeSamples
#' @export
setGeneric("mergeCols",
            signature = "x",
            function(x, f, archetype = 1L, ...)
                standardGeneric("mergeCols"))

#' @rdname deprecate
#' @aliases mergeSamples
#' @export
setMethod("mergeCols", signature = c(x = "SummarizedExperiment"),
            function(x, f, archetype = 1L, ...){
                .Deprecated(msg = paste0("'mergeCols' is deprecated.",
                                        "Use AgglomerateByVariable with 
                                        parameter MARGIN = 'cols' instead"))
                AgglomerateByVariable(x, MARGIN = "cols", f, 
                                    archetype = archetype, ...)
            }
)

#' @rdname deprecate
#' @aliases mergeSamples
#' @export
setMethod("mergeCols", signature = c(x = "TreeSummarizedExperiment"),
          function(x, f, archetype = 1L, mergeTree = FALSE, 
                   mergeRefSeq = FALSE, ...){
              .Deprecated(msg = paste0("'mergeCols' is deprecated.",
                                       "Use AgglomerateByVariable with 
                                        parameter MARGIN = 'cols' instead"))
              AgglomerateByVariable(x, MARGIN = 'cols', f, 
                                    archetype = archetype,
                                    mergeTree = mergeTree,
                                    mergeRefSeq = mergeRefSeq,
                                    ...)
          }
)

#' @rdname deprecate
#' @aliases mergeRows
#' @export
setGeneric("mergeFeatures",
            signature = "x",
            function(x, f, archetype = 1L, ...)
                standardGeneric("mergeFeatures"))

#' @rdname deprecate
#' @export
setMethod("mergeFeatures", signature = c(x = "SummarizedExperiment"),
            function(x, f, archetype = 1L, ...){
                .Deprecated(msg = paste0("'mergeFeatures' is deprecated.",
                                        "Use AgglomerateByVariable with 
                                        parameter MARGIN = 'rows' instead"))
                AgglomerateByVariable(x, MARGIN = 'rows', f, 
                                    archetype = archetype, ...)
            }
)

#' @rdname deprecate
#' @export
setMethod("mergeFeatures", signature = c(x = "TreeSummarizedExperiment"),
          function(x, f, archetype = 1L, mergeTree = FALSE, 
                   mergeRefSeq = FALSE, ...){
              .Deprecated(msg = paste0("'mergeFeatures' is deprecated.",
                                       "Use AgglomerateByVariable with 
                                        parameter MARGIN = 'rows' instead"))
              AgglomerateByVariable(x, MARGIN = 'rows', f, 
                                    archetype = archetype,
                                    mergeTree = mergeTree,
                                    mergeRefSeq = mergeRefSeq,
                                    ...)
          }
)

#' @rdname deprecate
#' @aliases mergeCols
#' @export
setGeneric("mergeSamples",
            signature = "x",
            function(x, f, archetype = 1L, ...)
                standardGeneric("mergeSamples"))

#' @rdname deprecate
#' @aliases mergeSamples
#' @export
setMethod("mergeSamples", signature = c(x = "SummarizedExperiment"),
            function(x, f, archetype = 1L, ...){
                .Deprecated(msg = paste0("'mergeSamples' is deprecated.",
                                        "Use AgglomerateByVariable with 
                                        parameter MARGIN = 'cols' instead"))
                AgglomerateByVariable(x, MARGIN = "cols", f, 
                                    archetype = archetype, ...)
            }
)

#' @rdname deprecate
#' @aliases mergeSamples
#' @export
setMethod("mergeSamples", signature = c(x = "TreeSummarizedExperiment"),
          function(x, f, archetype = 1L, mergeTree = FALSE, 
                   mergeRefSeq = FALSE, ...){
              .Deprecated(msg = paste0("'mergeSamples' is deprecated.",
                                       "Use AgglomerateByVariable with 
                                        parameter MARGIN = 'cols' instead"))
              AgglomerateByVariable(x, MARGIN = 'cols', f, 
                                    archetype = archetype,
                                    mergeTree = mergeTree,
                                    mergeRefSeq = mergeRefSeq,
                                    ...)
          }
)
