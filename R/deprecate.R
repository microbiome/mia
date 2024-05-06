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
setGeneric("mergeRows",
            signature = "x",
            function(x, ...)
                standardGeneric("mergeRows"))

#' @rdname deprecate
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
#' @export
setGeneric("mergeCols",
            signature = "x",
            function(x, ...)
                standardGeneric("mergeCols"))

#' @rdname deprecate
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
#' @export
setGeneric("mergeSamples",
            signature = "x",
            function(x, ...)
                standardGeneric("mergeSamples"))

#' @rdname deprecate
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

#' @rdname deprecate
#' @export
loadFromBiom <- function(...) {
    .Deprecated(msg = paste0("'loadFromBiom' is deprecated.",
                            " Use 'importBIOM' instead."))
    importBIOM(...)
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

#' @rdname deprecate
#' @export
setGeneric("full_join", signature = c("x"),
           function(x, ...)
               standardGeneric("full_join"))

#' @rdname deprecate
#' @export
setMethod("full_join", signature = c(x = "ANY"),
            function(x, ...){
                .Deprecated(msg = paste0("'full_join' is deprecated. ",
                                        "Use 'mergeSEs' with 'join = full' ",
                                        "instead."))
                mergeSEs(x, join = "full", ...)
            }
)

#' @rdname deprecate
#' @export
setGeneric("inner_join", signature = c("x"),
           function(x, ...)
               standardGeneric("inner_join"))

#' @rdname deprecate
#' @export
setMethod("inner_join", signature = c(x = "ANY"),
            function(x, ...){
                .Deprecated(msg = paste0("'inner_join' is deprecated. ",
                                        "Use 'mergeSEs' with 'join = inner' ",
                                        "instead."))
                mergeSEs(x, join = "inner", ...)
            }
)

#' @rdname deprecate
#' @export
setGeneric("left_join", signature = c("x"),
           function(x, ...)
               standardGeneric("left_join"))

#' @rdname deprecate
#' @export
setMethod("left_join", signature = c(x = "ANY"),
          function(x, ...){
              .Deprecated(msg = paste0("'left_join' is deprecated. ",
                                        "Use 'mergeSEs' with 'join = left' ",
                                        "instead."))
              mergeSEs(x, join = "left", ...)
          }
)

#' @rdname deprecate
#' @export
setGeneric("right_join", signature = c("x"),
           function(x, ...)
               standardGeneric("right_join"))

#' @rdname deprecate
#' @export
setMethod("right_join", signature = c(x = "ANY"),
            function(x, ...){
                .Deprecated(msg = paste0("'right_join' is deprecated. ",
                                        "Use 'mergeSEs' with 'join = right' ",
                                        "instead."))
                mergeSEs(x, join = "right", ...)
            }
)

#' @rdname deprecate
#' @export    
plotNMDS <- function(x, ...){
    .Deprecated(msg = paste0("'plotNMDS' is deprecated. ",
                             "Use 'scater::plotReducedDim' with ",
                             "dimred = 'NMDS' instead."))
    plotReducedDim(x, ncomponents = 2, dimred = "NMDS",...)
}

#' @rdname deprecate
#' @export
setGeneric("estimateDivergence",signature = c("x"),
            function(x, ...)
                standardGeneric("estimateDivergence"))

#' @rdname deprecate
#' @export
setMethod("estimateDivergence", signature = c(x="SummarizedExperiment"),
            function(x, ...){
                .Deprecated(msg = paste0("'estimateDivergence' is deprecated. ",
                                        "Use 'addDivergence' instead."))
                addDivergence(x, ...)
            }
)
