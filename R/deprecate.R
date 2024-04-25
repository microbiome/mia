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
setGeneric("getUniqueFeatures",
            signature = c("x"),
            function(x, ...)
                standardGeneric("getUniqueFeatures"))

#' @rdname deprecate
#' @export
setMethod("getUniqueFeatures", signature = c(x = "SummarizedExperiment"),
            function(x,...){
                .Deprecated(msg = paste0("'getUniqueFeatures' is deprecated. ",
                                        "Use 'getUnique' instead."))
                getUnique(x,...)
            }
)

#' @rdname deprecate
#' @export
setGeneric("getUniqueTaxa",
            signature = c("x"),
            function(x, ...)
                standardGeneric("getUniqueTaxa"))

#' @rdname deprecate
#' @export
setMethod("getUniqueTaxa", signature = c(x = "SummarizedExperiment"),
            function(x,...){
                .Deprecated(msg = paste0("'getUniqueTaxa' is deprecated. ",
                                        "Use 'getUnique' instead."))
                getUnique(x,...)
            }
)

#' @rdname deprecate
#' @export
setGeneric("getTopFeatures", signature = "x",
            function(x,...)
                standardGeneric("getTopFeatures"))

#' @rdname deprecate
#' @export
setMethod("getTopFeatures", signature = c(x = "SummarizedExperiment"),
            function(x,...){
                .Deprecated(msg = paste0("'getTopFeatures' is deprecated. ",
                                        "Use 'getTop' instead."))
                getTop(x,...)
            }
)

#' @rdname deprecate
#' @export
setGeneric("getTopTaxa", signature = "x",
            function(x,...)
                standardGeneric("getTopTaxa"))

#' @rdname deprecate
#' @export
setMethod("getTopTaxa", signature = c(x = "SummarizedExperiment"),
            function(x,...){
                .Deprecated(msg = paste0("'getTopTaxa' is deprecated. ",
                                        "Use 'getTop' instead."))
                getTop(x,...)
            }
)

#' @rdname deprecate
#' @export
setGeneric("getRareFeatures", signature = "x",
            function(x, ...)
                standardGeneric("getRareFeatures"))

#' @rdname deprecate
#' @export
setMethod("getRareFeatures", signature = c(x = "ANY"),
            function(x,...){
                .Deprecated(msg = paste0("'getRareFeatures' is deprecated. ",
                                        "Use 'getRare' instead."))
                getRare(x,...)
            }
)

#' @rdname deprecate
#' @export
setMethod("getRareFeatures", signature = c(x = "SummarizedExperiment"),
            function(x,...){
                .Deprecated(msg = paste0("'getRareFeatures' is deprecated. ",
                                        "Use 'getRare' instead."))
                getRare(x,...)
            }
)

#' @rdname deprecate
#' @export
setGeneric("getRareTaxa", signature = "x",
            function(x, ...)
                standardGeneric("getRareTaxa"))

#' @rdname deprecate
#' @export
setMethod("getRareTaxa", signature = c(x = "ANY"),
            function(x,...){
                .Deprecated(msg = paste0("'getRareTaxa' is deprecated. ",
                                        "Use 'getRare' instead."))
                getRare(x,...)
            }
)

#' @rdname deprecate
#' @export
setMethod("getRareTaxa", signature = c(x = "SummarizedExperiment"),
            function(x,...){
                .Deprecated(msg = paste0("'getRareTaxa' is deprecated. ",
                                        "Use 'getRare' instead."))
                getRare(x,...)
            }
)

#' @rdname deprecate
#' @export
setGeneric("getPrevalentFeatures", signature = "x",
            function(x, ...)
                standardGeneric("getPrevalentFeatures"))

#' @rdname deprecate
#' @export
setMethod("getPrevalentFeatures", signature = c(x = "ANY"),
            function(x,...){
                .Deprecated(msg = paste0("'getPrevalentFeatures' is deprecated. ",
                                        "Use 'getPrevalent' instead."))
                getPrevalent(x,...)
            }
          
)

#' @rdname deprecate
#' @export
setMethod("getPrevalentFeatures", signature = c(x = "SummarizedExperiment"),
            function(x,...){
                .Deprecated(msg = paste0("'getPrevalentFeatures' is deprecated. ",
                                        "Use 'getPrevalent' instead."))
                getPrevalent(x,...)
            }
)

#' @rdname deprecate
#' @export
setGeneric("getPrevalentTaxa", signature = "x",
            function(x, ...)
                standardGeneric("getPrevalentTaxa"))

#' @rdname deprecate
#' @export
setMethod("getPrevalentTaxa", signature = c(x = "ANY"),
            function(x,...){
                .Deprecated(msg = paste0("'getPrevalentTaxa' is deprecated. ",
                                        "Use 'getPrevalent' instead."))
                getPrevalent(x,...)
            }
           
)

#' @rdname deprecate
#' @export
setMethod("getPrevalentTaxa", signature = c(x = "SummarizedExperiment"),
            function(x,...){
                .Deprecated(msg = paste0("'getPrevalentTaxa' is deprecated. ",
                                        "Use 'getPrevalent' instead."))
                getPrevalent(x,...)
            }
)

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
