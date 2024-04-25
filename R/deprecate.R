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

#' @importFrom bluster clusterRows
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
setGeneric(
    "estimateEvenness", signature = c("x"),
    function(x, ...) standardGeneric("estimateEvenness"))

#' @rdname deprecate
#' @export
setMethod(
    "estimateEvenness", signature = c(x="ANY"),
    function(x, ...){
        .Deprecated(msg = paste0("'estimateEvenness' is deprecated. ",
                                "Use 'addAlpha' instead."))
        .estimate_evenness(x, ...)
    }
)

#' @rdname deprecate
#' @export
setGeneric(
    "estimateRichness", signature = c("x"),
    function(x, ...) standardGeneric("estimateRichness"))

#' @rdname deprecate
#' @export
setMethod(
    "estimateRichness", signature = c(x="ANY"),
    function(x, ...){
        .Deprecated(msg = paste0("'estimateRichness' is deprecated. ",
                                "Use 'addAlpha' instead."))
        .estimate_richness(x, ...)
    }
)

#' @rdname deprecate
#' @export
setGeneric(
    "estimateDiversity", signature = c("x"),
    function(x, ...) standardGeneric("estimateDiversity"))

#' @rdname deprecate
#' @export
setMethod(
    "estimateDiversity", signature = c(x="ANY"),
    function(x, ...){
        .Deprecated(msg = paste0("'estimateDiversity' is deprecated. ",
                                "Use 'addAlpha' instead."))
        .estimate_diversity(x, ...)
    }
)

#' @rdname deprecate
#' @export
setGeneric(
    "estimateFaith", signature = c("x"),
    function(x, ...) standardGeneric("estimateFaith"))

#' @rdname deprecate
#' @export
setMethod(
    "estimateFaith", signature = c(x="ANY"),
    function(x, ...){
        .Deprecated(msg = paste0("'estimateFaith' is deprecated. ",
                                "Use 'addAlpha' instead."))
        .estimate_faith(x, ...)
    }
)

#' @rdname deprecate
#' @export
setGeneric(
    "estimateDominance", signature = c("x"),
    function(x, ...) standardGeneric("estimateDominance"))

#' @rdname deprecate
#' @export
setMethod(
    "estimateDominance", signature = c(x="ANY"),
    function(x, ...){
        .Deprecated(msg = paste0("'estimateDominance' is deprecated. ",
                                "Use 'addAlpha' instead."))
        .estimate_dominance(x, ...)
    })
