#' These functions will be deprecated. Please use other functions instead.
#' 
#' @param x A
#'   \code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#'   object.
#'   
#' @param y A
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

#' @rdname deprecate
#' @export
setGeneric("calculateDPCoA", signature = c("x", "y"),
            function(x,y,...)
                standardGeneric("calculateDPCoA"))

#' @rdname deprecate
#' @export
setMethod("calculateDPCoA", c("ANY","ANY"), 
            function(x,y,...) {
                .Deprecated(msg = paste0("'calculateDPCoA' is deprecated. ",
                                        "Use 'getDPCoA' instead."))
                getDPCoA(x,y,...)
            }
)

#' @rdname deprecate
#' @importFrom ape cophenetic.phylo
#' @export
setMethod("calculateDPCoA", signature = c("TreeSummarizedExperiment","missing"),
            function(x,y,...) {
                .Deprecated(msg = paste0("'calculateDPCoA' is deprecated. ",
                                        "Use 'getDPCoA' instead."))
                getDPCoA(x,y,...)
            }
)

#' @rdname deprecate
#' @export
setGeneric("calculateNMDS", function(x, ...) standardGeneric("calculateNMDS"))

#' @rdname deprecate
#' @export
setMethod("calculateNMDS", "ANY", 
            function(x,...) {
                .Deprecated(msg = paste0("'calculateNMDS' is deprecated. ",
                                        "Use 'getNMDS' instead."))
                getNMDS(x,...)
            }
)

#' @rdname deprecate
#' @importFrom SummarizedExperiment assay
#' @export
setMethod("calculateNMDS", "SummarizedExperiment",
            function(x, ...) {
                .Deprecated(msg = paste0("'calculateNMDS' is deprecated. ",
                                        "Use 'getNMDS' instead."))
                getNMDS(x,...)
            }
)

#' @rdname deprecate
#' @export
setMethod("calculateNMDS", "SingleCellExperiment",
            function(x, ...){
                .Deprecated(msg = paste0("'calculateNMDS' is deprecated. ",
                                        "Use 'getNMDS' instead."))
                getNMDS(x,...)
            }
)

#' @rdname deprecate
#' @export 
setGeneric("calculateRDA", signature = c("x"),
            function(x, ...)
                standardGeneric("calculateRDA"))

#' @rdname deprecate
#' @export
setMethod("calculateRDA", "ANY",
            function(x, ...){
                .Deprecated(msg = paste0("'calculateRDA' is deprecated. ",
                                        "Use 'getRDA' instead."))
                getRDA(x,...)
            }
)

#' @rdname deprecate
#' @export
setMethod("calculateRDA", "SummarizedExperiment",
            function(x,...) {
                .Deprecated(msg = paste0("'calculateRDA' is deprecated. ",
                                        "Use 'getRDA' instead."))
                getRDA(x,...)
            }
)

#' @rdname deprecate
#' @export
setGeneric("calculateCCA", signature = c("x"),
            function(x,...)
                standardGeneric("calculateCCA"))

#' @rdname deprecate
#' @export
setMethod("calculateCCA", "ANY",
            function(x,...){
                .Deprecated(msg = paste0("'calculateCCA' is deprecated. ",
                                        "Use 'getCCA' instead."))
                getCCA(x,...)
            }
)

#' @rdname deprecate
#' @export
setMethod("calculateCCA", "SummarizedExperiment",
            function(x,...){
                .Deprecated(msg = paste0("'calculateCCA' is deprecated. ",
                                        "Use 'getCCA' instead."))
                getCCA(x,...)
            }
)
