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
loadFromQZA <- function(...) {
    .Deprecated(msg = paste0("'loadFromQZA' is deprecated.",
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
