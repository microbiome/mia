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
loadFromBiom <- function(file, ...) {
    .Deprecated(msg = paste0("'loadFromBiom' is deprecated.",
                            " Use 'importBIOM' instead."))
    importBiom(file,...)
}

#' @rdname deprecate
#' @export
loadFromQIIME2 <- function(featureTableFile,
                         taxonomyTableFile = NULL,
                         sampleMetaFile = NULL,
                         featureNamesAsRefSeq = TRUE,
                         refSeqFile = NULL,
                         phyTreeFile = NULL,
                         ...) {
    .Deprecated(msg = paste0("'loadFromQIIME2' is deprecated.",
                            " Use 'importQIIME2' instead."))
    importQIIME2(featureTableFile,
                 taxonomyTableFile = NULL,
                 sampleMetaFile = NULL,
                 featureNamesAsRefSeq = TRUE,
                 refSeqFile = NULL,
                 phyTreeFile = NULL,
                 ...)
}

#' @rdname deprecate
#' @export
loadFromQZA <- function(file, temp = tempdir(), ...) {
    .Deprecated(msg = paste0("'loadFromQZA' is deprecated.",
                            " Use 'importQZA' instead."))
    importQZA(file, temp = tempdir(), ...)
}

#' @rdname deprecate
#' @export
loadFromMothur <- function(sharedFile, taxonomyFile = NULL, designFile = NULL) {
    .Deprecated(msg = paste0("'loadFromMothur' is deprecated.",
                             " Use 'importMothur' instead."))
    importMothur(sharedFile, taxonomyFile = NULL, designFile = NULL)
}

#' @rdname deprecate
#' @export   
loadFromMetaphlan <- function(
        file, colData = sample_meta, sample_meta = NULL, phy_tree = NULL, ...) {
    .Deprecated(msg = paste0("'loadFromMetaphlan' is deprecated.",
                             " Use 'importMetaPhlAn' instead."))
    importMetaPhlAn(file, colData = sample_meta, sample_meta = NULL, 
                    phy_tree = NULL, ...)
}

#' @rdname deprecate
#' @export    
loadFromHumann <- function(file, colData = NULL, ...) {
    .Deprecated(msg = paste0("'loadFromHumann' is deprecated.",
                             " Use 'importHUMAnN' instead."))
    importHUMAnN(file, colData = NULL, ...)
}
