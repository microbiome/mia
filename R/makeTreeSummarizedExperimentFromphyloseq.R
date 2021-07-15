#' Coerce a \code{phyloseq} object to a \code{TreeSummarizedExperiment}
#'
#' \code{makeTreeSummarizedExperimentFromphyloseq} converts \code{phyloseq}
#' objects into \code{TreeSummarizedExperiment} objects.
#'
#' All data stored in a \code{phyloseq} object is transfered.
#'
#' @param obj a \code{phyloseq} object
#'
#' @return An object of class \code{TreeSummarizedExperiment}
#'
#' @importFrom S4Vectors SimpleList DataFrame
#' @importFrom SummarizedExperiment colData colData<-
#'
#' @export
#'
#' @name makeTreeSummarizedExperimentFromphyloseq
#' @seealso
#' \code{\link[=makeSummarizedExperimentFromBiom]{makeSummarizedExperimentFromBiom}}
#' \code{\link[=makeTreeSummarizedExperimentFromDADA2]{makeTreeSummarizedExperimentFromDADA2}}
#' \code{\link[=loadFromQIIME2]{loadFromQIIME2}}
#' \code{\link[=loadFromMothur]{loadFromMothur}}
#'
#' @examples
#' if (requireNamespace("phyloseq")) {
#'     data(GlobalPatterns, package="phyloseq")
#'     makeTreeSummarizedExperimentFromphyloseq(GlobalPatterns)
#'     data(enterotype, package="phyloseq")
#'     makeTreeSummarizedExperimentFromphyloseq(enterotype)
#'     data(esophagus, package="phyloseq")
#'     makeTreeSummarizedExperimentFromphyloseq(esophagus)
#' }
makeTreeSummarizedExperimentFromphyloseq <- function(obj) {
    # input check
    .require_package("phyloseq")
    if(!is(obj,"phyloseq")){
        stop("'obj' must be a 'phyloseq' object")
    }
    #
    assays <- SimpleList(counts = obj@otu_table@.Data)
    rowData <- S4Vectors:::make_zero_col_DataFrame(nrow(assays$counts))
    colData <- S4Vectors:::make_zero_col_DataFrame(ncol(assays$counts))
    if(!is.null(obj@tax_table@.Data)){
        rowData <- DataFrame(data.frame(obj@tax_table@.Data))
    }
    if(!is.null(obj@sam_data)){
        colData <- DataFrame(data.frame(obj@sam_data))
    }
    if(!is.null(obj@phy_tree)){
        rowTree <- obj@phy_tree
    } else {
        rowTree <- NULL
    }
    if (!is.null(obj@refseq)) {
        referenceSeq <- obj@refseq
    } else {
        referenceSeq <- NULL
    }
    TreeSummarizedExperiment(assays = assays,
                             rowData = obj@tax_table@.Data,
                             colData = colData,
                             rowTree = rowTree,
                             referenceSeq = referenceSeq)
}
