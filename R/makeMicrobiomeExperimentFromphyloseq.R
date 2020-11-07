#' Coerce phyloseq object
#'
#' @param obj a \code{phyloseq} object
#'
#' @return An object of class MicrobiomeExperiment
#'
#' @importFrom S4Vectors SimpleList DataFrame
#' @importFrom SummarizedExperiment colData colData<-
#'
#' @export
#'
#' @examples
#' if (requireNamespace("phyloseq")) {
#'     data(GlobalPatterns, package="phyloseq")
#'     makeMicrobiomeExperimentFromphyloseq(GlobalPatterns)
#'     data(enterotype, package="phyloseq")
#'     makeMicrobiomeExperimentFromphyloseq(enterotype)
#'     data(esophagus, package="phyloseq")
#'     makeMicrobiomeExperimentFromphyloseq(esophagus)
#'
#'     # MovingPictures data requires minimal manual fix
#'     data(MovingPictures, package = "MicrobeDS")
#'     MovingPictures@phy_tree$node.label[[50]] <- ""
#'     makeMicrobiomeExperimentFromphyloseq(MovingPictures)
#'
#' }
makeMicrobiomeExperimentFromphyloseq <- function(obj) {
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
    MicrobiomeExperiment(assays = assays,
                         rowData = obj@tax_table@.Data,
                         colData = colData,
                         rowTree = rowTree,
                         referenceSeq = referenceSeq)
}
