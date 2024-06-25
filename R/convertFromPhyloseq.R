#' @param phy a \code{phyloseq} object
#' 
#' @details 
#' \code{convertFromPhyloseq} converts \code{phyloseq}
#' objects into 
#' \code{\link[TreeSummarizedExperiment:TreeSummarizedExperiment-class]{TreeSummarizedExperiment}} objects.
#' All data stored in a \code{phyloseq} object is transferred.
#'
#' @return 
#' \code{convertFromPhyloseq} returns an object of class 
#' \code{\link[TreeSummarizedExperiment:TreeSummarizedExperiment-class]{TreeSummarizedExperiment}}
#'
#' @importFrom S4Vectors SimpleList DataFrame make_zero_col_DFrame
#' @importFrom SummarizedExperiment colData colData<-
#'
#' @export
#'
#' @rdname convert
#' 
#' @seealso
#' \code{\link[=convert]{convertFromBIOM}}
#' \code{\link[=convert]{convertFromDADA2}}
#' \code{\link[=importQIIME2]{importQIIME2}}
#' \code{\link[=importMothur]{importMothur}}
#'
#' @examples
#' 
#' ### Coerce a phyloseq object to a TreeSE object
#' if (requireNamespace("phyloseq")) {
#'     data(GlobalPatterns, package="phyloseq")
#'     convertFromPhyloseq(GlobalPatterns)
#'     data(enterotype, package="phyloseq")
#'     convertFromPhyloseq(enterotype)
#'     data(esophagus, package="phyloseq")
#'     convertFromPhyloseq(esophagus)
#' }
convertFromPhyloseq <- function(x) {
    # input check
    .require_package("phyloseq")
    if(!is(x,"phyloseq")){
        stop("'x' must be a 'phyloseq' object")
    }
    #
    # Get the assay
    counts <- x@otu_table@.Data
    # Check the orientation, and transpose if necessary
    if( !x@otu_table@taxa_are_rows ){
        counts <- t(counts)
    }
    # Create a list of assays
    assays <- SimpleList(counts = counts)
    
    if(!is.null(x@tax_table@.Data)){
        rowData <- DataFrame(data.frame(x@tax_table@.Data))
    } else{
        rowData <- S4Vectors::make_zero_col_DFrame(nrow(assays$counts))
        rownames(rowData) <- rownames(assays$counts)
    }
    if(!is.null(x@sam_data)){
        colData <- DataFrame(data.frame(x@sam_data))
    } else{
        colData <- S4Vectors::make_zero_col_DFrame(ncol(assays$counts))
        rownames(colData) <- colnames(assays$counts)
    }
    if(!is.null(x@phy_tree)){
        rowTree <- x@phy_tree
    } else {
        rowTree <- NULL
    }
    if (!is.null(x@refseq)) {
        referenceSeq <- x@refseq
    } else {
        referenceSeq <- NULL
    }
    TreeSummarizedExperiment(assays = assays,
                            rowData = rowData,
                            colData = colData,
                            rowTree = rowTree,
                            referenceSeq = referenceSeq)
}
