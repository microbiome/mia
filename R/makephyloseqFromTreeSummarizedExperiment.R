#' Create a phyloseq object from a TreeSummarizedExperiment object
#'
#' @param x a \code{TreeSummarizedExperiment} object
#'
#' @param abund_values
#' A single character value for selecting the
#'   \code{\link[SummarizedExperiment:SummarizedExperiment-class]{assay}}
#'   to be included in the phyloseq object that is created.
#'
#'  @details
#' \code{makePhyloseqFromTreeSummarizedExperiment} is used for creating a phyloseq
#' object from TreeSummarizedExperiment object.
#'
#' @return
#' An object of class \code{Phyloseq} object.
#'
#' @name makePhyloseqFromTreeSummarizedExperiment
#' @export
#'
#' @author Leo Lahti and Tuomas Borman. Contact: \url{microbiome.github.io}
#'
#' @examples
#' # Get TSE object
#' data("GlobalPatterns")
#' TSE <- GlobalPatterns
#'
#' # Create a phyloseq object from it
#' phy <- makePhyloseqFromTreeSummarizedExperiment(TSE)
#'
#' # By default the chosen table is counts, but if there is other tables,
#' # they can be chosen with abund_values.
#'
#' # Counts relative abundances table
#' TSE <- transformCounts(TSE, method = "relabundance")
#' phy2 <- makePhyloseqFromTreeSummarizedExperiment(TSE, abund_values = "relabundance")
#'
NULL

#' @rdname makePhyloseqFromTreeSummarizedExperiment
#' @export
setGeneric("makePhyloseqFromTreeSummarizedExperiment", signature = c("x"),
           function(x, abund_values = "counts")
               standardGeneric("makePhyloseqFromTreeSummarizedExperiment"))


#' @rdname makePhyloseqFromTreeSummarizedExperiment
#' @export
setMethod("makePhyloseqFromTreeSummarizedExperiment", signature = c(x = "SummarizedExperiment"),
          function(x, abund_values = "counts"){

            # Input check
            # Check object
            if(!is(x,"TreeSummarizedExperiment")){
                stop("'x' must be a 'TreeSummarizedExperiment' object")
            }

            # Check abund_values
            .check_abund_values(abund_values, x)

            # Get the assay from x
            otu_table <- assay(x, abund_values)
            # Convert to otu_table
            otu_table <- otu_table(otu_table, taxa_are_rows = TRUE)

            # If rowData has information, stores it to tax_table and converts is to phyloseq's tax_table.
            # Otherwise, tax_table is NULL, and it is converted to phyloseq's tax_table.
            if(nrow(rowData(x))>0 && ncol(rowData(x))>0){
                tax_table <- as.matrix(rowData(x))
                tax_table <- tax_table(tax_table)
            } else {
                tax_table <- NULL
                tax_table <- tax_table(tax_table, errorIfNULL = FALSE)
            }

            # If colData has information, stores it to sample_data and converts is to phyloseq's sample_data.
            # Otherwise, sample_data is NULL, and it is converted to phyloseq's sample_data.
            if(nrow(colData(x))>0 && ncol(colData(x))>0){
                sample_data <- data.frame(colData(x))
                sample_data <- sample_data(sample_data)
            } else {
                sample_data <- NULL
                sample_data <- sample_data(sample_data, errorIfNULL = FALSE)
            }

            # If rowTree has information, stores it to phy_tree and converts is to phyloseq's phy_tree.
            # Otherwise, phy_tree is NULL, and it is converted to phyloseq's phy_tree.
            if(!is.null(rowTree(x))){
                phy_tree <- rowTree(x)
                phy_tree <- phy_tree(phy_tree)
            } else {
                phy_tree <- NULL
                phy_tree <- phy_tree(phy_tree, errorIfNULL = FALSE)
            }

            # If referenceSeq has information, stores it to refseq and converts is to phyloseq's refseq.
            # Otherwise, refseq is NULL, and it is converted to phyloseq's refseq.
            if(!is.null(referenceSeq(x))){
                refseq <- referenceSeq(x)
                refseq <- refseq(refseq)
            } else {
                refseq <- NULL
                refseq <- refseq(refseq, errorIfNULL = FALSE)
            }

            # Creates a phyloseq object
            phyloseq <- phyloseq(otu_table,
                                 tax_table,
                                 sample_data,
                                 phy_tree,
                                 refseq)

            return(phyloseq)
})
