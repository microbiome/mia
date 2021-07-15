#' Create a phyloseq object from a TreeSummarizedExperiment object
#'
#' This function creates a phyloseq object from a TreeSummarizedExperiment
#' object. By using \code{abund_values}, it is possible to specify which table
#' from \code{assay} is added to the phyloseq object.
#'
#' @param x a \code{TreeSummarizedExperiment} object
#'
#' @param abund_values A single character value for selecting the
#'   \code{\link[SummarizedExperiment:SummarizedExperiment-class]{assay}} to be
#'   included in the phyloseq object that is created. By default, it is counts
#'   table.
#'
#' @param ... additional arguments
#'
#' @details
#' \code{makePhyloseqFromTreeSummarizedExperiment} is used for creating a
#' phyloseq object from TreeSummarizedExperiment object.
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
#' # Get tse object
#' data(GlobalPatterns)
#' tse <- GlobalPatterns
#'
#' # Create a phyloseq object from it
#' phy <- makePhyloseqFromTreeSummarizedExperiment(tse)
#' phy
#'
#' # By default the chosen table is counts, but if there are other tables,
#' # they can be chosen with abund_values.
#'
#' # Counts relative abundances table
#' tse <- transformCounts(tse, method = "relabundance")
#' phy2 <- makePhyloseqFromTreeSummarizedExperiment(tse, abund_values = "relabundance")
#' phy2
NULL

#' @rdname makePhyloseqFromTreeSummarizedExperiment
#' @export
setGeneric("makePhyloseqFromTreeSummarizedExperiment", signature = c("x"),
           function(x, ...)
               standardGeneric("makePhyloseqFromTreeSummarizedExperiment"))


#' @rdname makePhyloseqFromTreeSummarizedExperiment
#' @export
setMethod("makePhyloseqFromTreeSummarizedExperiment",
          signature = c(x = "SummarizedExperiment"),
    function(x, abund_values = "counts"){
        # Input check
        .require_package("phyloseq")
        # Check that tse do not have zero rows
        if(!all(dim(x) > 0)){
            stop("'x' contains zero rows. 'x' can not be converted
                 to a phyloseq object.",
                 call. = FALSE)
        }
        # Check abund_values
        .check_assay_present(abund_values, x)
        
        # phyloseq object requires nonduplicated rownames. If there are 
        # duplicated rownames, they are converted so that they are unique
        if( any(duplicated(rownames(x))) ){
            rownames(x) <- getTaxonomyLabels(x)
        }
        # List of arguments
        args = list()
        # Gets the abundance data from assay, and converts it to otu_table
        otu_table <- as.matrix(assay(x, abund_values))
        otu_table <- phyloseq::otu_table(otu_table, taxa_are_rows = TRUE)
        # Adds to the list
        args[["otu_table"]] <- otu_table

        # If rowData includes information
        if(!( length(rowData(x)[,taxonomyRanks(x)]) == 0 ||
              is.null((rowData(x)[,taxonomyRanks(x)])) )){
            # Gets the taxonomic data from rowData, and converts it to tax_table
            tax_table <- as.matrix(rowData(x)[,taxonomyRanks(x)])
            tax_table <- phyloseq::tax_table(tax_table)
            # Adds to the list
            args[["tax_table"]] <- tax_table
        }
        
        # If colData includes information
        if(!( length(colData(x)) == 0 || is.null(ncol(colData(x))) )){
            # Gets the feature_data from colData and converts it to sample_data
            sample_data <- as.data.frame(colData(x))
            sample_data <- phyloseq::sample_data(sample_data)
            # Adds to the list
            args[["sample_data"]] <- sample_data
        }
        
        # Creates a phyloseq object
        phyloseq <- do.call(phyloseq::phyloseq, args)
        return(phyloseq)
    }
)

#' @rdname makePhyloseqFromTreeSummarizedExperiment
#' @export
setMethod("makePhyloseqFromTreeSummarizedExperiment",
          signature = c(x = "TreeSummarizedExperiment"),
    function(x, ...){
        # phyloseq and tree objects require nonduplicated rownames. If there are 
        # duplicated rownames, they are converted so that they are unique
        if( any(duplicated(rownames(x))) ){
            rownames(x) <- getTaxonomyLabels(x)
        }
        # Gets otu_table object from the function above this, if tse contains
        # only abundance table.
        # Otherwise, gets a phyloseq object with otu_table, and tax_table
        # and/or sample_data
        obj <- callNextMethod()
        # List of arguments
        args = list()
        # Adds to the list of arguments, if 'obj' is not a phyloseq object
        # i.e. is an otu_table
        if(!is(obj,"phyloseq")){
            # Adds otu_table to the list
            args[["otu_table"]] <- obj
        }

        # If rowTree has information, stores it to phy_tree and converts is to
        # phyloseq's phy_tree.
        if(!( length(rowTree(x)) == 0 || is.null(rowTree(x)) )){
            # If rowTree exists, checks if the rowTree match with rownames:
            # tips labels are found from rownames
            if( any(!( rowTree(x)$tip) %in% rownames(x)) ){
                # If rowtree do not match, tree is pruned
                x <- .get_x_with_pruned_tree(x)
            }
            phy_tree <- rowTree(x)
            phy_tree <- phyloseq::phy_tree(phy_tree)
            # If the object is a phyloseq object, adds phy_tree to it
            if(is(obj,"phyloseq")){
                phyloseq::phy_tree(obj) <- phy_tree
            } else{
                # Adds to the list
                args[["phy_tree"]] <- phy_tree
            }
        }

        # If referenceSeq has information, stores it to refseq and converts is
        # to phyloseq's refseq.
        if(!( length(referenceSeq(x)) == 0 || is.null(referenceSeq(x)) )){
            refseq <- referenceSeq(x)
            refseq <- phyloseq::refseq(refseq)
            # If the object is a phyloseq object, adds refseq to it
            if(is(obj,"phyloseq")){
                obj <- phyloseq::merge_phyloseq(obj, refseq)
            } else{
                # Adds to the list
                args[["refseq"]] <- refseq
            }
        }
        
        # If 'obj' is not a phyloseq object, creates one.
        if(!is(obj,"phyloseq")){
            # Creates a phyloseq object
            phyloseq <- do.call(phyloseq::phyloseq, args)
        } else{
            phyloseq <- obj
        }
        phyloseq
    }
)

################################ HELP FUNCTIONS ################################

.get_x_with_pruned_tree <- function(x){
    # Gets node labels
    node_labs <- rowLinks(x)$nodeLab
    # Gets the corresponding rownames
    node_labs_rownames <- rownames(rowLinks(x))
    # Prunes the tree
    tree_pruned <- ape::keep.tip(rowTree(x), node_labs)
    # Replace tip labels with corresponding rownames
    tree_pruned$tip.label <- node_labs_rownames
    # Assigns the pruned tree back to TSE object
    rowTree(x) <- tree_pruned
    warning("rowTree is pruned to match rownames.")
    return(x)
}
