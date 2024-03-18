#' @export
setGeneric("convert", signature = c("x"),
           function(x,...)
               standardGeneric("convert"))

#' @export
setMethod("convert", signature = c(x = "SummarizedExperiment"),
          function(x,...){
              .makePhyloseqFromSE(x,...)
          }
)

#' @export
setMethod("convert", signature = c(x = "TreeSummarizedExperiment"),
          function(x,...){
              .makePhyloseqFromTreeSE(x,...)
          }
)

#' @export
setMethod("convert", signature = c(x = "dada"),
          function(x,...){
              .makeTreeSEFromDADA2(x,...)
          }
)

#' @export
setMethod("convert", signature = c(x = "phyloseq"),
          function(x,...){
              .makeTreeSEFromPhyloseq(x,...)
          }
)

#' @export
setMethod("convert", signature = c(x = "biom"),
          function(x,...){
              .makeTreeSEFromBiom(x,...)
          }
)

##################### makeTreeSummarizedExperimentFromPhyloseq #################

.makeTreeSEFromPhyloseq <- function(obj) {
    # input check
    .require_package("phyloseq")
    if(!is(obj,"phyloseq")){
        stop("'obj' must be a 'phyloseq' object")
    }
    #
    # Get the assay
    counts <- obj@otu_table@.Data
    # Check the orientation, and transpose if necessary
    if( !obj@otu_table@taxa_are_rows ){
        counts <- t(counts)
    }
    # Create a list of assays
    assays <- SimpleList(counts = counts)
    
    if(!is.null(obj@tax_table@.Data)){
        rowData <- DataFrame(data.frame(obj@tax_table@.Data))
    } else{
        rowData <- S4Vectors::make_zero_col_DFrame(nrow(assays$counts))
        rownames(rowData) <- rownames(assays$counts)
    }
    if(!is.null(obj@sam_data)){
        colData <- DataFrame(data.frame(obj@sam_data))
    } else{
        colData <- S4Vectors::make_zero_col_DFrame(ncol(assays$counts))
        rownames(colData) <- colnames(assays$counts)
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
                             rowData = rowData,
                             colData = colData,
                             rowTree = rowTree,
                             referenceSeq = referenceSeq)
}

################### makeTreeSummarizedExperimentFromDADA2 ######################

.makeTreeSEFromDADA2 <- function(x,...) {
    # input checks
    .require_package("dada2")
    .require_package("stringr")
    #
    mergers <- dada2::mergePairs(x,...)
    seqtab <- dada2::makeSequenceTable(mergers)
    seqtab <- t(seqtab)
    # generate row and col names
    rName <- paste0("ASV",
                    stringr::str_pad(seq.int(1L,nrow(seqtab)),
                                     nchar(nrow(seqtab)) + 1L,
                                     pad="0"))
    cName <- colnames(seqtab)
    # retrieve count data and reference sequence
    assays <- S4Vectors::SimpleList(counts = unname(seqtab))
    refseq <- Biostrings::DNAStringSet(rownames(seqtab))
    # construct ME an name rows and cols
    output <- TreeSummarizedExperiment(assays = assays,
                                       referenceSeq = refseq)
    colnames(output) <- cName
    rownames(output) <- rName
    output
}

################### makeTreeSummarizedExperimentFromBiom #######################

loadFromBiom <- function(file, ...) {
    .require_package("biomformat")
    biom <- biomformat::read_biom(file)
    makeTreeSEFromBiom(biom, ...)
}

.makeTreeSEFromBiom <- function(
        obj, removeTaxaPrefixes = FALSE, rankFromPrefix = FALSE,
        remove.artifacts = FALSE, ...){
    # input check
    .require_package("biomformat")
    if(!is(obj,"biom")){
        stop("'obj' must be a 'biom' object", call. = FALSE)
    }
    if( !.is_a_bool(removeTaxaPrefixes) ){
        stop("'removeTaxaPrefixes' must be TRUE or FALSE.", call. = FALSE)
    }
    if( !.is_a_bool(rankFromPrefix) ){
        stop("'rankFromPrefix' must be TRUE or FALSE.", call. = FALSE)
    }
    if( !.is_a_bool(remove.artifacts) ){
        stop("'remove.artifacts' must be TRUE or FALSE.", call. = FALSE)
    }
    #
    counts <- as(biomformat::biom_data(obj), "matrix")
    sample_data <- biomformat::sample_metadata(obj)
    feature_data <- biomformat::observation_metadata(obj)
    
    # colData is initialized with empty tables with rownames if it is NULL
    if( is.null(sample_data) ){
        sample_data <- S4Vectors::make_zero_col_DFrame(ncol(counts))
        rownames(sample_data) <- colnames(counts)
        # Otherwise convert it into correct format if it is a list
    } else if( is(sample_data, "list") ){
        # Merge list of data.frames into one
        sample_data <- bind_rows(sample_data)
        sample_data < as.data.frame(sample_data)
    }
    # rowData is initialized with empty tables with rownames if it is NULL
    if( is.null(feature_data) ){
        feature_data <- S4Vectors::make_zero_col_DFrame(nrow(counts))
        rownames(feature_data) <- rownames(counts)
        # Otherwise convert it into correct format if it is a list
    } else if( is(feature_data, "list") ){
        # Feature data is a list of taxa info. Dfs are merged together 
        # differently than sample metadata since the column names are only 
        # "Taxonomy". If there is only one taxonomy level, the column name does 
        # not get a suffix.
        # --> bind rows based on the index of column.
        
        # Get the maximum length of list
        max_length <- max( lengths(feature_data) )
        # Get the column names from the taxa info that has all the levels that 
        # occurs in the data
        colnames <- names( head(
            feature_data[ lengths(feature_data) == max_length ], 1)[[1]])
        # Convert the list so that all individual taxa info have the max length
        # of the list objects. All vectors are appended with NAs, if they do not
        # have all the levels. E.g., if only Kingdom level is found, all lower
        # ranks are now NA
        feature_data <- lapply(feature_data, function(x){
            length(x) <- max_length 
            return(x)
        })
        # Create a data.frame from the list
        feature_data <- do.call(rbind, feature_data)
        # Transposing feature_data and make it df object
        feature_data <- as.data.frame(feature_data)
        # Add column that includes all the data
        feature_data[["taxonomy_unparsed"]] <- apply(feature_data, 1, paste0, 
                                                    collapse = ";")
        # Add correct colnames
        colnames(feature_data) <- c(colnames, "taxonomy_unparsed")
    }
    # If there is only one column in the feature data, it is the most probable
    # that the taxonomy is not parsed. Try to parse it.
    if( ncol(feature_data) == 1 ){
        colnames(feature_data) <- "taxonomy_unparsed"
        tax_tab <- .parse_taxonomy(feature_data, 
                                    column_name = colnames(feature_data))
        feature_data <- cbind(tax_tab, feature_data)
        feature_data <- as.data.frame(feature_data)
    }
    
    # Clean feature_data from possible character artifacts if specified
    if( remove.artifacts ){
        feature_data <- .detect_taxa_artifacts_and_clean(feature_data, ...)
    }
    
    # Replace taxonomy ranks with ranks found based on prefixes
    if( rankFromPrefix && all(
        unlist(lapply(colnames(feature_data),
                      function(x) !x %in% TAXONOMY_RANKS)))){
        # Find ranks
        ranks <- lapply(colnames(feature_data),
                        .replace_colnames_based_on_prefix, x=feature_data)
        # Replace old ranks with found ranks
        colnames(feature_data) <- unlist(ranks)
    }
    
    # Remove prefixes if specified and rowData includes info
    if(removeTaxaPrefixes && ncol(feature_data) > 0){
        feature_data <- .remove_prefixes_from_taxa(feature_data, ...)
    }
    
    # Adjust row and colnames
    rownames(counts) <- rownames(feature_data) <- biomformat::rownames(obj)
    colnames(counts) <- rownames(sample_data) <- biomformat::colnames(obj)
    
    # Convert into DataFrame
    sample_data <- DataFrame(sample_data)
    feature_data <- DataFrame(feature_data)
    # Convert into list
    assays <- SimpleList(counts = counts)
    
    # Create TreeSE
    tse <- TreeSummarizedExperiment(
        assays = assays,
        colData = sample_data,
        rowData = feature_data)
    return(tse)
}

################### makePhyloseqFromSummarizedExperiment #######################

.makePhyloseqFromSE <- function(x, assay.type = "counts", assay_name = NULL,
                                ...){
    # Input check
    .require_package("phyloseq")
    # Check that tse do not have zero rows
    if(!all(dim(x) > 0)){
        stop("'x' contains zero rows. 'x' can not be converted
                 to a phyloseq object.",
             call. = FALSE)
    }
    
    if (!is.null(assay_name)) {
        .Deprecated(old="assay_name", new="assay.type", "Now assay_name is 
                    deprecated. Use assay.type instead.")
    }
    
    # Check assay.type
    .check_assay_present(assay.type, x)
    
    # phyloseq object requires nonduplicated rownames. If there are 
    # duplicated rownames, they are converted so that they are unique
    if( any(duplicated(rownames(x))) ){
        rownames(x) <- getTaxonomyLabels(x)
    }
    # List of arguments
    args = list()
    # Gets the abundance data from assay, and converts it to otu_table
    otu_table <- as.matrix(assay(x, assay.type))
    otu_table <- phyloseq::otu_table(otu_table, taxa_are_rows = TRUE)
    # Adds to the list
    args[["otu_table"]] <- otu_table
    
    # If rowData includes information
    if(!( length(rowData(x)[,taxonomyRanks(x)]) == 0 ||
          is.null((rowData(x)[,taxonomyRanks(x)])) )){
        # Converts taxonomy table to characters if it's not already
        rowData(x) <- DataFrame(lapply(rowData(x), as.character))
        # Gets the taxonomic data from rowData, and converts it to tax_table
        tax_table <- as.matrix(rowData(x)[,taxonomyRanks(x),drop=FALSE])
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

################# makePhyloseqFromTreeSummarizedExperiment #####################

.makePhyloseqFromTreeSE <- function(x, tree_name = "phylo", ...){
    # If rowTrees exist, check tree_name
    if( length(x@rowTree) > 0 ){
        .check_rowTree_present(tree_name, x)
        # Subset the data based on the tree
        x <- x[ rowLinks(x)$whichTree == tree_name, ]
        add_phy_tree <- TRUE
    } else{
        add_phy_tree <- FALSE
    }
    #
    
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
    
    # Add phylogenetic tree
    if( add_phy_tree ){
        phy_tree <- .get_rowTree_for_phyloseq(x, tree_name)
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
    if( !is.null(referenceSeq(x)) ){
        # Get referenceSeqs
        refseq <- .get_referenceSeq_for_phyloseq(x, ...)
        # IF refSeq passed the test, add it
        if( !is.null(refseq) ){
            # Convert it to phyloseq object
            refseq <- phyloseq::refseq(refseq)
            # If the object is a phyloseq object, adds refseq to it
            if(is(obj,"phyloseq")){
                obj <- phyloseq::merge_phyloseq(obj, refseq)
            } else{
                # Adds to the list
                args[["refseq"]] <- refseq
            }
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

################################ HELP FUNCTIONS ################################

# This function removes prefixes from taxonomy names
.remove_prefixes_from_taxa <- function(
        feature_tab, prefixes = "sk__|([dkpcofgs]+)__",
        only.taxa.col = TRUE, ...){
    if( !.is_a_bool(only.taxa.col) ){
        stop("'only.taxa.col' must be TRUE or FALSE.", call. = FALSE)
    }
    #
    # Subset by taking only taxonomy info if user want to remove the pattern 
    # only from those. (Might be too restricting, e.g., if taxonomy columns are 
    # not detected in previous steps. That is way the default is FALSE)
    if( only.taxa.col ){
        ind <- tolower(colnames(feature_tab)) %in% TAXONOMY_RANKS
        temp <- feature_tab[, ind, drop = FALSE]
    } else{
        ind <- rep(TRUE, ncol(feature_tab))
        temp <- feature_tab
    }
    
    # If there are columns left for removing the pattern
    if( ncol(temp) > 0 ){
        # Remove patterns
        temp <- lapply(
            temp, gsub, pattern = prefixes, replacement = "")
        temp <- as.data.frame(temp)
        # If cell had only prefix, it is now empty string. Convert to NA
        temp[ temp == "" ] <- NA
        # Combine table
        feature_tab[, ind] <- temp
    }
    return(feature_tab)
}

# Find taxonomy rank based on prefixes. If found, return
# corresponding rank. Otherwise, return the original
# rank that is fed to function.
.replace_colnames_based_on_prefix <- function(colname, x){
    # Get column
    col = x[ , colname]
    # List prefixes
    prefixes <- c(
        "^d__",
        "^k__",
        "^p__",
        "^c__",
        "^o__",
        "^f__",
        "^g__",
        "^s__"
    )
    # Find which prefix is found from each column value, if none.
    found_rank <- lapply(prefixes, FUN = function(pref){
        all(grepl(pattern = pref, col) | is.na(col)) && !all(is.na(col))
    })
    found_rank <- unlist(found_rank)
    # If only one prefix was found (like it should be), get the corresponding
    # rank name.
    if( sum(found_rank) == 1 ){
        colname <- TAXONOMY_RANKS[found_rank]
        # Make it capitalized
        colname <- paste0(toupper(substr(colname, 1, 1)),
                          substr(colname, 2, nchar(colname)))
    }
    return(colname)    
}

# Detect and clean non wanted characters from Taxonomy data if needed.
.detect_taxa_artifacts_and_clean <- function(x, pattern = "auto", ...) {
    #
    if( !.is_non_empty_character(pattern) ){
        stop("'pattern' must be a single character value.", call. = FALSE)
    }
    #
    row_names <- rownames(x)
    # Remove artifacts
    if( pattern == "auto" ){
        .require_package("stringr")
        # Remove all but these characters
        pattern <- "[[:alnum:]]|-|_|\\[|\\]|,|;\\||[[:space:]]"
        x <- lapply(x, function(col){
            # Take all specified characters as a matrix where each column is a 
            # character
            temp <- stringr::str_extract_all(col, pattern = pattern, 
                                            simplify = TRUE)
            # Collapse matrix to strings
            temp <- apply(temp, 1, paste, collapse = "")
            # Now NAs are converted into characters. Convert them back
            temp[ temp == "NA" ] <- NA
            # Convert also empty strings to NA
            temp[ temp == "" ] <- NA
            return(temp)
        })
    } else{
        # Remove pattern specified by user
        x <- lapply(x, gsub, pattern = pattern, replacement = "")
    }
    x <- as.data.frame(x)
    # Add rownames because they are dropped while removing artifacts
    rownames(x) <- row_names
    return(x)
}

# If tips do not match with rownames, prune the tree
.get_x_with_pruned_tree <- function(x, tree_name){
    # Get rowLinks
    row_links <- rowLinks(x)
    # Gets node labels
    node_labs <- row_links[ , "nodeLab"]
    # Prunes the tree
    tree_pruned <- ape::keep.tip(rowTree(x), node_labs)
    # Replace tip labels with corresponding rownames
    tree_pruned$tip.label <- rownames(x)
    # Assigns the pruned tree back to TSE object
    rowTree(x) <- tree_pruned
    warning("Tips of rowTree are renamed to match rownames.", call. = FALSE)
    return(x)
}

# In phyloseq, tips and rownames must match
.get_rowTree_for_phyloseq <- function(x, tree_name){
    # Check if the rowTree's tips match with rownames:
    # tips labels are found from rownames
    if( any(!( rowTree(x, tree_name)$tip.label) %in% rownames(x)) ){
        # If rowtree do not match, tree is pruned
        x <- .get_x_with_pruned_tree(x, tree_name)
    }
    # Get rowTree
    phy_tree <- rowTree(x, tree_name)
    # Convert rowTree to phyloseq object
    phy_tree <- phyloseq::phy_tree(phy_tree)
    
    return(phy_tree)
}

.get_referenceSeq_for_phyloseq <- function(x, referenceSeq = 1, ...){
    # Get reference seqs
    refSeqs <- referenceSeq(x)
    # Is referenceSeq a list / does it contain multiple DNA sets
    is_list <- is(refSeqs, "DNAStringSetList")
    
    # Take only one set, if it is a list
    if( is_list ){
        # Check referenceSeq
        if( !( (.is_non_empty_string(referenceSeq) && 
                referenceSeq %in% names(refSeqs)) ||
               (.is_an_integer(referenceSeq) &&
                (referenceSeq>0 && referenceSeq<=length(refSeqs))) )
        ){
            stop("'referenceSeq' must be a non-empty single character value or 
                 an integer specifying the DNAStringSet from DNAStringSetList.",
                 call. = FALSE)
        }
        # Get specified referenceSeq
        refSeqs <- refSeqs[[referenceSeq]]
        warning("Use 'referenceSeq' to specify DNA set from DNAStringSetList. ",
                "Current choice is '", referenceSeq, "'.", 
                call. = FALSE)
    }
    # Check if all rownames have referenceSeqs
    if( !(all(rownames(x) %in% names(refSeqs)) &&
          all(names(refSeqs) %in% rownames(x) )) ){
        warning("referenceSeq does not match with rownames so they are 
                discarded.",
                call. = FALSE)
        refSeqs <- NULL
    }
    return(refSeqs)
}
