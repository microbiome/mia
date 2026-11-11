# This function can be used to unify the group id vector. It can be any
# kind of vector, but this converts it to factor.
.norm_f <- function(
        i, f, dim.type = c("rows","columns"), empty.rm = FALSE, ...){
    if(!.is_a_bool(empty.rm)){
        stop("'empty.rm' must be TRUE or FALSE.", call. = FALSE)
    }
    dim.type <- match.arg(dim.type)
    if(!is.character(f) && !is.factor(f)){
        stop("'f' must be a factor or character vector coercible to a ",
            "meaningful factor.", call. = FALSE)
    }
    if(i != length(f)){
        stop("'f' must have the same number of ",dim.type," as 'x'",
            call. = FALSE)
    }
    # This is done otherwise we lose NA values
    if( !empty.rm && any(is.na(f)) ){
        f <- as.character(f)
        f[ is.na(f) ] <- "NA"
    }
    if(is.character(f)){
        f <- factor(f)
    }
    f
}

# If the group is specified as a column name of rowData or colData, this
# function returns the values of that column. Otherwise, the group is returned
# as it is.
.get_group_vector <- function(x, f, by){
    rowData_FUN <- switch(by, rowData, colData)
    nrow_FUN <- switch(by, nrow, ncol)
    if( .is_a_string(f) && f %in% colnames(rowData_FUN(x)) ){
        f <- rowData_FUN(x)[[ f ]]
    } else if( .is_a_string(f) && !f %in% colnames(rowData_FUN(x)) &&
                nrow_FUN(x) != 1L){
        stop("Specified grouping variable ('group' or 'rank' argument) must ",
            "specify a column from ", switch(by, "row", "col"), "Data(x).",
            call. = FALSE)
    }
    return(f)
}

# When we merge rows or columns, first member of group is kept as default
# (in colData or rowData). This function controls this and allows user to
# specify some other element than the first one.
.norm_archetype <- function(f, archetype){
    if(length(archetype) > 1L){
        if(length(levels(f)) != length(archetype)){
            stop("length of 'archetype' must have the same length as ",
                "levels('f')",
                call. = FALSE)
        }
    }
    f_table <- table(f)
    if(!is.null(names(archetype))){
        if(anyNA(names(archetype)) || anyDuplicated(names(archetype))){
            stop("If 'archetype' is named, names must be non-NA and unqiue.",
                call. = FALSE)
        }
        archetype <- archetype[names(f_table)]
    }
    if(any(f_table < archetype)){
        stop("'archetype' out of bounds for some levels of 'f'. The maximum of",
            " 'archetype' is defined as table('f')", call. = FALSE)
    }
    if(length(archetype) == 1L){
        archetype <- rep(archetype,length(levels(f)))
    }
    archetype
}

# This function returns the index/position of rows/columns that are kept
# after merging.
#' @importFrom S4Vectors splitAsList
.get_element_pos <- function(f, archetype){
    archetype <- as.list(archetype)
    f_pos <- seq.int(1L, length(f))
    f_pos_split <- S4Vectors::splitAsList(f_pos, f)
    f_pos <- unlist(f_pos_split[archetype])
    f_pos
}

# This function merges assays and row/colData.
#' @importFrom S4Vectors SimpleList
#' @importFrom BiocParallel bpmapply SerialParam
.merge_rows_or_cols <- function(
        x, f, by, archetype = 1L, average = FALSE, BPPARAM = SerialParam(),
        check.assays = TRUE, na.rm = FALSE, ...){
    # input check
    if( !.is_a_bool(average) ){
        stop("'average' must be TRUE or FALSE.", call. = FALSE)
    }
    if( !.is_a_bool(check.assays) ){
        stop("'check.assays' must be TRUE or FALSE.", call. = FALSE)
    }
    if( !.is_a_bool(na.rm) ){
        stop("'na.rm' must be TRUE or FALSE.", call. = FALSE)
    }
    # Get correct functions based on whether we agglomerate rows or cols
    rowData_FUN <- switch(by, rowData, colData)
    nrow_FUN <- switch(by, nrow, ncol)
    rownames_FUN <- switch(by, rownames, colnames)
    rownames_ass_FUN <- switch(by, `rownames<-`, `colnames<-`)
    # If user specified column name from row/colData, get the values
    f <- .get_group_vector(x, f, by)
    # Check that the group ID vector is specifying groups for each element
    f <- .norm_f(nrow_FUN(x), f, ...)
    # If the data is already agglomerated at each group
    if(length(levels(f)) == nrow_FUN(x)){
        return(x)
    }
    # In merging, first element of certain group is kept by default. archetype,
    # can control this behavior; it can specify the preserved rows for every
    # group or index.
    archetype <- .norm_archetype(f, archetype)
    # Now we have agglomerated assays, but TreeSE has still the original form.
    # We take specified rows/columns from the TreeSE.
    idx <- .get_element_pos(f, archetype = archetype)
    # Retrieve experiment assays
    assays <- assays(x)
    # Merge assays
    assays <- bpmapply(
        .agglomerate_assay, assayNames(x), assays,
        MoreArgs = list(
            by = by, ids = f, na.rm = na.rm, average = average,
            check.assay = check.assays
        ),
        SIMPLIFY = FALSE, BPPARAM = BPPARAM
    )
    # Convert to SimpleList
    assays <- assays |> SimpleList()
    if( by == 1L ){
        x <- x[idx, ]
    } else{
        x <- x[ , idx]
    }
    # Add assays back to TreeSE
    assays(x, withDimnames = FALSE) <- assays
    # Change row/colnames. Currently, they have same names as in original data
    # but just certain rows. Change them to represent groups
    x <- rownames_ass_FUN(x, rownames_FUN(assays[[1]]))
    return(x)
}

# This function works similarly to scrapper::aggregateAcrossGenes but it returns
# same-class matrix instead of list of vectors
#' @importFrom DelayedArray DelayedArray type rowsum
.sum_counts_across_features <- function(x, ids, average, na.rm){
    # Remove rows with NA group from assay. This behavior is similar to what
    # aggregateAcrossGenes does; rows without a group are dropped.
    to_remove <- ids |> is.na() |> which()
    if( length(to_remove) != 0L ){
        # Remove NA groups
        ids <- ids[-to_remove]
        # Remove assay rows corresponding to NA groups
        x <- x[-to_remove, , drop = FALSE]
    }
    # Which cell is not NA?
    is_not_na <- !is.na(x)
    type(is_not_na) <- "integer"
    # Aggregate data to certain groups
    x <- rowsum(x, ids, na.rm = na.rm)
    # Calculate average if specified
    if( average ){
        x <- x/rowsum(is_not_na, ids)
    }
    return(x)
}

# This functions checks if assay has negative or binary values. It does not
# make sense to sum them, so we give warning to user.
.check_assay_for_merge <- function(assay.type, assay){
    # Check if assays include binary or negative values
    if( all(assay == 0 | assay == 1) ){
        warning("'", assay.type, "'", " includes binary values.",
                "\nAgglomeration of it might lead to meaningless values.",
                "\nCheck the assay, and consider doing transformation again",
                "manually with agglomerated data.",
                call. = FALSE)
    }
    if( !all( assay >= 0 | is.na(assay) ) ){
        warning("'", assay.type, "'", " includes negative values.",
                "\nAgglomeration of it might lead to meaningless values.",
                "\nCheck the assay, and consider doing transformation again",
                "manually with agglomerated data.",
                call. = FALSE)
    }
    return(assay)
}

#' @importFrom DelayedArray DelayedArray type rowsum
.agglomerate_assay <- function(
    assay.type, assay, by, ids, na.rm, average, check.assay
    ){
    # Check assay
    if( check.assay ){
        .check_assay_for_merge(assay.type, assay)
    }
    # Transpose if we are merging columns
    if( by == 2L ){
        assay <- t(assay)
    }
    # Sum counts across features
    assay <- .sum_counts_across_features(assay, ids, average, na.rm)
    # Transpose back to original orientation
    if( by == 2L ){
        assay <- t(assay)
    }
    return(assay)
}

# This function checks that sets are in correct format, i.e., there should
# be set information for each row/column. Columns should represent sets
# and rows features. Each cell specifies the membership of feature to the
# specific set.
.check_and_process_sets <- function(sets){
    # Check validity of set type
    if( !is.logical(sets) && !is.numeric(sets) ){
        stop("'group' variables must be numeric or logical.", call. = FALSE)
    }
    # Find NAs
    is_na <- which(is.na(sets))
    # Replace NAs
    if( length(is_na) != 0L ){
        warning(
            "NAs in 'group' variables were removed before agglomeration.",
            call. = FALSE
        )
        # Zero out NA sets
        sets[is_na] <- 0
    }
    return(sets)
}

# This function agglomerates an abundance table based on sets.
.agglomerate_set_assay <- function(assay.type, assay, by, sets, na.rm){
    # Check assay
    .check_assay_for_merge(assay.type, assay)
    # Replace NAs with 0
    if( na.rm ){
        assay[is.na(assay)] <- 0
    }
    # Transpose if we are merging columns
    if( by == 2L ){
        assay <- t(assay)
    }
    # Compute set-wise assay by cross-product
    assay <- as.matrix(crossprod(sets, assay))
    # Transpose back to original orientation
    if( by == 2L ){
        assay <- t(assay)
    }
    return(assay)
}

#' @importFrom Biostrings DNAStringSetList
.merge_refseq_list <- function(sequences_list, f, names, ...){
    threshold <- list(...)[["threshold"]]
    if(is.null(threshold)){
        threshold <- 0.05
    }
    if(!is(sequences_list,"DNAStringSetList")){
        return(.merge_refseq(sequences_list, f, threshold))
    }
    names <- names(sequences_list)
    seqs <- DNAStringSetList(lapply(sequences_list, .merge_refseq, f,
                                    threshold))
    names(seqs) <- names
    seqs
}

#' @importFrom Biostrings DNAStringSetList
#' @importFrom DECIPHER ConsensusSequence
.merge_refseq <- function(sequences, f, threshold){
    sequences <- split(sequences,f)
    seq <- unlist(DNAStringSetList(lapply(sequences, ConsensusSequence,
                                            threshold = threshold)))
    seq
}

# Merged rows (or columns) are named by their groups. However, if each row was
# its own group, the data was not merged, and rows keep their original names.
# This function returns the group of each row of the merged data.
.get_groups_of_merged <- function(names, f){
    groups <- names
    if( length(names) == length(f) ){
        groups <- as.character(f)
    }
    return(groups)
}

# This function returns the row or column tree. Agglomeration supports only a
# single tree. If the data has multiple trees, NULL is returned with a warning.
.get_single_tree <- function(x, by){
    tree_names <- switch(by, rowTreeNames, colTreeNames)(x)
    if( length(tree_names) > 1L ){
        warning(
            "The tree is removed because the data has multiple trees. Tree ",
            "agglomeration supports only a single tree.", call. = FALSE)
        return(NULL)
    }
    tree <- switch(by, rowTree, colTree)(x)
    return(tree)
}

# This function removes all row or column trees.
.remove_trees <- function(x, by){
    if( by == 1L ){
        rowTree(x) <- NULL
    } else{
        colTree(x) <- NULL
    }
    return(x)
}

# This function agglomerates the tree based on groups and adds it to the
# merged data. If the tree is not provided or it cannot be agglomerated, the
# tree is removed.
.update_tree_after_merge <- function(x, tree, tree_name, links, f, by){
    res <- NULL
    if( !is.null(tree) ){
        # Features without a group are not merged, so they are not included in
        # the tree either. Instead of feature names, nodes are identified by
        # node numbers.
        keep <- !is.na(f)
        nodes <- f[keep]
        names(nodes) <- links[["nodeNum"]][keep]
        res <- .agglomerate_tree(tree, nodes)
    }
    if( is.null(res) ){
        x <- .remove_trees(x, by)
        return(x)
    }
    # Link each merged row to the node of its group
    names <- switch(by, rownames, colnames)(x)
    groups <- .get_groups_of_merged(names, f)
    node_lab <- res[["links"]][ match(groups, names(res[["links"]])) ]
    x <- .set_tree(x, res[["tree"]], node_lab, tree_name, by)
    return(x)
}

# This function sets the tree. Each row (or column) i is linked to the node
# node_lab[i].
#' @importFrom TreeSummarizedExperiment changeTree
.set_tree <- function(x, tree, node_lab, tree_name, by){
    # All rows must be linked to the tree. Otherwise, changeTree() would drop
    # them.
    if( anyNA(node_lab) ){
        warning(
            "The tree is removed because some features could not be ",
            "linked to the tree.", call. = FALSE)
        x <- .remove_trees(x, by)
        return(x)
    }
    names(node_lab) <- switch(by, rownames, colnames)(x)
    # The data does not always have a tree. For example,
    # agglomerateByPrevalence() removes the tree before merging non-prevalent
    # features into "Other" group, and adds the tree back with this function.
    # If the data does not have a tree, changeTree() adds it with the default
    # name "phylo"; whichRowTree/whichColTree only selects an existing tree to
    # replace. That is why the original name is set explicitly.
    if( by == 1L ){
        x <- changeTree(x, rowTree = tree, rowNodeLab = node_lab)
        rowTreeNames(x) <- tree_name
    } else{
        x <- changeTree(x, colTree = tree, colNodeLab = node_lab)
        colTreeNames(x) <- tree_name
    }
    return(x)
}

#' @importFrom TreeSummarizedExperiment changeTree
.merge_rows_TSE <- function(x, f, archetype = 1L, update.tree = FALSE,
    update.refseq = mergeRefSeq, mergeRefSeq = FALSE, ...){
    # input check
    if(!.is_a_bool(update.tree)){
        stop("'update.tree' must be TRUE or FALSE.", call. = FALSE)
    }
    if(!.is_a_bool(update.refseq)){
        stop("'update.refseq' must be TRUE or FALSE.", call. = FALSE)
    }
    # for optionally merging referenceSeq and rowTree
    refSeq <- NULL
    if(update.refseq){
        refSeq <- referenceSeq(x)
    }
    # If user specified column name from rowData, get the values
    f <- .get_group_vector(x, f, by = 1L)
    row_links <- NULL
    tree <- NULL
    tree_name <- NULL
    # Tree can be updated only if there is a tree
    if( update.tree && !is.null(rowLinks(x)) ){
        row_links <- rowLinks(x)
        tree <- .get_single_tree(x, by = 1L)
        tree_name <- rowTreeNames(x)[[1L]]
        # Groups for the tree must match the groups used in merging
        f <- .norm_f(nrow(x), f, dim.type = "rows", ...)
    }
    #
    x <- .merge_rows_or_cols(x, f, by = 1L, archetype = 1L, ...)
    # Optionally merge rowTree. If it is not possible, the tree is removed.
    x <- .update_tree_after_merge(x, tree, tree_name, row_links, f, by = 1L)
    # optionally merge referenceSeq
    if(!is.null(refSeq)){
        referenceSeq(x) <- .merge_refseq_list(refSeq, f, rownames(x), ...)
    }
    x
}

#' @importFrom TreeSummarizedExperiment changeTree
.merge_cols_TSE <- function(x, f, archetype = 1L, update.tree = FALSE, ...){
    # input check
    if(!.is_a_bool(update.tree)){
        stop("'update.tree' must be TRUE or FALSE.", call. = FALSE)
    }
    #
    # If user specified column name from colData, get the values
    f <- .get_group_vector(x, f, by = 2L)
    col_links <- NULL
    tree <- NULL
    tree_name <- NULL
    # Tree can be updated only if there is a tree
    if( update.tree && !is.null(colLinks(x)) ){
        col_links <- colLinks(x)
        tree <- .get_single_tree(x, by = 2L)
        tree_name <- colTreeNames(x)[[1L]]
        # Groups for the tree must match the groups used in merging
        f <- .norm_f(ncol(x), f, dim.type = "columns", ...)
    }
    x <- .merge_rows_or_cols(x, f, by = 2L, archetype = 1L, ...)
    # Optionally merge colTree. If it is not possible, the tree is removed.
    x <- .update_tree_after_merge(x, tree, tree_name, col_links, f, by = 2L)
    return(x)
}
