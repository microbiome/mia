.get_unifrac <- function(x, assay.type = assay_name, assay_name = exprs_values, 
                          exprs_values = "counts", tree.name = tree_name, 
                          tree_name = "phylo", ...){
        # Get functions and parameters based on direction
        tree_present_FUN <- if (transposed) .check_colTree_present
            else .check_rowTree_present
        tree_FUN <- if (transposed) colTree else rowTree
        links_FUN <- if (transposed) colLinks else rowLinks
        margin_name <- if (transposed) "col" else "row"
        # Check assay.type
        .check_assay_present(assay.type, x)
        # Check tree.name
        tree_present_FUN(tree.name, x)
        #
        # Select only those features/samples that are in the tree
        links <- links_FUN(x)
        present_in_tree <- links[, "whichTree"] == tree.name
        if( any(!present_in_tree) ){
            warning(
                "Not all ", margin_name, "s were present in the ", margin_name,
                "Tree specified by 'tree.name'. 'x' is subsetted.",
                call. = FALSE)
            # Subset the data
            if( transposed ){
                x <- x[, present_in_tree]
            } else{
                x <- x[present_in_tree, ]
            }
        }
        # Get assay and transpose it if specified. Features must be in columns
        # and samples in rows.
        mat <- assay(x, assay.type)
        # Get tree
        tree <- tree_FUN(x, tree.name)
        # Get links and take only nodeLabs
        links <- links_FUN(x)
        links <- links[ , "nodeLab" ]
        # Calculate unifrac
        res <- getUnifrac(x, tree = tree, node.label = links, ...)
        return(res)
}

################################################################################
#' @importFrom ape drop.tip
#' @importFrom rbiom unifrac
.calculate_unifrac <- function(
        x, tree, weighted = FALSE, node.label = nodeLab, nodeLab = NULL, ...){
    # Check x
    if( !is.matrix(as.matrix(x)) ){
        stop("'x' must be a matrix", call. = FALSE)
    }
    # input check
    if(!.is_a_bool(weighted)){
        stop("'weighted' must be TRUE or FALSE.", call. = FALSE)
    }
    if(is.null(colnames(x)) || is.null(rownames(x))){
        stop("colnames and rownames must not be NULL", call. = FALSE)
    }
    # node.label should be NULL or character vector specifying links between 
    # rows and tree labels
    if( !(is.null(node.label) ||
            (is.character(node.label) && length(node.label) == nrow(x) &&
            all(node.label[ !is.na(node.label) ] %in% c(tree$tip.label)))) ){
        stop(
            "'node.label' must be NULL or character specifying links between ",
            "abundance table and tree labels.", call. = FALSE)
    }
    # check that matrix and tree are compatible
    if( is.null(node.label) && !all(rownames(x) %in% c(tree$tip.label)) ) {
        stop(
            "Incompatible tree and abundance table! Please try to provide ",
            "'node.label'.", call. = FALSE)
    }
    # Merge rows, so that rows that are assigned to same tree node are agglomerated
    # together. If nodeLabs were provided, merge based on those. Otherwise merge
    # based on rownames
    if( is.null(node.label) ){
        node.label <- rownames(x)
    }
    # Prune tree if there are nodes that cannot be found from tips or if there
    # are tips that cannot be found from abundance matrix. It might be
    # that certain row is linked to internal node or that the tree has extra
    # tips that do not match with rows (e.g. after subsetting).
    if( any( !node.label %in% tree$tip.label ) ||
            any( !tree$tip.label %in% node.label) ){
        tree <- .prune_tree(tree, node.label, ...)
        warning("Pruning tree...", call. = FALSE)
    }
    # If node labels cannot be found from tips even after pruning, give error.
    # This kind of tree cannot be used in unifrac since it expects that every
    # row is linked to tips.
    if( any( !node.label %in% tree$tip.label ) ){
        stop(
            "Unifrac cannot be calculated since tree is not compatible. ",
            "Each row must be linked to tip of the tree.", call. = FALSE)
    }
    # Merge assay so that each row represent single tip. It might be that
    # multiple rows are linked to single tip.
    x <- .merge_assay_by_rows(x, node.label, ...)
    # Remove those tips that are not present in the data
    if( any(!tree$tip.label %in% rownames(x)) ){
        tree <- drop.tip(
            tree, tree$tip.label[!tree$tip.label %in% rownames(x)])
        warning(
            "The tree is pruned so that tips that cannot be found from ", 
            "the abundance matrix are removed.", call. = FALSE)
    }
    # Calculate unifrac. Use implementation from rbiom package
    res <- unifrac(x, tree = tree, weighted = weighted)
    return(res)
}

# Aggregate matrix based on nodeLabs. At the same time, rename rows based on node.label
# --> each row represent specific node of tree
#' @importFrom scuttle sumCountsAcrossFeatures
.merge_assay_by_rows <- function(x, node.label, average = FALSE, ...){
    if( !.is_a_bool(average) ){
        stop("'average' must be TRUE or FALSE.", call. = FALSE)
    }
    # Merge assay based on nodeLabs
    x <- sumCountsAcrossFeatures(
        x, ids = node.label, subset.row = NULL, subset.col = NULL,
        average = average)
    # Remove NAs from node.label
    node.label <- node.label[ !is.na(node.label) ]
    # Get the original order back
    x <- x[ node.label, ]
    return(x)
}

#' @importFrom ape is.rooted root
.norm_tree_to_be_rooted <- function(tree, names){
    if( !is.rooted(tree) ){
        randoroot <- sample(names, 1)
        warning("Randomly assigning root as -- ", randoroot, " -- in the",
                " phylogenetic tree in the data you provided.", call. = FALSE)
        tree <- root(phy = tree, outgroup = randoroot,
                     resolve.root = TRUE, interactive = FALSE)
        if( !is.rooted(tree) ){
            stop("Problem automatically rooting tree. Make sure your tree ",
                 "is rooted before attempting Unifrac calculation. See ",
                 "?ape::root", call. = FALSE)
        }
    }
    tree
}
