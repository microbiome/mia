# Helper function to extract tree and align matrix with tree nodes
.get_tree_and_linked_mat <- function(
        x, mat, tree = NULL, tree.name = "phylo", node.label = node_lab,
        node_lab = NULL, index_name = "phylogenetic", ...){
    # Input check
    # If tree is NULL, then we must have TreeSE object that has rowTree
    if( is.null(tree) &&
        !(is(x, "TreeSummarizedExperiment") && !is.null(rowTree(x))) ){
        stop("'tree' must be provided.", call. = FALSE)
    }
    # If tree is not specified, then we get rowTree
    if( is.null(tree) ){
        .check_rowTree_present(tree.name, x)
        tree <- rowTree(x, tree.name)
        # When we get rowTree, we know the linking between rows and nodes.
        node.label <- rowLinks(x)[ , "nodeLab" ]
        node.label[ rowLinks(x)[, "whichTree"] != tree.name ] <- NA
    }
    # Check that tree is correct
    if( is.null(tree) || is.null(tree$edge.length) ){
        stop(
            "'tree' is NULL or it does not have any branches. ",
            "The ", index_name, " alpha diversity index is not possible to ",
            "calculate.", call. = FALSE)
    }

    # Check that node.label is NULL or it specifies links between rownames and
    # node labs
    node_for_each <- is.character(node.label) &&
        length(node.label) == nrow(x) &&
        all(node.label[ !is.na(node.label) ] %in% tree$tip.label)
    named_vector <- is.character(node.label) && !is.null(names(node.label)) &&
        all(rownames(x) %in% names(node.label))
    if( !(is.null(node.label) || node_for_each || named_vector ) ){
        stop(
            "'node.label' must be NULL or character specifying links between ",
            "abundance table and tree labels.", call. = FALSE)
    }
    # If the labels were provided as named vector where names represent
    # original rows and values represent tips.
    if( named_vector ){
        node.label <- node.label[ match(rownames(x), names(node.label)) ]
    }

    # Subset rows of the assay to correspond node_labs (if there are any NAs
    # in node labels)
    if( !is.null(node.label) && any(is.na(node.label)) ){
        warning(
            "The tree named does not include all the ",
            "rows. 'x' is subsetted.", call. = FALSE)
        mat <- mat[ !is.na(node.label), , drop = FALSE]
        node.label <- node.label[ !is.na(node.label) ]
    }
    # If there are node labels (any TreeSE should have because they are rowLinks
    # by default), rename the features in matrix to match with labels found in
    # tree.
    if( !is.null(node.label) ){
        rownames(mat) <- node.label
    }

    # To calculate tree diversity, the assay must have rownames. TreeSE has
    # always rownames at this point, but if the object is SE, it might be that
    # it is missing rownames.
    if( is.null(rownames(mat)) ){
        stop("'x' must have rownames.", call. = FALSE)
    }
    return(list(mat = mat, tree = tree))
}

.estimate_faith <- function(x, mat, ...){
    temp <- .get_tree_and_linked_mat(x, mat, index_name = "Faith's", ...)
    args <- list(...)
    args <- args[ !names(args) %in% c(
        "tree", "tree.name", "node.label", "node_lab", "index_name", "index") ]
    res <- do.call(
        .calc_faith,
        c(list(mat = temp$mat, tree = temp$tree), args))
    return(res)
}

.estimate_allen <- function(x, mat, ...){
    temp <- .get_tree_and_linked_mat(x, mat, index_name = "Allen's", ...)
    args <- list(...)
    args <- args[ !names(args) %in% c(
        "tree", "tree.name", "node.label", "node_lab", "index_name", "index") ]
    res <- do.call(
        .calc_tree_diversity,
        c(list(mat = temp$mat, tree = temp$tree, index = "allen"), args))
    return(res)
}

.estimate_rao <- function(x, mat, ...){
    temp <- .get_tree_and_linked_mat(x, mat, index_name = "Rao's", ...)
    args <- list(...)
    args <- args[ !names(args) %in% c(
        "tree", "tree.name", "node.label", "node_lab", "index_name", "index") ]
    res <- do.call(
        .calc_tree_diversity,
        c(list(mat = temp$mat, tree = temp$tree, index = "rao"), args))
    return(res)
}

.calc_shannon <- function(mat, ...){
    vegan::diversity(t(mat), index="shannon")
}

# NOTE: vegan::diversity(x, index = "simpson")
# gives Simpson diversity, also called Gini-Simpson
# index: 1-lambda, where lambda is the Simpson index
# (lambda). This may cause confusion if your familiarity
# with diversity indices is limited.
# Moreover, Simpson's lambda is simply the
# squared sum of relative abundances so we can
# just use that for clarity and simplicity.
#.get_simpson <- function(x, ...){
.simpson_lambda <- function(mat, ...){

    # Convert table to relative values
    rel <- .calc_rel_abund(mat)

    # Squared sum of relative abundances
    colSums2(rel^2)
}

.calc_gini_simpson <- function(mat, ...){
    1 - .simpson_lambda(mat, ...)
}

.calc_inverse_simpson <- function(mat, ...){
    1 / .simpson_lambda(mat, ...)
}

.calc_coverage <- function(mat, threshold = 0.9, ...){

    # Threshold must be a numeric value between 0-1
    if( !( is.numeric(threshold) && (threshold >= 0 && threshold <= 1) ) ){
        stop("'threshold' must be a numeric value between 0-1.",
            call. = FALSE)
    }

    # Convert table to relative values
    rel <- .calc_rel_abund(mat)

    # Number of groups needed to have threshold (e.g. 50 %) of the
    # ecosystem occupied
    coverage <- apply(rel, 2, function(x) {
        min(which(cumsum(rev(sort(x/sum(x)))) >= threshold))
    })
    names(coverage) <- colnames(rel)
    coverage
}

.calc_fisher <- function(mat, ...){
    vegan::fisher.alpha(t(mat))
}

# These tags are required to enable the use of Rcpp in the package
#' @useDynLib mia
#' @importFrom Rcpp sourceCpp
NULL

#' @importFrom ape reorder.phylo
.calc_faith <- function(mat, tree, only.tips = FALSE, ...){
    # Input check
    if( !.is_a_bool(only.tips) ){
        stop("'only.tips' must be TRUE or FALSE.", call. = FALSE)
    }
    # Remove internal nodes if specified
    if( only.tips ){
        mat <- mat[ rownames(mat) %in% tree$tip.label, ]
    }
    
    # To ensure that the function works with NA also, convert NAs to 0.
    # Zero means that the taxon is not present --> same as NA (no information)
    mat[ is.na(mat) ] <- 0
    
    # The tree must be in cladewise order for the algorithm to work correctly
    tree <- reorder.phylo(tree, "cladewise")
    
    # Call the C++ code
    return(.faith_cpp(mat, tree))
}

.calc_tree_diversity <- function(mat, tree, index = c("allen", "rao"), ...){
    index <- match.arg(index)
    # Ensure NAs are 0
    mat[ is.na(mat) ] <- 0
    # Relative abundances per sample
    rel <- .calc_rel_abund(mat)
    # Ensure matrix format
    if( is.vector(rel) ){
        rel <- matrix(
            rel, ncol = 1, dimnames = list(rownames(mat), colnames(mat)))
    }
    # In case any sample has sum 0, rel can have NaNs
    rel[ is.nan(rel) ] <- 0

    # Reorder tree in postorder for single-pass bottom-up accumulation
    tree <- reorder.phylo(tree, "postorder")
    n_tips <- length(tree$tip.label)
    n_nodes <- tree$Nnode
    n_samples <- ncol(rel)

    node_abund <- matrix(0, nrow = n_tips + n_nodes, ncol = n_samples)
    m <- match(tree$tip.label, rownames(rel))
    valid_tips <- !is.na(m)
    node_abund[which(valid_tips), ] <- rel[m[valid_tips], , drop = FALSE]

    edge <- tree$edge
    edge_len <- tree$edge.length

    if( any(is.na(edge_len)) ){
        stop("'tree$edge.length' contains NA values.", call. = FALSE)
    }

    res <- numeric(n_samples)

    for( i in seq_len(nrow(edge)) ){
        child <- edge[i, 2]
        parent <- edge[i, 1]
        a_i <- node_abund[child, ]
        L_i <- edge_len[i]

        # Accumulate descendant abundance to parent node
        node_abund[parent, ] <- node_abund[parent, ] + a_i

        # Calculate contribution
        pos <- a_i > 0
        if( any(pos) ){
            if( index == "allen" ){
                # Allen: - L_i * a_i * log(a_i)
                res[pos] <- res[pos] - L_i * a_i[pos] * log(a_i[pos])
            } else if( index == "rao" ){
                # Rao: 2 * L_i * a_i * (1 - a_i)
                res[pos] <- res[pos] + 2 * L_i * a_i[pos] * (1 - a_i[pos])
            }
        }
    }
    names(res) <- colnames(mat)
    return(res)
}

.calc_log_modulo_skewness <- function(mat, quantile = 0.5,
    nclasses = num_of_classes, num_of_classes = 50, ...){
    # quantile must be a numeric value between 0-1
    if( !( is.numeric(quantile) && (quantile >= 0 && quantile <= 1) ) ){
        stop("'quantile' must be a numeric value between 0-1.",
            call. = FALSE)
    }
    # nclasses must be a positive numeric value
    if( !( is.numeric(nclasses) && nclasses > 0 ) ){
        stop("'nclasses' must be a positive numeric value.",
            call. = FALSE)
    }
    # Determine the quantile point.
    quantile_point <- quantile(max(mat), quantile)
    # Tabulate the arithmetic abundance classes. Use the same classes
    # for all samples for consistency
    cutpoints <- c(seq(0, quantile_point, length=nclasses), Inf)
    # Calculates sample-wise frequencies. How many taxa in each interval?
    freq_table <- table(cut(mat, cutpoints), col(mat))
    # Calculates the skewness of frequency table. Returns skewness for each
    # sample
    r <- .calc_skewness(freq_table)
    # Return log-modulo
    log(1 + r)
}

#' @importFrom DelayedMatrixStats rowSums2 rowMeans2
.calc_skewness <- function(x) {
    # Transposes the table
    x <- t(x)
    # Each value is subtracted by sample-wise mean, which is raised to the
    # power of 3.
    # Then the sample-wise sum is taken from these values.
    numerator <- rowSums2((x - rowMeans2(x))^3)
    # Sample-wise sum is divided by number of taxa that are not NA.
    numerator <- numerator/rowSums2(!is.na(x))
    # Each value is subtracted by sample-wise mean, which is raises to the
    # power of 2.
    # Then the sample-wise sum is taken from these values.
    denominator <- rowSums2((x - rowMeans2(x))^2)
    # Sample-wise sum is divided by number of taxa that are not NA. Then
    # these values
    # are raised to the power of 3/2.
    denominator <- (denominator/rowSums2(!is.na(x)))^(3/2)
    # Result
    result <- numerator/denominator
    return(result)
}

#' @importFrom SummarizedExperiment assay assays
.get_diversity_values <- function(index, x, mat, ...){
    FUN <- switch(index,
        allen = .estimate_allen,
        shannon = .calc_shannon,
        gini_simpson = .calc_gini_simpson,
        inverse_simpson = .calc_inverse_simpson,
        coverage = .calc_coverage,
        fisher = .calc_fisher,
        faith = .estimate_faith,
        log_modulo_skewness = .calc_log_modulo_skewness,
        rao = .estimate_rao
        )
    res <- FUN(x = x, mat = mat, ...)
    res <- unname(res)
    return(res)
}
