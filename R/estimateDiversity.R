
.estimate_faith <- function(
        x, mat, tree = NULL, tree.name = "phylo", node.label = node_lab,
        node_lab = NULL, ...){
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
            "The Faith's alpha diversity index is not possible to ",
            "calculate.", call. = FALSE)
    }

    # Check that node.label is NULL or it specifies links between rownames and
    # node labs
    if( !( is.null(node.label) ||
            is.character(node.label) && length(node.label) == nrow(x) ) ){
        stop(
            "'node.label' must be NULL or a vector specifying links between ",
            "rownames and node labs of 'tree'.",
            call. = FALSE)
    }

    # Subset rows of the assay to correspond node_labs (if there are any NAs
    # in node labels)
    if( !is.null(node.label) && any(is.na(node.label)) ){
        warning(
            "The tree named does not include all the ",
            "rows. 'x' is subsetted.", call. = FALSE)
        mat <- mat[ !is.na(node.label), ]
        node.label <- node.label[ !is.na(node.label) ]
    }
    # If there are node labels (any TreeSE should have because they are rowLinks
    # by default), rename the features in matrix to match with labels found in
    # tree.
    if( !is.null(node.label) ){
        rownames(mat) <- node.label
    }

    # To calculate faith, the assay must have rownames. TreeSE has always
    # rownames at this point, but if the object is SE, it might be that it is
    # missing rownames.
    if( is.null(rownames(mat)) ){
        stop("'x' must have rownames.", call. = FALSE)
    }
    # Calculates Faith index
    res <- .calc_faith(mat, tree, ...)
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
    temp <- reorder.phylo(tree, "cladewise")
    
    # Call the C++ code
    return(faith_cpp(mat, temp))
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
        shannon = .calc_shannon,
        gini_simpson = .calc_gini_simpson,
        inverse_simpson = .calc_inverse_simpson,
        coverage = .calc_coverage,
        fisher = .calc_fisher,
        faith = .estimate_faith,
        log_modulo_skewness = .calc_log_modulo_skewness
        )
    res <- FUN(x = x, mat = mat, ...)
    res <- unname(res)
    return(res)
}
