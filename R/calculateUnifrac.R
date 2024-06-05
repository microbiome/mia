#' Calculate weighted or unweighted (Fast) Unifrac distance
#'
#' This function calculates the Unifrac distance for all sample-pairs
#' in a \code{\link[TreeSummarizedExperiment:TreeSummarizedExperiment-class]{TreeSummarizedExperiment}}
#' object. The function utilizes \code{\link[rbiom:unifrac]{rbiom:unifrac()}}.
#'
#' Please note that if \code{calculateUnifrac} is used as a \code{FUN} for
#' \code{runMDS}, the argument \code{ntop} has to be set to \code{nrow(x)}.
#'
#' @param x a numeric matrix or a
#'   \code{\link[TreeSummarizedExperiment:TreeSummarizedExperiment-class]{TreeSummarizedExperiment}}
#'   object containing a tree.
#'
#'   Please  note that \code{runUnifrac} expects a matrix with samples per row
#'   and not per column. This is implemented to be compatible with other
#'   distance calculations such as \code{\link[stats:dist]{dist}} as much as
#'   possible.
#'
#' @param tree if \code{x} is a matrix, a
#'   \code{\link[TreeSummarizedExperiment:phylo]{phylo}} object matching the
#'   matrix. This means that the phylo object and the columns should relate
#'   to the same type of features (aka. microorganisms).
#'   
#' @param node.lab if \code{x} is a matrix, 
#'   a \code{character} vector specifying links between rows/columns and tips of \code{tree}.
#'   The length must equal the number of rows/columns of \code{x}. Furthermore, all the 
#'   node labs must be present in \code{tree}.
#'   
#' @param nodeLab Deprecated. Use \code{node.lab} instead.
#'
#' @param assay.type a single \code{character} value for specifying which
#'   assay to use for calculation.
#'   
#' @param exprs_values a single \code{character} value for specifying which
#'   assay to use for calculation.
#'   (Please use \code{assay.type} instead.)
#'   
#' @param assay_name a single \code{character} value for specifying which
#'   assay to use for calculation.
#'   (Please use \code{assay.type} instead. At some point \code{assay_name}
#'   will be disabled.)
#'
#' @param tree.name a single \code{character} value for specifying which
#'   tree will be used in calculation. 
#'   (By default: \code{tree.name = "phylo"})
#'   
#' @param tree_name Deprecated. Use \code{tree.name} instead.
#'   
#' @param weighted \code{TRUE} or \code{FALSE}: Should use weighted-Unifrac
#'   calculation? Weighted-Unifrac takes into account the relative abundance of
#'   species/taxa shared between samples, whereas unweighted-Unifrac only
#'   considers presence/absence. Default is \code{FALSE}, meaning the
#'   unweighted-Unifrac distance is calculated for all pairs of samples.
#'
#' @param transposed Logical scalar, is x transposed with cells in rows, i.e., 
#'   is Unifrac distance calculated based on rows (FALSE) or columns (TRUE).
#'   (By default: \code{transposed = FALSE})
#'   
#' @param ... optional arguments not used.
#'
#' @return a sample-by-sample distance matrix, suitable for NMDS, etc.
#'
#' @references
#' \url{http://bmf.colorado.edu/unifrac/}
#'
#' See also additional descriptions of Unifrac in the following articles:
#'
#' Lozupone, Hamady and Knight, ``Unifrac - An Online Tool for Comparing
#' Microbial Community Diversity in a Phylogenetic Context.'', BMC
#' Bioinformatics 2006, 7:371
#'
#' Lozupone, Hamady, Kelley and Knight, ``Quantitative and qualitative (beta)
#' diversity measures lead to different insights into factors that structure
#' microbial communities.'' Appl Environ Microbiol. 2007
#'
#' Lozupone C, Knight R. ``Unifrac: a new phylogenetic method for comparing
#' microbial communities.'' Appl Environ Microbiol. 2005 71 (12):8228-35.
#'
#' @name calculateUnifrac
#'
#' @export
#'
#' @examples
#' data(esophagus)
#' library(scater)
#' calculateUnifrac(esophagus, weighted = FALSE)
#' calculateUnifrac(esophagus, weighted = TRUE)
#' # for using calculateUnifrac in conjunction with runMDS the tree argument
#' # has to be given separately. In addition, subsetting using ntop must
#' # be disabled
#' esophagus <- runMDS(esophagus, FUN = calculateUnifrac, name = "Unifrac",
#'                     tree = rowTree(esophagus),
#'                     assay.type = "counts",
#'                     ntop = nrow(esophagus))
#' reducedDim(esophagus)
NULL

#' @rdname calculateUnifrac
#' @export
setGeneric("calculateUnifrac", signature = c("x", "tree"),
            function(x, tree, ... )
                standardGeneric("calculateUnifrac"))

#' @rdname calculateUnifrac
#' @export
setMethod("calculateUnifrac", signature = c(x = "ANY", tree = "phylo"),
    function(x, tree, weighted = FALSE, ...){
        if(is(x,"SummarizedExperiment")){
            stop("When providing a 'tree', please provide a matrix-like as 'x'",
                " and not a 'SummarizedExperiment' object. Please consider ",
                "combining both into a 'TreeSummarizedExperiment' object.",
                call. = FALSE) 
        }
        .calculate_distance(x, FUN = runUnifrac, tree = tree,
                            weighted = weighted, ...)
    }
)

#' @rdname calculateUnifrac
#'
#' @importFrom SummarizedExperiment assay
#'
#' @export
setMethod("calculateUnifrac",
        signature = c(x = "TreeSummarizedExperiment",
                    tree = "missing"),
    function(x, assay.type = assay_name, assay_name = exprs_values, exprs_values = "counts", 
             tree.name = tree_name, tree_name = "phylo", transposed = FALSE, ...){
        # Check assay.type and get assay
        .check_assay_present(assay.type, x)
        mat <- assay(x, assay.type)
        if(!transposed){
            # Check tree.name
            .check_rowTree_present(tree.name, x)
            # Get tree
            tree <- rowTree(x, tree.name)
            # Select only those features that are in the rowTree
            whichTree <- rowLinks(x)[, "whichTree"] == tree.name
            if( any(!whichTree) ){
                warning("Not all rows were present in the rowTree specified by 'tree.name'.",
                        "'x' is subsetted.", call. = FALSE)
                # Subset the data
                x <- x[ whichTree, ]
                mat <- mat[ whichTree, ]
            }
            mat <- t(mat)
            tree <- .norm_tree_to_be_rooted(tree, rownames(x))
            # Get links
            links <- rowLinks(x)
        } else {
            # Check tree.name
            .check_colTree_present(tree.name, x)
            # Get tree
            tree <- colTree(x, tree.name)
            # Select only those samples that are in the colTree
            whichTree <- colLinks(x)[, "whichTree"] == tree.name
            if( any(!whichTree) ){
                warning("Not all columns were present in the colTree specified by 'tree.name'.",
                        "'x' is subsetted.", call. = FALSE)
                # Subset the data
                x <- x[ , whichTree ]
                mat <- mat[ , whichTree ]
            }
            tree <- .norm_tree_to_be_rooted(tree, colnames(x))
            # Get links
            links <- colLinks(x)
        }
        # Remove those links (make them NA) that are not included in this tree
        links[ links$whichTree != tree.name, ] <- NA
        # Take only nodeLabs
        links <- links[ , "nodeLab" ]
        calculateUnifrac(mat, tree = tree, node.lab = links, ...)
    }
)

################################################################################
# Fast Unifrac for R.
# Adapted from The ISME Journal (2010) 4, 17-27; doi:10.1038/ismej.2009.97
#
# adopted from original implementation in phyloseq implemented by
# Paul J. McMurdie (https://github.com/joey711/phyloseq)
################################################################################
#' @rdname calculateUnifrac
#'
#' @importFrom ape drop.tip
#' @importFrom rbiom unifrac
#' @export
runUnifrac <- function(x, tree, weighted = FALSE, normalized = TRUE,
    node.lab = nodeLab, nodeLab = NULL, BPPARAM = SerialParam(), ...){
    # Check x
    if( !is.matrix(as.matrix(x)) ){
        stop("'x' must be a matrix", call. = FALSE)
    }
    # x has samples as row. Therefore transpose. This benchmarks faster than
    # converting the function to work with the input matrix as is
    x <- try(t(x), silent = TRUE)
    if(is(x,"try-error")){
        stop("The input to 'runUnifrac' must be a matrix-like object: ", 
             as.character(x), call. = FALSE)
    }
    # input check
    if(!.is_a_bool(weighted)){
        stop("'weighted' must be TRUE or FALSE.", call. = FALSE)
    }
    if(is.null(colnames(x)) || is.null(rownames(x))){
        stop("colnames and rownames must not be NULL", call. = FALSE)
    }
    # node.lab should be NULL or character vector specifying links between 
    # rows and tree labels
    if( !(is.null(node.lab) ||
        (is.character(node.lab) && length(node.lab) == nrow(x) &&
        all(node.lab[ !is.na(node.lab) ] %in% c(tree$tip.label)))) ){
        stop("'node.lab' must be NULL or character specifying links between ",
             "abundance table and tree labels.", call. = FALSE)
    }
    # check that matrix and tree are compatible
    if( is.null(node.lab) && 
        !all(rownames(x) %in% c(tree$tip.label)) ) {
        stop("Incompatible tree and abundance table! Please try to provide ",
             "'node.lab'.", call. = FALSE)
    }
    # Merge rows, so that rows that are assigned to same tree node are agglomerated
    # together. If nodeLabs were provided, merge based on those. Otherwise merge
    # based on rownames
    if( is.null(node.lab) ){
        node.lab <- rownames(x)
    }
    # Prune tree if there are nodes that cannot be found from tips or if there
    # are tips that cannot be found from abundance matrix. It might be
    # that certain row is linked to internal node or that the tree has extra
    # tips that do not match with rows (e.g. after subsetting).
    if( any( !node.lab %in% tree$tip.label ) ||
            any( !tree$tip.label %in% node.lab) ){
        tree <- .prune_tree(tree, node.lab)
        warning("Pruning tree...", call. = FALSE)
    }
    # If node labels cannot be found from tips even after pruning, give error.
    # This kind of tree cannot be used in unifrac since it expects that every
    # row is linked to tips.
    if( any( !node.lab %in% tree$tip.label ) ){
        stop(
            "Unifrac cannot be calculated since tree is not compatible. ",
            "Each row must be linked to tip of the tree.", call. = FALSE)
    }
    # Merge assay so that each row represent single tip. It might be that
    # multiple rows are linked to single tip.
    x <- .merge_assay_by_rows(x, node.lab, ...)
    # Modify tree so that it will become rooted.
    tree <- .norm_tree_to_be_rooted(tree, rownames(x))
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

# Aggregate matrix based on nodeLabs. At the same time, rename rows based on node.lab
# --> each row represent specific node of tree
#' @importFrom scuttle sumCountsAcrossFeatures
.merge_assay_by_rows <- function(x, node.lab, average = FALSE, ...){
    if( !.is_a_bool(average) ){
        stop("'average' must be TRUE or FALSE.", call. = FALSE)
    }
    # Merge assay based on nodeLabs
    x <- sumCountsAcrossFeatures(
        x, ids = node.lab, subset.row = NULL, subset.col = NULL,
        average = average)
    # Remove NAs from node.lab
    node.lab <- node.lab[ !is.na(node.lab) ]
    # Get the original order back
    x <- x[ node.lab, ]
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
