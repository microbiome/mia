#' Calculate weighted or unweighted (Fast) Unifrac distance
#'
#' This function calculates the (Fast) Unifrac distance for all sample-pairs
#' in a \code{\link[TreeSummarizedExperiment:TreeSummarizedExperiment-class]{TreeSummarizedExperiment}}
#' object.
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
#' @param nodeLab if \code{x} is a matrix, 
#'   a \code{character} vector specifying links between rows/columns and nodes of \code{tree}.
#'   The length must equal the number of rows/columns of \code{x}. Furthermore, all the 
#'   node labs must be present in \code{tree}.
#'
#' @param assay_name a single \code{character} value for specifying which
#'   assay to use for calculation.
#'   
#' @param exprs_values a single \code{character} value for specifying which
#'   assay to use for calculation.
#'   (Please use \code{assay_name} instead.)
#'   
#' @param abund_values a single \code{character} value for specifying which
#'   assay to use for calculation.
#'   (Please use \code{assay_name} instead. At some point \code{abund_values}
#'   will be disabled.)
#'
#' @param tree_name a single \code{character} value for specifying which
#'   tree will be used in calculation. 
#'   (By default: \code{tree_name = "phylo"})
#'   
#' @param weighted \code{TRUE} or \code{FALSE}: Should use weighted-Unifrac
#'   calculation? Weighted-Unifrac takes into account the relative abundance of
#'   species/taxa shared between samples, whereas unweighted-Unifrac only
#'   considers presence/absence. Default is \code{FALSE}, meaning the
#'   unweighted-Unifrac distance is calculated for all pairs of samples.
#'
#' @param normalized \code{TRUE} or \code{FALSE}: Should the output be
#'   normalized such that values range from 0 to 1 independent of branch length
#'   values? Default is \code{TRUE}. Note that (unweighted) \code{Unifrac} is
#'   always normalized by total branch-length, and so this value is ignored when
#'   \code{weighted == FALSE}.
#'
#' @param BPPARAM A
#'   \code{\link[BiocParallel:BiocParallelParam-class]{BiocParallelParam}}
#'   object specifying whether the Unifrac calculation should be parallelized.
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
#' The main implementation (Fast Unifrac) is adapted from the algorithm's
#' description in:
#'
#' Hamady, Lozupone, and Knight,
#' ``\href{http://www.nature.com/ismej/journal/v4/n1/full/ismej200997a.html}{Fast
#' UniFrac:} facilitating high-throughput phylogenetic analyses of microbial
#' communities including analysis of pyrosequencing and PhyloChip data.'' The
#' ISME Journal (2010) 4, 17--27.
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
#' @author
#' Paul J. McMurdie.
#' Adapted for mia by Felix G.M. Ernst
#'
#' @examples
#' data(esophagus)
#' library(scater)
#' calculateUnifrac(esophagus, weighted = FALSE)
#' calculateUnifrac(esophagus, weighted = TRUE)
#' calculateUnifrac(esophagus, weighted = TRUE, normalized = FALSE)
#' # for using calculateUnifrac in conjunction with runMDS the tree argument
#' # has to be given separately. In addition, subsetting using ntop must
#' # be disabled
#' esophagus <- runMDS(esophagus, FUN = calculateUnifrac, name = "Unifrac",
#'                     tree = rowTree(esophagus),
#'                     exprs_values = "counts",
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
    function(x, tree, weighted = FALSE, normalized = TRUE,
             BPPARAM = SerialParam(), ...){
        if(is(x,"SummarizedExperiment")){
           stop("When providing a 'tree', please provide a matrix-like as 'x'",
                " and not a 'SummarizedExperiment' object. Please consider ",
                "combining both into a 'TreeSummarizedExperiment' object.",
                call. = FALSE) 
        }
        .calculate_distance(x, FUN = runUnifrac, tree = tree,
                            weighted = weighted, normalized = normalized,
                            BPPARAM = BPPARAM, ...)
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
    function(x, assay_name = abund_values, abund_values = exprs_values, exprs_values = "counts", 
             tree_name = "phylo", transposed = FALSE, ...){
        # Check assay_name and get assay
        .check_assay_present(assay_name, x)
        mat <- assay(x, assay_name)
        if(!transposed){
            # Check tree_name
            .check_rowTree_present(tree_name, x)
            # Get tree
            tree <- rowTree(x, tree_name)
            # Select only those features that are in the rowTree
            whichTree <- rowLinks(x)[, "whichTree"] == tree_name
            if( any(!whichTree) ){
                warning("Not all rows were present in the rowTree specified by 'tree_name'.",
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
            # Check tree_name
            .check_colTree_present(tree_name, x)
            # Get tree
            tree <- colTree(x, tree_name)
            # Select only those samples that are in the colTree
            whichTree <- colLinks(x)[, "whichTree"] == tree_name
            if( any(!whichTree) ){
                warning("Not all columns were present in the colTree specified by 'tree_name'.",
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
        links[ links$whichTree != tree_name, ] <- NA
        # Take only nodeLabs
        links <- links[ , "nodeLab" ]
        calculateUnifrac(mat, tree = tree, nodeLab = links, ...)
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
#' @importFrom ape prop.part reorder.phylo node.depth node.depth.edgelength
#' @importFrom utils combn
#' @importFrom stats as.dist
#' @importFrom BiocParallel SerialParam register bplapply bpisup bpstart bpstop
#' @importFrom DelayedArray getAutoBPPARAM setAutoBPPARAM
#'
#' @export
runUnifrac <- function(x, tree, weighted = FALSE, normalized = TRUE,
                       nodeLab = NULL, BPPARAM = SerialParam(), ...){
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
    if(!.is_a_bool(normalized)){
        stop("'normalized' must be TRUE or FALSE.", call. = FALSE)
    }
    if(is.null(colnames(x)) || is.null(rownames(x))){
        stop("colnames and rownames must not be NULL", call. = FALSE)
    }
    # nodeLab should be NULL or character vector specifying links between 
    # rows and tree labels
    if( !(is.null(nodeLab) ||
        (is.character(nodeLab) && length(nodeLab) == nrow(x) &&
        all(nodeLab[ !is.na(nodeLab) ] %in% c(tree$tip.label, tree$node.label)))) ){
        stop("'nodeLab' must be NULL or character specifying links between ",
             "abundance table and tree labels.", call. = FALSE)
    }
    # check that matrix and tree are compatible
    if( is.null(nodeLab) && 
        !all(rownames(x) %in% c(tree$tip.label, tree$node.label)) ) {
        stop("Incompatible tree and abundance table! Please try to provide ",
             "'nodeLab'.", call. = FALSE)
    }
    # Merge rows, so that rows that are assigned to same tree node are agglomerated
    # together. If nodeLabs were provided, merge based on those. Otherwise merge
    # based on rownames
    if( is.null(nodeLab) ){
        nodeLab <- rownames(x)
    }
    # Merge assay
    x <- .merge_assay_by_rows(x, nodeLab, ...)
    # Modify tree
    tree <- .norm_tree_to_be_rooted(tree, rownames(x))
    #
    old <- getAutoBPPARAM()
    setAutoBPPARAM(BPPARAM)
    on.exit(setAutoBPPARAM(old))
    if (!(bpisup(BPPARAM) || is(BPPARAM, "MulticoreParam"))) {
        bpstart(BPPARAM)
        on.exit(bpstop(BPPARAM), add = TRUE)
    }
    #
    # create N x 2 matrix of all pairwise combinations of samples.
    spn <- utils::combn(colnames(x), 2, simplify = FALSE)

    ########################################
    # Build the requisite matrices as defined
    # in the Fast Unifrac article.
    ########################################
    ## This only needs to happen once in a call to Unifrac.
    ## Notice that A and B do not appear in this section.
    # Begin by building the edge descendants matrix (edge-by-sample)
    # `edge_array`
    #
    
    # Remove those tips that are not present in the data
    tree <- ape::drop.tip(
        tree, tree$tip.label[!tree$tip.label %in% rownames(x)])
    # Create a list of descendants, starting from the first internal node (root)
    ntip <- length(tree$tip.label)
    # Create a matrix that maps each internal node to its 2 descendants
    # This matrix doesn't include the tips, so must use node#-ntip to index into
    # it
    ## to suppress the warning if the number of node edges is uneven, we add the
    ## first node again
    node.desc <- tree$edge[order(tree$edge[,1]),][,2]
    if(length(node.desc) %% 2 == 1){
        node.desc <- c(node.desc,node.desc[1])
    }
    node.desc <- matrix(node.desc, byrow = TRUE, ncol = 2)
    # Define the edge_array object
    # Right now this is a node_array object, each row is a node (including tips)
    # It will be subset and ordered to match tree$edge later
    edge_array <- matrix(0, nrow = ntip+tree$Nnode, ncol = ncol(x),
                         dimnames = list(NULL, sample_names = colnames(x)))
    # Load the tip counts in directly
    x_tip <- x[ rownames(x) %in% tree$tip.label, ]
    edge_array[seq_len(ntip),] <- x_tip
    # Get a list of internal nodes ordered by increasing depth
    ord.node <- order(node.depth(tree))[(ntip+1):(ntip+tree$Nnode)]
    # Loop over internal nodes, summing their descendants to get that nodes
    # count
    for(i in ord.node){
        edge_array[i,] <- colSums(edge_array[node.desc[i-ntip,], , drop=FALSE],
                                  na.rm = TRUE)
    }
    # Keep only those with a parental edge (drops root) and order to match
    # tree$edge
    edge_array <- edge_array[tree$edge[,2],]
    # calculate the sums per sample
    samplesums <- colSums(x)
    # Remove unneeded variables.
    rm(node.desc)
    ############################################################################
    # calculate the distances
    ############################################################################
    if(weighted){
        if(!normalized){
            distlist <- BiocParallel::bplapply(spn,
                                               unifrac_weighted_not_norm,
                                               tree = tree,
                                               samplesums = samplesums,
                                               edge_array = edge_array,
                                               BPPARAM = BPPARAM)
        } else {
            # This is only relevant to weighted-Unifrac.
            # For denominator in the normalized distance, we need the age of each
            # tip.
            # 'z' is the tree in postorder order used in calls to .C
            # Descending order of left-hand side of edge (the ancestor to the node)
            z <- ape::reorder.phylo(tree, order="postorder")
            # Call phyloseq-internal function that in-turn calls ape's internal
            # horizontal position function, in C, using the re-ordered phylo object,
            # `z`
            tipAges = node.depth.edgelength(tree)
            # Keep only the tips, and add the tip labels in case `z` order differs
            # from `tree`
            tipAges <- tipAges[seq.int(1L, length(tree$tip.label))]
            names(tipAges) <- z$tip.label
            # Explicitly re-order tipAges to match x
            tipAges <- tipAges[rownames(x)]
            distlist <- BiocParallel::bplapply(spn,
                                               unifrac_weighted_norm,
                                               mat = x,
                                               tree = tree,
                                               samplesums = samplesums,
                                               edge_array = edge_array,
                                               tipAges = tipAges,
                                               BPPARAM = BPPARAM)
        }
    } else {
        # For unweighted Unifrac, convert the edge_array to an occurrence
        # (presence/absence binary) array
        edge_occ <- (edge_array > 0) - 0
        distlist <- BiocParallel::bplapply(spn,
                                           unifrac_unweighted,
                                           tree = tree,
                                           samplesums = samplesums,
                                           edge_occ = edge_occ,
                                           BPPARAM = BPPARAM)
    }
    # Initialize UnifracMat with NAs
    UnifracMat <- matrix(NA_real_, ncol(x), ncol(x))
    rownames(UnifracMat) <- colnames(UnifracMat) <- colnames(x)
    # Matrix-assign lower-triangle of UnifracMat. Then coerce to dist and
    # return.
    matIndices <- matIndices <- matrix(c(vapply(spn,"[",character(1),2L),
                                         vapply(spn,"[",character(1),1L)),
                                       ncol = 2)
    UnifracMat[matIndices] <- unlist(distlist)
    #
    stats::as.dist(UnifracMat)
}

# Aggregate matrix based on nodeLabs. At the same time, rename rows based on nodeLab
# --> each row represent specific node of tree
.merge_assay_by_rows <- function(x, nodeLab, average = FALSE, ...){
    if( !.is_a_bool(average) ){
        stop("'average' must be TRUE or FALSE.", call. = FALSE)
    }
    # Merge assay based on nodeLabs
    x <- scuttle::sumCountsAcrossFeatures(x, ids = nodeLab, 
                                          subset.row = NULL, subset.col = NULL, 
                                          average = average)
    # Remove NAs from nodeLab
    nodeLab <- nodeLab[ !is.na(nodeLab) ]
    # Get the original order back
    x <- x[ nodeLab, ]
    return(x)
}

unifrac_unweighted <- function(i, tree, samplesums, edge_occ){
    A  <- i[1]
    B  <- i[2]
    AT <- samplesums[A]
    BT <- samplesums[B]
    # Unweighted Unifrac
    # Subset matrix to just columns A and B
    edge_occ_AB <- edge_occ[, c(A, B)]
    edge_occ_AB_rS <- rowSums(edge_occ_AB, na.rm = TRUE)
    # Keep only the unique branches. Sum the lengths
    edge_uni_AB_sum <- sum((tree$edge.length * edge_occ_AB)[edge_occ_AB_rS < 2,],
                           na.rm=TRUE)
    # Normalize this sum to the total branches among these two samples, A and
    # B
    uwUFpairdist <- edge_uni_AB_sum /
        sum(tree$edge.length[edge_occ_AB_rS > 0])
    uwUFpairdist
}
# if not-normalized weighted Unifrac, just return "numerator";
# the u-value in the w-Unifrac description
unifrac_weighted_not_norm <- function(i, tree, samplesums, edge_array){
    A  <- i[1]
    B  <- i[2]
    AT <- samplesums[A]
    BT <- samplesums[B]
    # weighted Unifrac
    wUF_branchweight <- abs(edge_array[, A]/AT - edge_array[, B]/BT)
    # calculate the w-UF numerator
    numerator <- sum((tree$edge.length * wUF_branchweight), na.rm = TRUE)
    # if not-normalized weighted Unifrac, just return "numerator";
    # the u-value in the w-Unifrac description
    numerator
}
unifrac_weighted_norm <- function(i, mat, tree, samplesums, edge_array,
                                  tipAges){
    A  <- i[1]
    B  <- i[2]
    AT <- samplesums[A]
    BT <- samplesums[B]
    # weighted Unifrac
    wUF_branchweight <- abs(edge_array[, A]/AT - edge_array[, B]/BT)
    # calculate the w-UF numerator
    numerator <- sum((tree$edge.length * wUF_branchweight), na.rm = TRUE)
    # denominator (assumes tree-indices and matrix indices are same order)
    denominator <- sum((tipAges * (mat[, A]/AT + mat[, B]/BT)), na.rm = TRUE)
    # return the normalized weighted Unifrac values
    numerator / denominator
}

################################################################################

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
