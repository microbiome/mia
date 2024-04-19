
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
#' @importFrom rbiom unifrac
setMethod("calculateUnifrac", signature = c(x = "ANY", tree = "phylo"),
          function(x, tree, weighted = FALSE){
              if(is(x,"SummarizedExperiment")){
                  stop("When providing a 'tree', please provide a matrix-like as 'x'",
                       " and not a 'SummarizedExperiment' object. Please consider ",
                       "combining both into a 'TreeSummarizedExperiment' object.",
                       call. = FALSE) 
              }
              .calculate_distance(x, FUN = rbiom::unifrac, tree = tree, weighted = weighted)
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
                   tree_name = "phylo", ...){
              # Check assay.type and get assay
              .check_assay_present(assay.type, x)
              mat <- assay(x, assay.type)
              # Check tree_name
              .check_rowTree_present(tree_name, x)
              # Get tree
              tree <- rowTree(x, tree_name)
              
              calculateUnifrac(mat, tree = tree,...)
              
          }
)
