#' Clustering wrapper
#'
#' This function returns a \code{SummarizedExperiment} with clustering 
#'   information in its colData or rowData
#'
#' @param x A
#'   \code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#'   object.
#'
#' @param BLUSPARAM A 
#'   \code{\link[https://rdrr.io/github/LTLA/bluster/man/BlusterParam-class.html]{BlusterParam}}
#'   object.
#'   
#' @param altExp A single character value indicating which alternate
#'   experiment to use for the clustering (if any). If the clustering is done
#'   on the main experiment, set \code{altExp = NULL}. This is the default 
#'   choice.
#'   
#' @param assay.type A single character value for selecting the
#'   \code{\link[SummarizedExperiment:SummarizedExperiment-class]{assay}}
#'   to use for identifying dominant taxa.
#'
#' @param assay_name a single \code{character} value for specifying which
#'   assay to use for calculation.
#'   (Please use \code{assay.type} instead. At some point \code{assay_name}
#'   will be disabled.)
#'  
#'  @param MARGIN A character parameter for selecting if the clustering is done
#'   row-wise / for features ("features") or column-wise / for samples 
#'   ("samples"). Must be \code{1} or  \code{2}. 
#'   (By default: \code{MARGIN = "features"}) 
#'   
#' @param full A logical scalar indicating whether the clustering statistics 
#'   should be put in metadata.
#'   
#' @param name A single character value indicating the name of the clustering 
#'   info in the metadata.
#'   
#' @details
#' This is a wrapper for the 
#' \link[https://bioconductor.org/packages/release/bioc/html/bluster.html]{bluster} 
#' package.
#' 
#' It returns a \code{SummarizedExperiment} object with clustering information 
#' named \code{clusters} stored in \code{colData} or \code{rowData}. 
#' 
#' When setting \code{full=TRUE}, the clustering information will be stored in
#' the metadata of the object.
#' 
#' By default, clustering is done on the features.
#'
#' @author Basil Courbayre
#'
#' @examples
#' data(GlobalPatterns)
#' x <- GlobalPatterns
#'
#' # Cluster on rows using Kmeans
#' tse <- cluster(tse, KmeansParam(centers = 3))
#' 
#' # Clustering done on the samples using Hclust
#' tse <- cluster(tse, 
#'                MARGIN="samples", 
#'                HclustParam(metric = "bray", dist.fun = vegan::vegdist)
#' 
#' # Getting the clusters
#' colData(tse)$clusters
#' 
NULL

#' @rdname cluster
#' @export
setGeneric("cluster", signature = c("x"),
           function(x, BLUSPARAM, altExp=NULL,
                    assay.type = assay_name, assay_name = "counts", 
                    MARGIN="features", full = FALSE, name="clusters")
               standardGeneric("cluster"))

#' @rdname cluster
#' @export
#' @importFrom bluster clusterRows
setMethod("cluster", signature = c(x = "SummarizedExperiment"),
          function(x, BLUSPARAM, altExp=NULL,  
                   assay.type = assay_name, assay_name = "counts", 
                   MARGIN="features", full = FALSE, name="clusters",) {
        # Checking parameters
        .check_margin(MARGIN)
        se <- .check_and_get_altExp(x, altExp)
        if (full) {
            .check_name(se, name)
        }
        .check_assay_present(assay.type, se)
        assay <- assay(se, assay.type)
        
        # Transpose if clustering on the columns
        if( MARGIN != "features" ){
            assay <- t(assay)
        }
        
        result <- clusterRows(assay, BLUSPARAM, full)
        
        # Getting the clusters and setting metadata
        if (full) {
            clusters <- result$clusters
            metadata(se)[["clusters"]] <- result$objects
        } else {
            clusters <- result
        }
        
        # Setting clusters in the object
        if (MARGIN == "features") {
            rowData(se)$clusters <- clusters
        } else {
            colData(se)$clusters <- clusters
        }
        
        if (!is.null(altExp) && any(dim(se) != dim(x))) {
            altExp(x, altExp) <- se
        } else {
            x <- se
        }
        
        x
    }
)

.check_margin <- function(MARGIN) {
    if( !.is_non_empty_string(MARGIN) && !MARGIN %in% c("features", "samples") ){
        stop("'MARGIN' must be 'features' or 'samples'.", call. = FALSE)
    }
}

.check_name <- function(x, name) {
    if (name %in% names(metadata(x))) {
        stop("'name' must be unexisting in the metadata of the object.", call. = FALSE)
    }
}