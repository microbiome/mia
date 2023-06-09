#' Clustering wrapper
#'
#' This function returns a \code{SummarizedExperiment} with clustering 
#'   information in its colData or rowData
#'
#' @param x A
#'   \code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#'   object.
#' @inheritParams bluster::clusterRows
#' 
#' @inheritParams runDMN
#' @inheritParams transformCounts
#'   
#' @param altexp A single character value indicating which alternate
#'   experiment to use for the clustering (if any). If the clustering is done
#'   on the main experiment, set \code{altexp = NULL}. This is the default 
#'   choice.
#'   
#' @details
#' This is a wrapper for the \code{clusterRows} function from the 
#' \link[https://bioconductor.org/packages/release/bioc/html/bluster.html]{bluster} 
#' package.
#'
#' When setting \code{full=TRUE}, the clustering information will be stored in
#' the metadata of the object.
#' 
#' By default, clustering is done on the features.
#'
#' @return
#' \code{cluster} returns an object of the same type as the \code{x} parameter 
#' with clustering information named \code{clusters} stored in \code{colData} 
#' or \code{rowData}. 
#'
#' @name cluster
#' @export
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
           function(x, BLUSPARAM, altexp=NULL,
                    assay.type = assay_name, assay_name = "counts", 
                    MARGIN="features", full = FALSE, name="clusters")
               standardGeneric("cluster"))

#' @rdname cluster
#' @export
#' @importFrom bluster clusterRows
setMethod("cluster", signature = c(x = "SummarizedExperiment"),
          function(x, BLUSPARAM, altexp=NULL,  
                   assay.type = assay_name, assay_name = "counts", 
                   MARGIN="features", full = FALSE, name="clusters") {
        # Checking parameters
        .check_margin(MARGIN)
        se <- .check_and_get_altExp(x, altexp)
        if (full) {
            .check_name(se, name)
        }
        .check_assay_present(assay.type, se)
        assay <- assay(se, assay.type)
        
        # Transpose if clustering on the columns
        if(MARGIN %in% c("cols", "samples")){
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
        if (MARGIN %in% c("rows", "features")) {
            rowData(se)$clusters <- clusters
        } else {
            colData(se)$clusters <- clusters
        }
        
        if (!is.null(altexp) && any(dim(se) != dim(x))) {
            altExp(x, altexp) <- se
        } else {
            x <- se
        }
        
        x
    }
)

.check_margin <- function(MARGIN) {
    if( !.is_non_empty_string(MARGIN) 
        && !MARGIN %in% c("features", "samples", "rows", "cols") ){
        stop("'MARGIN' must be either 'features', 'samples', 'rows', or 'cols'.",
             call. = FALSE)
    }
}

.check_name <- function(x, name) {
    if (name %in% names(metadata(x))) {
        stop("The 'name' must not exist in the metadata of the object.", call. = FALSE)
    }
}