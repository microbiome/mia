#' Clustering wrapper
#'
#' This function returns a \code{SummarizedExperiment} with clustering 
#'   information in its colData or rowData
#'
#' @param x A
#'   \code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#'   object.
#'   
#' @param clust.col A single character value indicating the name of the 
#'   \code{rowData} (or \code{colData}) where the data will be stored.
#'   
#' @param ... Additional parameters to use altExps for example
#' @inheritParams bluster::clusterRows
#' @inheritParams runDMN
#' @inheritParams transformAssay
#'   
#'   
#' @details
#' This is a wrapper for the \code{clusterRows} function from the 
#' \link[https://bioconductor.org/packages/release/bioc/html/bluster.html]{bluster} 
#' package.
#'
#' When setting \code{full = TRUE}, the clustering information will be stored in
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
#' library(bluster)
#' data(GlobalPatterns, package = "mia")
#' tse <- GlobalPatterns
#'
#' # Cluster on rows using Kmeans
#' tse <- cluster(tse, KmeansParam(centers = 3))
#' 
#' # Clustering done on the samples using Hclust
#' tse <- cluster(tse, 
#'                MARGIN = "samples", 
#'                HclustParam(metric = "bray", dist.fun = vegan::vegdist))
#' 
#' # Getting the clusters
#' colData(tse)$clusters
#' 
NULL

#' @rdname cluster
#' @export
setGeneric("cluster", signature = c("x"),
           function(x, BLUSPARAM, assay.type = assay_name, 
                    assay_name = "counts", MARGIN = "features", full = FALSE, 
                    name = "clusters", clust.col = "clusters", ...)
               standardGeneric("cluster"))

#' @rdname cluster
#' @export
#' @importFrom bluster clusterRows
setMethod("cluster", signature = c(x = "SummarizedExperiment"),
          function(x, BLUSPARAM, assay.type = assay_name, 
                   assay_name = "counts", MARGIN = "features", full = FALSE, 
                   name = "clusters", clust.col = "clusters", ...) {
        # Checking parameters
        MARGIN <- .check_MARGIN(MARGIN)
        se <- .check_and_get_altExp(x, ...)
        .check_assay_present(assay.type, se)
        if( !.is_a_string(name) ){
            stop("'name' must be a non-empty single character value.",
                 call. = FALSE)
        }
        #
        mat <- assay(se, assay.type)
        
        # Transpose if clustering on the columns
        if(MARGIN == 2){
            mat <- t(mat)
        }
        
        result <- clusterRows(mat, BLUSPARAM, full)
        # Getting the clusters and setting metadata
        if (full) {
            clusters <- result$clusters
            se <- .add_values_to_metadata(se, name, result$objects)
        } else {
            clusters <- result
        }
        clusters <- list(clusters)
        # Setting clusters in the object
        se <- .add_values_to_colData(se, clusters, name, MARGIN = 1, ...)
        return(se)
    }
)
