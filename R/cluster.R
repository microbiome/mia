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
    function(
            x, BLUSPARAM, assay.type = assay_name, 
            assay_name = "counts", MARGIN = "features", full = FALSE, 
            name = "clusters", clust.col = "clusters", ...)
    standardGeneric("cluster"))

#' @rdname cluster
#' @export
#' @importFrom bluster clusterRows
setMethod("cluster", signature = c(x = "SummarizedExperiment"),
    function(
            x, BLUSPARAM, assay.type = assay_name, 
            assay_name = "counts", MARGIN = "features", full = FALSE, 
            name = "clusters", clust.col = "clusters", ...) {
        .require_package("bluster")
        # Checking parameters
        MARGIN <- .check_MARGIN(MARGIN)
        se <- .check_and_get_altExp(x, ...)
        .check_assay_present(assay.type, se)
        if( !.is_a_string(name) ){
            stop("'name' must be a non-empty single character value.",
                call. = FALSE)
        }
        if( !.is_a_string(clust.col) ){
            stop("'clust.col' must be a non-empty single character value.",
                call. = FALSE)
        }
        if( !.is_a_bool(full) ){
            stop("'full' must be TRUE or FALSE.", call. = FALSE)
        }
        #
        # Get assay
        mat <- assay(se, assay.type)
        # Transpose if clustering on the columns
        if(MARGIN == 2){
            mat <- t(mat)
        }
        # Get clusters
        result <- clusterRows(mat, BLUSPARAM, full)
        # If user has specified full=TRUE, result includes additional info
        # that will be stored to metadata.
        if (full) {
            clusters <- result$clusters
            x <- .add_values_to_metadata(x, name, result$objects, ...)
        } else {
            clusters <- result
        }
        # Setting clusters in the object. The adding function requires data as
        # list
        clusters <- list(clusters)
        x <- .add_values_to_colData(
            x, clusters, clust.col, MARGIN = MARGIN, colname = "clust.col", ...)
        return(x)
    }
)
