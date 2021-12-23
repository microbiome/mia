#' Clusters rowData using \code{fastcluster::hclust}
#'
#' \code{rowClust} Performs hierarchical clustering on rowData of the 
#' \code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}, 
#' by using the dissimilarity function provided at \code{dissimilarity_FUN}.
#'
#' @param x a
#'   \code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#'   object.
#'   
#' @param abund_values a \code{character} vector, name of the assay containing the data.
#' 
#' @param name a \code{character} vector, name of the new column at \code{rowData} where to place cluster ids.
#' 
#' @param clustering_method a \code{character} vector representing the clustering 
#' method to be used, same as methods available at \code{\link[hclust]{hclust}}.
#' 
#' @param transposed \code{logical} whether the data should be transposed.
#' 
#' @param dissimilarity_FUN \code{stats::cor} function to compute the dissimilarity matrix
#' 
#' @param ... other arguments passed onto \code{\link[cutree]{cutree}}
#'
#' @return \code{x} with additional \code{\link{rowData}} named \code{*name*}
#'
#' @export
#'
#' @examples
#' # Load example data
#' library(mia)
#' data("GlobalPatterns")
#' tse <- GlobalPatterns
#' # Agglomerate to Class rank; for speeding up the example.
#' tse <- agglomerateByRank(tse, rank="Class")
#' # Cluster using Spearman correlation as a dissimilarity and ward.D2 
#' # for clustering; with results having 6 clusters
#' tse <- rowClust(tse, clustering_method = "ward.D2",
#'                 dissimilarity_FUN = function(x) stats::cor(x, method="spearman"),
#'                 k = 6)
#' # Aggregate clusters as a sum of each cluster values
#' # (CAG: Co-Abundance Groups)
#' tse_CAG <- mergeRows(tse, rowData(tse)$clust_id)
#'
NULL

#' @importFrom fastcluster hclust
#' @importFrom stats cutree as.dist

#' @export
setGeneric("rowClust",signature = c("x"),
    function(x, abund_values = "counts", name = "clust_id",
             dissimilarity_FUN = stats::cor,
             clustering_method = "complete", transposed = TRUE,
             ...)
        standardGeneric("rowClust"))

#' @export
setMethod("rowClust", signature = c(x="SummarizedExperiment"), 
    function(x, abund_values = "counts", name = "clust_id",
             dissimilarity_FUN = stats::cor,
             clustering_method = "complete", transposed = TRUE,
             ...) {

        # Checking and getting data
        .check_assay_present(abund_values, x)
        mat <- assay(x, abund_values)
        if (transposed) mat <- t(mat)
        
        # Calculating dissimilarity
        if(!.is_function(dissimilarity_FUN))
            stop("'dissimilarity_FUN' must be a function.",
                 call. = FALSE)
        mat <- do.call(dissimilarity_FUN, list(mat))
        
        # Clustering
        hc <- hclust(as.dist(mat), method = clustering_method)
        
        # Checking and setting cut-off
        clusters <- cutree(hc, ...)
        rowData(x)[,name] <- as.factor(clusters)
        return(x)
    }
)