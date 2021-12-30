#' Cluster assay rows
#'
#' \code{rowClust} Performs clustering on assay rows of the 
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
#' @param clustering_FUN a \code{stats::hclust} a function defined and provided by
#'  user, to compute the clustering. Note! The function should return an \code{hclust}
#'  object. Example functions: 'stats::hclust' or 'fastcluster::hclust'
#' 
#' @param dissimilarity_FUN \code{function(x) 1 - stats::cor(t(x))} a function defined and provided by
#'  user, to compute the dissimilarity matrix. Note! the function expects an input
#'  of type (features x samples).
#' 
#' @param ... other arguments passed onto \code{\link[stats:cutree]{cutree}}
#'
#' @return \code{x} with additional column at \code{rowData} named \code{*name*}
#'
#' @export
#'
#' @references
#' Claesson, M., Jeffery, I., Conde, S. et al. Gut microbiota composition
#'  correlates with diet and health in the elderly. Nature 488, 178–184 (2012)
#' \url{https://doi.org/10.1038/nature11319}
#'
#' Minot, S. S., & Willis, A. D. (2019). Clustering co-abundant genes identifies
#' components of the gut microbiome that are reproducibly associated with 
#' colorectal cancer and inflammatory bowel disease. Microbiome, 7(1), 110.
#' \url{https://doi.org/10.1186/s40168-019-0722-6}
#'
#' @name rowClust
#' @export
#' 
#' @examples
#' # Load example data
#' data("GlobalPatterns")
#' tse <- GlobalPatterns
#' # Agglomerate to Class rank; for speeding up the example.
#' tse <- agglomerateByRank(tse, rank="Class")
#' # Transforming abundance data with clr
#' tse <- transformCounts(tse, "counts", "clr", pseudocount = 1)
#' # Clustering with default methods using clr transform and 6 clusters as output
#' tse <- rowClust(tse, abund_values = "clr", k = 6)
#' # Spearman correlation is used as a dissimilarity and ward.D2 
#' # for clustering; with results having 6 clusters:
#' tse <- rowClust(tse,
#'                 abund_values = "clr",
#'                 clustering_FUN = function(x) stats::hclust(x, method="ward.D2"),
#'                 dissimilarity_FUN = function(x) 1 - stats::cor(t(x),method="spearman"),
#'                 k = 6)
#' # OR, using fastcluster package
#' \dontrun{
#' tse <- rowClust(tse,
#'                 abund_values = "clr",
#'                 clustering_FUN = function(x) fastcluster::hclust(x, method="ward.D2"),
#'                 dissimilarity_FUN = function(x) 1 - stats::cor(t(x),method="spearman"),
#'                 k = 6)
#' }
#' # Removing clr assay, since the next aggregation of values would be 
#' # meaningless on clr transform.
#' assay(tse, "clr") <- NULL
#' # Aggregate clusters as a sum of each cluster values
#' # (CAG: Co-Abundance Groups)
#' tse_CAG <- mergeRows(tse, rowData(tse)$clust_id)
#'
NULL



#' @rdname rowClust
#' @export
setGeneric("rowClust",signature = c("x"),
    function(x, abund_values = "counts", name = "clust_id",
             dissimilarity_FUN = function(x) 1 - stats::cor(t(x)),
             clustering_FUN = stats::hclust, ...)
        standardGeneric("rowClust"))

#' @rdname rowClust
#' @importFrom stats hclust cutree as.dist
#' @export
setMethod("rowClust", signature = c(x="SummarizedExperiment"), 
    function(x, abund_values = "counts", name = "clust_id",
             dissimilarity_FUN = function(x) 1 - stats::cor(t(x)),
             clustering_FUN = stats::hclust, ...) {

        # Checking and getting data
        .check_assay_present(abund_values, x)
        mat <- assay(x, abund_values)
        
        # Calculating dissimilarity
        if(!.is_function(dissimilarity_FUN))
            stop("'dissimilarity_FUN' must be a function.",
                 call. = FALSE)
        mat <- do.call(dissimilarity_FUN, list(mat))
        
        # Clustering
        if (!.is_function(clustering_FUN))
            stop("'", clustering_FUN, "' is not a function", call. = FALSE)
        hc <- do.call(clustering_FUN, list(as.dist(mat)))
        .check_clust_obj(hc)
        
        # Checking and setting cut-off
        clusters <- cutree(hc, ...)
        rowData(x)[,name] <- as.factor(clusters)
        return(x)
    }
)

################################ HELPER FUNCTIONS ##############################

.check_clust_obj <- function(obj){
    if (!is(obj, "hclust"))
        stop("The clustering function should return an object of class 'hclust'.
        Example functions to be used: 'stats::hclust' or 'fastcluster::hclust'",
             call. = FALSE)
}