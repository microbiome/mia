#' Clusters rowData using \code{fastcluster::hclust}
#'
#' @details
#' \code{rowClust} Performs hierarchical clustering on rowData of the 
#' \code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}, 
#' by using the dissimilarity function provided at \code{dissimilarity_FUN}.
#'
#' @param x a
#'   \code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#'   object.
#'   
#' @param abund_values name of the assay containing the data.
#' 
#' @param name name of the new column at rowData where to place cluster ids.
#' 
#' @param clustering_method clustering method to be used, same as methods 
#' available at \code{\link[hclust]{hclust}}.
#' 
#' @param transposed whether the data should be transposed.
#' 
#' @param dissimilarity_FUN function to compute the dissimilarity matrix
#' 
#' @param k number of clusters (same as in \code{\link[hclust]{hclust}}).
#' 
#' @param h cut-off height (same as \code{\link[hclust]{hclust}})
#'
#' @return \code{x} with additional \code{\link{rowData}} named \code{*name*}
#'
#' 
#' @name rowClust
#' 
#' @export
#'
#' @examples
#' # Load example data
#' library(mia)
#' data("GlobalPatterns")
#' tse <- GlobalPatterns
#' # Agglomerate to Class rank
#' tse <- agglomerateByRank(tse, rank="Class")
#' # Cluster using spearman correlation as a dissimilarity and ward.D2 
#' # for clustering; with results having 6 clusters
#' tse <- rowClust(tse, k = 6)
#' # Using cosine dissimilarity, average clustering and 0.6 as the height
#' # of cut-off
#' tse <- rowClust(tse, 
#'     dissimilarity_FUN=function(x) 1 - coop::cosine(x),
#'     clustering_method = "average",
#'     h = 0.6)
#' # Aggregate clusters as a sum of each cluster values
#' # (CAG: Co-Abundance Groups)
#' tse_CAG <- mergeRows(tse, rowData(tse)$clust_id)
#'
NULL

#' @importFrom fastcluster hclust
#' @importFrom stats cutree as.dist

#' @rdname rowClust
#' @export
setGeneric("rowClust",signature = c("x"),
    function(x, abund_values = "counts", name = "clust_id",
             dissimilarity_FUN = function(x) stats::cor(x, method = "spearman"),
             clustering_method = "ward.D2", transposed = TRUE, k = NULL,
             h = NULL)
        standardGeneric("rowClust"))

#' @rdname rowClust
#' @export
setMethod("rowClust", signature = c(x="SummarizedExperiment"), 
    function(x, abund_values = "counts", name = "clust_id",
             dissimilarity_FUN = function(x) stats::cor(x, method = "spearman"),
             clustering_method = "ward.D2", transposed = TRUE, k = NULL,
             h = NULL) {

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
        method <- match.arg(clustering_method,
                            c("single", "complete", "average", "mcquitty",
                              "ward.D", "ward.D2", "centroid","median"))
        hc <- hclust(as.dist(mat), method = method)
        
        # Checking and setting cut-off
        .check_cutoffs(k,h)
        clusters <- cutree(hc, k=k, h=h)
        rowData(x)[,name] <- as.factor(clusters)
        return(x)
    }
)

############################## HELPER FUNCTIONS ##############################

.check_cutoffs <- function(k,h) {
    if(is.null(k) && is.null(h))
        stop("Specify a cut-off point; k(number of clusters) or h(height)",
             call. = FALSE)
    else if (!is.numeric(c(k,h)) ||
             (!is.null(k) && length(k)!=1L) ||
             (!is.null(h) && length(h)!=1L))
        stop("k or h must be a single numeric value",
             call. = FALSE)
}