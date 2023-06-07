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
#' @param assay.type A single character value for selecting the
#'   \code{\link[SummarizedExperiment:SummarizedExperiment-class]{assay}}
#'   to use for identifying dominant taxa.
#'
#' @param assay_name a single \code{character} value for specifying which
#'   assay to use for calculation.
#'   (Please use \code{assay.type} instead. At some point \code{assay_name}
#'   will be disabled.)
#'  
#'  @param MARGIN A single numeric value for selecting if the clustering is done
#'   row-wise / for features (1) or column-wise / for samples (2). Must be \code{1} or
#'   \code{2}. (By default: \code{MARGIN = 1}) 
#'   
#' @param full A logical scalar indicating whether the clustering statistics 
#'   from both steps should be put in metadata
#'   
#' @details
#' When setting \code{full=TRUE}, the clustering information will be stored in
#' the metadata of the object.
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
NULL

#' @rdname cluster
#' @export
setGeneric("cluster", signature = c("x"),
           function(x, BLUSPARAM, assay.type = assay_name, assay_name = "counts", 
                    MARGIN=1, full = FALSE)
               standardGeneric("cluster"))

#' @rdname cluster
#' @export
#' @importFrom bluster clusterRows
setMethod("cluster", signature = c(x = "SummarizedExperiment"),
    function(x, BLUSPARAM, assay.type = assay_name, assay_name = "counts", 
             MARGIN=1, full = FALSE) {
        # Checking parameters
        .check_assay_present(assay.type, x)
        .check_margin(MARGIN)
        
        assay <- assay(x, assay.type)
        
        # Transpose if clustering on the columns
        if( MARGIN != 1 ){
            assay <- t(assay)
        }
        
        result <- clusterRows(assay, BLUSPARAM, full)
        
        if (full) {
            clusters <- result$clusters
            metadata(x)[["clusters"]] <- result$objects
        } else {
            clusters <- result
        }
        
        if (MARGIN == 1) {
            rowData(x)$clusters <- clusters
        } else {
            colData(x)$clusters <- clusters
        }
        
        x
    }
)

.check_margin <- function(MARGIN) {
    if( !is.numeric(MARGIN) && !MARGIN %in% c(1, 2) ){
        stop("'MARGIN' must be 1 or 2.", call. = FALSE)
    }
}