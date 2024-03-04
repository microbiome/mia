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
#' \code{addCluster} returns an object of the same type as the \code{x} parameter 
#' with clustering information named \code{clusters} stored in \code{colData} 
#' or \code{rowData}. 
#'
#' @name addCluster
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
#' tse <- addCluster(tse, KmeansParam(centers = 3))
#' 
#' # Clustering done on the samples using Hclust
#' tse <- addCluster(tse, 
#'                MARGIN = "samples", 
#'                HclustParam(metric = "bray", dist.fun = vegan::vegdist))
#' 
#' # Getting the clusters
#' colData(tse)$clusters
#' 
NULL

#' @rdname addCluster
#' @export
setGeneric("addCluster", signature = c("x"),
            function(x, BLUSPARAM, assay.type = assay_name, 
                    assay_name = "counts", MARGIN = "features", full = FALSE, 
                    name = "clusters", clust.col = "clusters", ...)
                standardGeneric("addCluster"))

#' @rdname addCluster
#' @export
#' @importFrom bluster clusterRows
setMethod("addCluster", signature = c(x = "SummarizedExperiment"),
            function(x, BLUSPARAM, assay.type = assay_name, 
                    assay_name = "counts", MARGIN = "features", full = FALSE, 
                    name = "clusters", clust.col = "clusters", ...) {
        # Checking parameters
        MARGIN <- .check_margin(MARGIN)
        se_altexp <- .get_altExp(x, ...)
        se <- se_altexp$se
        # If there wasn't an altExp in the SE (or if the name was wrong), altexp 
        # is set to NULL, if it exists altexp contains the name of the altExp
        altexp <- se_altexp$altexp
        .check_data_name(se, clust.col, MARGIN)
        .check_assay_present(assay.type, se)
        
        if (full) {
            .check_name(se, name)
        }
        mat <- assay(se, assay.type)
        
        # Transpose if clustering on the columns
        if(MARGIN == 2){
            mat <- t(mat)
        }
        
        result <- clusterRows(mat, BLUSPARAM, full)
        # Getting the clusters and setting metadata
        if (full) {
            clusters <- result$clusters
            metadata(se)[[name]] <- result$objects
        } else {
            clusters <- result
        }
        
        # Setting clusters in the object
        if (MARGIN == 1) {
            rowData(se)[[clust.col]] <- clusters
        } else {
            colData(se)[[clust.col]] <- clusters
        }
        
        # If there was an altexp, update it in the mainExp
        if (!is.null(altexp)) {
            altExp(x, altexp) <- se
        } else {
            x <- se
        }
        
        x
    }
)

.get_altExp <- function(x, ...) {
    altexppos <- which(...names() == "altexp")
    if (length(altexppos) == 0) {
        altexp <- NULL
    } else {
        altexp <- ...elt(altexppos)
    }
    se <- .check_and_get_altExp(x, altexp)
    list(se = se, altexp = altexp)
}

.check_margin <- function(MARGIN) {
    if (.is_non_empty_string(MARGIN)) {
        MARGIN <- tolower(MARGIN)
    }
    if (length(MARGIN) != 1L 
        || !(MARGIN %in% c(1, 2, "features", "samples", "columns", 
                            "col", "row", "rows", "cols"))) {
        stop("'MARGIN' must equal to either 1, 2, 'features', 'samples', 'columns', 'col', 'row', 'rows', or 'cols'.",
            call. = FALSE)
    }
    MARGIN <- ifelse(MARGIN %in% c("samples", "columns", "col", 2, "cols"), 
                    2, 1)
    MARGIN
}

.check_name <- function(x, name) {
    if (name %in% names(metadata(x))) {
        stop("The 'name' must not exist in the metadata of the object.", call. = FALSE)
    }
}

.check_data_name <- function(x, clust.col, MARGIN) {
    if (MARGIN == 1) {
        if (clust.col %in% names(rowData(x))) {
            stop("The 'clust.col' parameter must not exist in the names of the rowData of the object.", 
                call. = FALSE)
        }
    } else {
        if (clust.col %in% names(colData(x))) {
            stop("The 'clust.col' parameter must not exist in the names of the colData of the object.", 
                call. = FALSE)
        }
    }
}
