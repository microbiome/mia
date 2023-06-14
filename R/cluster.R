#' Clustering wrapper
#'
#' This function returns a \code{SummarizedExperiment} with clustering 
#'   information in its colData or rowData
#'
#' @param x A
#'   \code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#'   object.
#'   
#' @param data.name A single character value indicating the name of the 
#'   \code{rowData} (or \code{colData}) where the data will be stored.
#'   
#' @param ... Additional parameters to use altExps for example
#' @inheritParams bluster::clusterRows
#' @inheritParams runDMN
#' @inheritParams transformCounts
#'   
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
#' library(bluster)
#' data(GlobalPatterns, package = "mia")
#' tse <- GlobalPatterns
#'
#' # Cluster on rows using Kmeans
#' tse <- cluster(tse, KmeansParam(centers = 3))
#' 
#' # Clustering done on the samples using Hclust
#' tse <- cluster(tse, 
#'                MARGIN="samples", 
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
                    assay_name = "counts", MARGIN="features", full = FALSE, 
                    name="clusters", data.name="clusters", ...)
               standardGeneric("cluster"))

#' @rdname cluster
#' @export
#' @importFrom bluster clusterRows
setMethod("cluster", signature = c(x = "SummarizedExperiment"),
          function(x, BLUSPARAM, assay.type = assay_name, 
                   assay_name = "counts", MARGIN="features", full = FALSE, 
                   name="clusters", data.name="clusters", ...) {
        # Checking parameters
        MARGIN <- .check_margin(MARGIN)
        altexp <- .get_altExp(...)
        se <- .check_and_get_altExp(x, altexp)
        .check_data_name(se, data.name, MARGIN)
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
            rowData(se)[[data.name]] <- clusters
        } else {
            colData(se)[[data.name]] <- clusters
        }
        
        if (!is.null(altexp) && any(dim(se) != dim(x))) {
            altExp(x, altexp) <- se
        } else {
            x <- se
        }
        
        x
    }
)

.get_altExp <- function(...) {
    altexppos <- which(...names() == "altexp")
    if (length(altexppos) == 0) {
        NULL
    } else {
        ...elt(altexppos)
    }
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
}

.check_name <- function(x, name) {
    if (name %in% names(metadata(x))) {
        stop("The 'name' must not exist in the metadata of the object.", call. = FALSE)
    }
}

.check_data_name <- function(x, data.name, MARGIN) {
    if (MARGIN == 1) {
        if (data.name %in% names(rowData(x))) {
            stop("The 'data.name' parameter must not exist in the names of the rowData of the object.", 
                 call. = FALSE)
        }
    } else {
        if (data.name %in% names(colData(x))) {
            stop("The 'data.name' parameter must not exist in the names of the colData of the object.", 
                 call. = FALSE)
        }
    }
}