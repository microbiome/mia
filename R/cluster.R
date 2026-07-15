#' Clustering wrapper
#'
#' This function returns a \code{SummarizedExperiment} with clustering
#'   information in its colData or rowData
#'
#' @param x A
#' \code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#' object.
#'
#' @param assay.type \code{Character scalar}. Specifies the name of assay
#' used in calculation.
#'
#' @param dimred \code{Character scalar} or \code{integer scalar}.
#' Specifies dimension reduction results used in calculation. Either
#' \code{dimred} or \code{assay.type} must be specified.
#'
#' @param by \code{Character scalar}. Determines if association is calculated
#'   row-wise / for features ('rows') or column-wise / for samples ('cols').
#'   Must be \code{'rows'} or \code{'cols'}.
#'
#' @param MARGIN Deprecated. Use \code{by} instead.
#'
#' @param clust.col \code{Character scalar}. Indicates the name of the
#'   \code{rowData} (or \code{colData}) where the data will be stored.
#'   (Default: \code{"clusters"})
#'
#' @param full \code{Logical scalar}. Specifies whether additional clustering
#'   information should be stored in the object metadata.
#'
#' @param ... Additional parameters to use altExps for example
#' @param BLUSPARAM A \link[bluster]{BlusterParam-class} object specifying the algorithm to use.
#' @inheritParams runDMN
#' @inheritParams transformAssay
#'
#'
#' @details
#' This is a wrapper for the \code{clusterRows} function from the
#' \link[bluster]{bluster} package.
#'
#' When setting \code{full = TRUE}, the clustering information will be stored in
#' the metadata of the object.
#'
#' By default, clustering is done on the features.
#'
#' @return
#' \code{addCluster} returns an object of the same type as the \code{x}
#' parameter with clustering information named \code{clusters} stored in
#' \code{colData} or \code{rowData}.
#'
#' @name addCluster
#' @export
#'
#' @examples
#' library(bluster)
#' data(GlobalPatterns, package = "mia")
#' tse <- GlobalPatterns
#'
#' # Cluster on rows using Kmeans
#' tse <- addCluster(
#'     tse,
#'     assay.type = "counts",
#'     BLUSPARAM = KmeansParam(centers = 3)
#' )
#'
#' # Clustering done on the samples using Hclust
#' res <- getCluster(
#'     tse,
#'     assay.type = "counts",
#'     by = "samples",
#'     BLUSPARAM = HclustParam(metric = "bray", dist.fun = vegan::vegdist)
#' )
#' res |> head()
#'
#' # Apply clustering to PCA results
#' library(scater)
#' tse <- transformAssay(tse, method = "rclr")
#' tse <- runPCA(tse, assay.type = "rclr")
#' tse <- addCluster(
#'     tse,
#'     dimred = "PCA",
#'     BLUSPARAM = KmeansParam(centers = 3),
#'     by = 2
#' )
#' tse$cluster |> head()
#'
NULL

#' @rdname addCluster
#' @export
setMethod("addCluster", signature = c(x = "SummarizedExperiment"),
    function(
            x, BLUSPARAM, assay.type = assay_name,
            assay_name = NULL, by = MARGIN, MARGIN = "rows",
            name = "clustering", clust.col = "cluster", full = FALSE, ...){
        by <- .check_MARGIN(by)
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
        result <- getCluster(
            x = x, BLUSPARAM = BLUSPARAM, assay.type = assay.type,
            by = by, full = full, ...)
        # If user has specified full=TRUE, result includes additional info
        # that will be stored to metadata.
        if( full ){
            clusters <- result$clusters
        x <- .add_values_to_metadata(x, name, result$objects, ...)
        } else {
            clusters <- result
        }
        # Setting clusters in the object. The adding function requires data as
        # list
        clusters <- list(clusters)
        x <- .add_values_to_colData(
            x, clusters, clust.col, MARGIN = by, colname = "clust.col", ...)
        return(x)
    }
)

#' @rdname addCluster
#' @export
setMethod("addCluster", signature = c(x = "SingleCellExperiment"),
    function(
            x, BLUSPARAM, assay.type = assay_name,
            assay_name = NULL, dimred = NULL, by = MARGIN, MARGIN = "rows",
            name = "clustering", clust.col = "cluster", full = FALSE, ...){
        by <- .check_MARGIN(by)
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
        # Hiddenly support altExp
        tse <- .check_and_get_altExp(x, ...)
        # Calculate indices
        args <- c(list(
            x = tse, BLUSPARAM = BLUSPARAM, assay.type = assay.type,
            dimred = dimred, by = by, full = full), list(...))
        args <- args[ !names(args) %in% c("altexp") ]
        result <- do.call(getCluster, args)
        # If user has specified full=TRUE, result includes additional info
        # that will be stored to metadata.
        if( full ){
            clusters <- result$clusters
            x <- .add_values_to_metadata(x, name, result$objects, ...)
        } else {
            clusters <- result
        }
        # Setting clusters in the object. The adding function requires data as
        # list
        clusters <- list(clusters)
        x <- .add_values_to_colData(
            x, clusters, clust.col, MARGIN = by, colname = "clust.col", ...)
        return(x)
    }
)

#' @rdname addCluster
#' @export
setMethod("getCluster", signature = c(x = "SummarizedExperiment"),
    function(
            x, BLUSPARAM, assay.type = assay_name,
            assay_name = NULL, by = MARGIN, MARGIN = "rows", ...){
        # Checking parameters
        by <- .check_MARGIN(by)
        x <- .check_and_get_altExp(x, ...)
        # Get assay
        .check_assay_present(assay.type, x)
        mat <- assay(x, assay.type)
        #
        # Get clusters
        result <- getCluster(
            x = mat, BLUSPARAM = BLUSPARAM, by = by, ...)
        return(result)
    }
)

#' @rdname addCluster
#' @export
setMethod("getCluster", signature = c(x = "SingleCellExperiment"),
    function(
            x, BLUSPARAM, assay.type = assay_name,
            assay_name = NULL, dimred = NULL, by = MARGIN, MARGIN = "rows",
            ...){
        # Checking parameters
        by <- .check_MARGIN(by)
        x <- .check_and_get_altExp(x, ...)
        # User can specify either abundance matrix or matrix from reducedDim.
        if( sum(c(is.null(assay.type), is.null(dimred))) != 1L ){
            stop("Either 'assay.type' or 'dimred' must be specified.",
                call. = FALSE)
        }
        # Clustering based on PCA or other reduced dim makes only sense, when
        # by is 2.
        if( by == 1 && !is.null(dimred) ){
            stop("Specify 'by=2' when 'dimred' is selected.", call. = FALSE)
        }
        # Get matrix
        if( !is.null(assay.type) ){
            .check_assay_present(assay.type, x)
            mat <- assay(x, assay.type)
        } else{
            .check_dimred_present(dimred, x)
            # Transpose because in reducedDim the samples are in rows while in
            # assay they are in columns.
            mat <- reducedDim(x, dimred) |> t()
        }
        #
        # Get clusters
        result <- getCluster(
            x = mat, BLUSPARAM = BLUSPARAM, by = by, ...)
        return(result)
    }
)

#' @rdname addCluster
#' @export
#' @importFrom bluster clusterRows
setMethod("getCluster", signature = c(x = "ANY"),
    function(x, BLUSPARAM, by = MARGIN, MARGIN = "rows", full = FALSE, ...){
        .require_package("bluster")
        if( !is.matrix(x) ){
            stop("'x' must be a matrix.", call. = FALSE)
        }
        by <- .check_MARGIN(by)
        if( !.is_a_bool(full) ){
            stop("'full' must be TRUE or FALSE.", call. = FALSE)
        }
        #
        # Transpose if clustering on the columns
        if(by == 2){
            x <- x |> t()
        }
        # Get clusters
        result <- clusterRows(x, BLUSPARAM, full)
        return(result)
    }
)
