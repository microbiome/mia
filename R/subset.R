#' Subset functions
#'
#' To make a transition from \code{phyloseq} easier, the `subsetSamples` and
#' `subsetFeatures` functions are implemented. To avoid name clashes they are
#' named differently.
#'
#' However, the use of these functions is discouraged since subsetting using
#' \code{\link[SummarizedExperiment:SummarizedExperiment-class]{[}} works on
#' both dimension at the same time, is more flexible and is used hroughout R to
#' subset data with two or more dimension. Therefore, these functions will be
#' removed in Bioconductor release 3.15 (April, 2022).
#'
#' @param x a
#'   \code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#'   object
#'
#' @param ... See \code{\link[BiocGenerics:subset]{subset}}. \code{drop} is
#'   not supported.
#'
#' @name subsetSamples
#'
#' @return
#' A subset of \code{x}
#'
#' @examples
#' data(GlobalPatterns)
#' subsetSamples(GlobalPatterns, colData(GlobalPatterns)$SampleType == "Soil")
#' subsetFeatures(GlobalPatterns, rowData(GlobalPatterns)$Kingdom == "Bacteria")
NULL

#' @rdname subsetSamples
#' @export
setGeneric("subsetSamples", signature = "x",
           function(x, ...)
               standardGeneric("subsetSamples"))
#' @rdname subsetSamples
#' @export
setGeneric("subsetFeatures", signature = "x",
           function(x, ...)
               standardGeneric("subsetFeatures"))
#' @rdname subsetSamples
#' @export
setGeneric("subsetTaxa", signature = "x",
           function(x, ...)
               standardGeneric("subsetTaxa"))

.get_subset_args <- function(x, subset = NULL, select = NULL, ...){
    rows <- subset
    columns <- select
    if(is.null(rows)){
        rows <- rep(TRUE,nrow(x))
    }
    if(is.null(columns)){
        columns <- rep(TRUE,ncol(x))
    }
    return(list(rows = rows, columns = columns))
}

#' @rdname subsetSamples
#' @export
setMethod("subsetSamples", signature = "SummarizedExperiment",
    function(x, ...){
        message("Please use '[]' for subsetting instead.", call. = FALSE)
        subset_args <- .get_subset_args(x, ...)
        x[subset_args$columns,subset_args$rows]
    }
)

#' @rdname subsetSamples
#' @export
setMethod("subsetFeatures", signature = "SummarizedExperiment",
    function(x, ...){
        message("Please use '[]' for subsetting instead.", call. = FALSE)
        subset_args <- .get_subset_args(x, ...)
        x[subset_args$rows, subset_args$columns]
    }
)

#' @rdname subsetSamples
#' @export
setMethod("subsetTaxa", signature = "SummarizedExperiment",
    function(x, ...){
        subsetFeatures(x, ...)
    }
)
