#' @export
setGeneric("subsetSamples", signature = "x",
           function(x, ...)
               standardGeneric("subsetSamples"))
#' @export
setGeneric("subsetFeatures", signature = "x",
           function(x, ...)
               standardGeneric("subsetFeatures"))
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

#' @export
setMethod("subsetSamples", signature = "SummarizedExperiment",
          function(x, ...){
              .Deprecated(msg = paste0("subsetSamples is deprecated. Please ",
                                       "use '[]' for subsetting instead."))
              subset_args <- .get_subset_args(x, ...)
              x[subset_args$columns,subset_args$rows]
          }
)

#' @export
setMethod("subsetFeatures", signature = "SummarizedExperiment",
          function(x, ...){
              .Deprecated(msg = paste0("subsetFeatures is deprecated. Please",
                                       " use '[]' for subsetting instead."))
              subset_args <- .get_subset_args(x, ...)
              x[subset_args$rows, subset_args$columns]
          }
)

#' @export
setMethod("subsetTaxa", signature = "SummarizedExperiment",
          function(x, ...){
              .Deprecated(msg = paste0("subsetFeatures is deprecated. Please",
                                       " use '[]' for subsetting instead."))
              subsetFeatures(x, ...)
          }
)
