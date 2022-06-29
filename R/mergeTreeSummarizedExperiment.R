#' Merge TreeSE objects into single TreeSE object.
#' 
#' @param x a \code{\link{SummarizedExperiment}} object or a list of 
#' \code{\link{SummarizedExperiment}} objects.
#' 
#' @param y a \code{\link{SummarizedExperiment}} object 
#' 
#' @param missing_values NA, 0, or a single character values specifying the notation
#' of missing values.
#'
#'
#' @param BPPARAM A
#'   \code{\link[BiocParallel:BiocParallelParam-class]{BiocParallelParam}}
#'   object specifying whether calculation of estimates should be parallelized.
#'
#' @param ... optional arguments (not used).
#'
#' @return A single \code{TreeSummarizedExperiment} object.
#'
#' @details
#' This function merges multiple \code{SummarizedExperiment} objects. It combines
#' \code{rowData}, \code{assays}, and \code{colData} so that the output includes
#' each unique row and column ones. If, for example, all rows are not shared with
#' individual objects, there are missing values in \code{assays}. The notation of missing
#' can be specified with the \code{missing_values} argument.
#'
#' @references
#'
#' @seealso
#'
#' @name mergeTreeSummarizedExperiment
#' @export
#'
#' @author Leo Lahti and Tuomas Borman. Contact: \url{microbiome.github.io}
#' 
#' @examples
#' data(GlobalPatterns)
#' tse <- GlobalPatterns
#' 
#' 
NULL

######################### Function for list of TreeSEs #########################

#' @rdname mergeTreeSummarizedExperiment
#' @export
setGeneric("mergeTreeSummarizedExperiment", signature = c("x"),
        function(x, missing_values = 0, ..., BPPARAM = SerialParam() )
            standardGeneric("mergeTreeSummarizedExperiment"))

#' @rdname mergeTreeSummarizedExperiment
#' @export
setMethod("mergeTreeSummarizedExperiment", signature = c(x = "list"),
        function(x, missing_values = 0, ..., BPPARAM = SerialParam() ){
              
        }
)

########################### Function for two TreeSEs ###########################

#' @rdname mergeTreeSummarizedExperiment
#' @export
setGeneric("mergeTreeSummarizedExperiment", signature = c("x", "y"),
        function(x, y, ... )
            standardGeneric("mergeTreeSummarizedExperiment"))

#' @rdname mergeTreeSummarizedExperiment
#' @export
setMethod("mergeTreeSummarizedExperiment", signature = c(x = "SummarizedExperiment", 
                                                         y = "SummarizedExperiment"),
        function(x, y, ...){
            # Create a list based on TreeSEs
            list <- list(x, y)
            # Call the function for list
            mergeTreeSummarizedExperiment(list, ...)
        }
)
