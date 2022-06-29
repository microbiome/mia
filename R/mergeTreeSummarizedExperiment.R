#' Merge TreeSE objects into single TreeSE object.
#' 
#' @param x a \code{\link{SummarizedExperiment}} object or a list of 
#' \code{\link{SummarizedExperiment}} objects.
#' 
#' @param y a \code{\link{SummarizedExperiment}} object 
#' 
#' @param abund_values A single character value for selecting the
#' \code{\link[SummarizedExperiment:SummarizedExperiment-class]{assay}}
#' to be merged.
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
        function(x, abund_values = "counts", missing_values = 0, ..., BPPARAM = SerialParam() )
            standardGeneric("mergeTreeSummarizedExperiment"))

#' @rdname mergeTreeSummarizedExperiment
#' @export
setMethod("mergeTreeSummarizedExperiment", signature = c(x = "list"),
        function(x, abund_values = "counts", missing_values = 0, ..., BPPARAM = SerialParam() ){
            ################## Input check ##################
            # Can the abund_value tbe found form all the objects
            abund_values_bool <- lapply(x, .assay_cannot_be_found, abund_values = abund_values)
            if( any(abund_values_bool) ){
                stop("'abund_values' must specify an assay from assays. 'abund_values' ",
                     "cannot be found at least in one TreeSE.",
                     call. = FALSE)
            }
            # Is missing_values one of the allowed ones
            missing_values_bool <- (is.numeric(missing_values) && missing_values = 0) ||
                .is_a_string(missing_values) || is.na(missing_values)
            # If not then give error
            if(  missing_values_bool ){
                stop("'missing_values' must be 0, NA, or a single character value.",
                     call. = FALSE)
            }
            ################ Input check end ################
            # Take first element and remove it from the list
            tse <- x[[1]]
            x[[1]] <- NULL
            # Remove all information but rowData, colData, and assay
            rowdata <- rowData(tse)
            colData(tse) <- colData(tse)
            assay <- assay(tse, abund_values)
            assays <- SimpleList(name = assay)
            names(assays) <- abund_values
            tse <- TreeSummarizedExperiment(assays = assays,
                                            rowData = rowdata,
                                            colData = colData
                                            )
            
            for(i in list ){
                tse <- .merge_TreeSE(i, tse_original = tse, missing_values = missing_values)
            }
              
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

################################ HELP FUNCTIONS ################################

.assay_cannot_be_found <- function(abund_values, obj){
    tryCatch(
        {
            .check_assay_present(abund_values, obj)
            return(FALSE)
            
        },
        error=function(cond) {
            
            return(TRUE)
        }
    )
}

.merge_TreeSE <- function(tse_original, tse, missing_values){
    rowdata <- .merge_rowdata(tse_original, tse)
    coldata <- .merge_coldata(tse_original, tse)
    assays <- .merge_assays(tse_original, tse)
}

.merge_rowdata <- function(tse_original, tse){
    
}

.merge_coldata <- function(tse_original, tse){
    
}

.merge_assays <- function(tse_original, tse){
    
}

