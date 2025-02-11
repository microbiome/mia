#' @name
#' utilization_functions
#'
#' @title
#' Utilization functions for \code{TreeSummarizedExperiment} objects
#'
#' @description
#' A set of utility functions designed to facilitate operations with
#' \code{TreeSummarizedExperiment} objects
#'
#' @details
#' \code{getReducedDimElement} is a utility function that retrieves specific
#' elements from the attributes of \code{reducedDim} in a
#' \code{TreeSummarizedExperiment} object. These attributes may contain
#' loadings, statistical test results, or other metadata, depending on the
#' methods used to generate the results.
#'
#' @return
#' The extracted element from the \code{reducedDim} attribute.
#'
#' @inheritParams addAlpha
#'
#' @param dimred \code{Character scalar} or \code{integer scalar}. A name or
#' index of dimension reduction results. (Default: \code{1L})
#'
#' @param name \code{Character vector}. A name of values retrieved from
#' attributes of \code{reducedDim(x, dimred)}. If \code{NULL}, all the values
#' are retrieved. (Default: \code{NULL})
#'
#' @param ... additional arguments, not used currently.
#'
#' @examples
#' data(GlobalPatterns)
#' tse <- GlobalPatterns
#'
#' # Reduce the number of features
#' tse <- agglomerateByPrevalence(tse, rank = "Phylum")
#'
#' # Run NMF and add the result to reducedDim(tse, "NMF").
#' tse <- addNMF(tse, k = 1, name = "NMF")
#'
#' # Extract feature loadings
#' res <- getReducedDimElement(tse, dimred = "NMF", name = "loadings")
#' res |> head()
#'
#' @seealso
#' \code{\link[=runCCA]{runCCA}}, \code{\link[=addNMF]{addNMF}}, and
#' \code{\link[=addLDA]{addLDA}}
#'
NULL

#'
#' @export
#' @rdname utilization_functions
#' @importFrom SingleCellExperiment reducedDim
setMethod("getReducedDimElement", "SingleCellExperiment",
    function(x, dimred = 1L, name = NULL, ...){
        if( !(is.null(name) || is.character(name) || .is_integer(name)) ){
            stop("'name' must be NULL, character or integer value.",
                call. = FALSE)
        }
        # Check and get attributes from reducedDim
        temp <- .check_dimred_present(dimred, x)
        mat <- reducedDim(x, dimred)
        values <- attributes(mat)
        # Extract the user-specified elements
        values <- .extract_elements_from_reduceddim(values, name)
        return(values)
    }
)

################################ HELP FUNCTIONS ################################

# From a list, extract those elements that user has specified
.extract_elements_from_reduceddim <- function(res, name){
    # Remove matrix-specific attributes
    rm <- c("dim", "dimnames")
    res <- res[ !names(res) %in% rm ]
    # Check that name is correct
    if( is.character(name) && !any(names(res) %in% name) ){
        stop("'name' must be from the following following options: '",
            paste0(names(res), collapse = "', '"), "'", call. = FALSE)
    }
    if( .is_integer(name) && !all(name>0L & name<=length(res)) ){
        stop("'name' must be from the following range (0, ", length(res), "]",
            call. = FALSE)
    }
    # Get values specified by name or integer
    if( is.character(name) || .is_integer(name) ){
        res <- res[ name ]
    }
    # If there was only single match (which should be in the most cases)
    # extract the element from the list
    if( length(res) == 1L ){
        res <- res[[1L]]
    }
    return(res)
}
