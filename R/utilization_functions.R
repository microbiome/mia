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
#' @param name \code{Character scalar}. A name of values retrieved from
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
#' res <- getReducedDimElement(tse, "NMF", "loadings")
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
        # Check and get attributes from reducedDim
        temp <- .check_dimred_present(dimred, x)
        mat <- reducedDim(x, dimred)
        values <- attributes(mat)
        # Remove matrix-specific attributes
        rm <- c("dim", "dimnames")
        values <- values[ !names(values) %in% rm ]
        # Check that name is correct
        if( !(is.null(name) || (.is_a_string(name) &&
                any(names(values) %in% name))) ){
            stop("'name' must be NULL or a single character value from the ",
                "following options: '",
                paste0(names(values), collapse = "', '"), "'", call. = FALSE)
        }
        # Get valus specified by name
        if( !is.null(name) ){
            values <- values[ names(values) %in% name ]
        }
        if( length(values) == 1L ){
            values <- values[[1L]]
        }
        return(values)
    }
)
