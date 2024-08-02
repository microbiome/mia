#' Non-negative Matrix Factorization
#'
#' These functions perform Non-negative Matrix Factorization on data stored in a 
#'  \code{\link[TreeSummarizedExperiment:TreeSummarizedExperiment-class]{TreeSummarizedExperiment}}
#'  object.
#' 
#' @param x a \code{\link[TreeSummarizedExperiment:TreeSummarizedExperiment-class]{TreeSummarizedExperiment}}
#'  object.
#'   
#' @param k \code{Integer vector}. A number of latent vectors/topics. 
#'  (Default: \code{2})
#' 
#' @param name \code{Character scalar}. The name to be used to store the result 
#'  in the reducedDims of the output. (Default: \code{"NMF"})
#'  
#' @param assay.type \code{Character scalar}. Specifies which assay to use for 
#'  NMF ordination. (Default: \code{"counts"})
#' 
#' @param ... optional arguments passed to \code{nmf::NMF}
#' 
#' @return 
#' For \code{getNMF}, the ordination matrix with feature loadings matrix
#'  as attribute \code{"loadings"}.
#'  
#' For \code{addNMF}, a \code{\link[TreeSummarizedExperiment:TreeSummarizedExperiment-class]{TreeSummarizedExperiment}}
#'  object is returned containing the ordination matrix in reducedDims(..., name)
#'  with feature loadings matrix as attribute \code{"loadings"}.
#'  
#' @details 
#' The functions \code{getNMF} and \code{addNMF} internally use \code{nmf::NMF} 
#'  compute the ordination matrix and 
#'  feature loadings.
#'  
#' All NMF scores for rank values from 2 to 10 are calculated and the rank value
#'  with highest explained variance is selected. The scores for this rank value 
#'  are returned with the NMF model and loadings as attributes.
#'  
#' @name addNMF
#' 
#' @examples
#' data(GlobalPatterns)
#' tse <- GlobalPatterns
#' 
#' # Reduce the number of features
#' tse <- agglomerateByPrevalence(tse, rank = "Phylum")
#' 
#' # Run NMF and add the result to reducedDim(tse, "NMF")
#' tse <- addNMF(tse)
#' 
#' # Extract feature loadings
#' loadings <- attr(reducedDim(tse, "NMF"), "loadings")
#' head(loadings)
#' 
NULL

#' @rdname addNMF
#' @export
setGeneric("getNMF", signature = c("x"),
            function(x, ...)
                standardGeneric("getNMF"))

#' @rdname addNMF
#' @export
setGeneric("addNMF", signature = c("x"),
            function(x, ...)
                standardGeneric("addNMF"))

#' @export
#' @rdname addNMF
setMethod("getNMF", "SummarizedExperiment",
    function(x, assay.type = "counts", ...){
        .require_package("NMF")
        # Both NmF and DelayedArray have method seed(). When running
        # NMF::nmf() an error occurs due to wrong method. That is why NMF
        # is first loaded into the session. 
        # Check if NMF package is loaded
        if("NMF" %in% (.packages())){
            detach("package:NMF", unload = TRUE)
        }
        library("NMF")
        .check_assay_present(assay.type, x)
        mat <- t(assay(x, assay.type))
        # Calculate nmf scores for different rank values
        k = c(2, 3, 4, 5, 6, 7, 8, 9, 10)
        nmf_rank <- NMF::nmf(mat, rank = k, ...)
        max_indice <- which.max(nmf_rank$measures$evar)
        # Calculate NMF model for best rank value
        nmf_model <- NMF::nmf(mat, rank = max_indice)
        # store scores
        scores <- nmf_model@fit@W
        # Add loadings as attribute of the scores matrix
        attr(scores, "loadings") <- t(nmf_model@fit@H)
        # Add NMF model as attribute of the scores matrix
        attr(scores, "model") <- nmf_model
        # The NMF package is unloaded
        detach("package:NMF", unload = TRUE)
        # Return scores with loadings, metrics and model as attribute
        return(scores)
    }
)

#' @export
#' @rdname addNMF
setMethod("addNMF", "SummarizedExperiment",
    function(x, k = 2, assay.type = "counts", name = "NMF", ...){
        # Input checks
        if( !.is_a_string(name) ){
            stop("'name' must be a non-empty single character value.",
                call. = FALSE)
        }
        # Fit the model
        nmf <- getNMF(x, k = k, assay.type = assay.type, ...)
        # Add scores matrix with loadings as attribute to reducedDims
        x <- .add_values_to_reducedDims(x, name = name, values = nmf)
        return(x)
    }
)
