#' Latent Dirichlet Allocation
#'
#' These functions perform Latent Dirichlet Allocation on data stored in a 
#'  \code{\link[TreeSummarizedExperiment:TreeSummarizedExperiment-class]{TreeSummarizedExperiment}}
#'  object.
#' 
#' @param x a \code{\link[TreeSummarizedExperiment:TreeSummarizedExperiment-class]{TreeSummarizedExperiment}}
#'  object.
#'   
#' @param k \code{Integer scalar}. A number of latent vectors/topics. 
#'  (Default: \code{2})
#' 
#' @param name \code{Character scalar}. The name to be used to store the result 
#'  in the reducedDims of the output. (Default: \code{"LDA"})
#'  
#' @param assay.type \code{Character scalar}. Specifies which assay to use for 
#'  LDA ordination. (Default: \code{"counts"})
#' 
#' @param ... optional arguments passed to \code{\link[topicmodels:LDA]{LDA}}
#' 
#' @return 
#' For \code{getLDA}, the ordination matrix with feature loadings matrix
#'  as attribute \code{"loadings"}.
#'  
#' For \code{addLDA}, a \code{\link[TreeSummarizedExperiment:TreeSummarizedExperiment-class]{TreeSummarizedExperiment}}
#'  object is returned containing the ordination matrix in reducedDims(..., name)
#'  with feature loadings matrix as attribute \code{"loadings"}.
#'  
#' @details 
#' The functions \code{getLDA} and \code{addLDA} internally use 
#'  \code{\link[topicmodels:LDA]{LDA}} to compute the ordination matrix and 
#'  feature loadings.
#'  
#' @name addLDA
#' 
#' @examples
#' data(GlobalPatterns)
#' tse <- GlobalPatterns
#' 
#' # Reduce the number of features 
#' tse <- agglomerateByPrevalence(tse, rank="Phylum")
#' 
#' # Run LDA and add the result to reducedDim(tse, "LDA")
#' tse <- addLDA(tse)
#' 
#' # Extract feature loadings
#' loadings <- attr(reducedDim(tse, "LDA"), "loadings")
#' head(loadings)
NULL

#' @rdname addLDA
#' @export
setGeneric("getLDA", signature = c("x"),
           function(x, ...)
             standardGeneric("getLDA"))

#' @rdname addLDA
#' @export
setGeneric("addLDA", signature = c("x"),
           function(x, ...)
             standardGeneric("addLDA"))

#' @export
#' @rdname addLDA
setMethod("getLDA", "SummarizedExperiment",
    function(x, k = 2, assay.type = "counts", ...){
        .require_package("topicmodels")
        # Input checks
        if( !.is_an_integer(k) ){
            stop("'k' must be an integer.", call. = FALSE)
        }
        .check_assay_present(assay.type, x) 
        df <- as.data.frame(t(assay(x, assay.type)))
        # Estimate LDA model using VEM algorithm
        lda_model <- topicmodels::LDA(df, k, ...)
        # Calculate scores and loadings
        posteriors <- topicmodels::posterior(lda_model, df)
        scores <- t(as.data.frame(posteriors$topics))
        loadings <- t(as.data.frame(posteriors$terms)) 
        # Add loadings as attribute of the scores matrix
        attr(scores, "loadings") <- loadings
        # Add LDA model as attribute of the scores matrix
        attr(scores, "model") <- lda_model
        # Return scores with loadings and model as attribute
        return(scores)
    }
)

#' @export
#' @rdname addLDA
setMethod("addLDA", "SummarizedExperiment",
    function(x, k = 2, assay.type = "counts", name = "LDA", ...){
        # Input checks
        if( !.is_an_integer(k) ){
            stop("'k' must be an integer.", call. = FALSE)
        }
        .check_assay_present(assay.type, x) 
        if( !.is_a_string(name) ){
            stop("'name' must be a non-empty single character value.",
                call. = FALSE)
        }
        scores <- t(getLDA(x, k, assay.type, ...))
        # Add scores matrix with loadings as attribute to reducedDims
        x <- .add_values_to_reducedDims(x, name, values = scores)
        return(x)
    }
)
