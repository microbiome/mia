#' Latent Dirichlet Allocation
#'
#' These functions perform Latent Dirichlet Allocation on data stored in a 
#'  \code{\link[TreeSummarizedExperiment:TreeSummarizedExperiment-class]{TreeSummarizedExperiment}}
#'  object.
#' 
#' @param x a \code{\link[TreeSummarizedExperiment:TreeSummarizedExperiment-class]{TreeSummarizedExperiment}}
#'  object.
#'   
#' @param k Integer, number of latent vectors/topics 
#'   (By default: \code{2}).
#' 
#' @param name The name to be used to store the result in the reducedDims of the
#'  output (By default: \code{name = "LDA"}).
#'  
#' @param assay.type a single \code{character} value for specifying which
#'   assay to use for LDA ordination.
#' 
#' @param ... optional arguments.
#' 
#' @return 
#' For \code{getLDA}, the ordination matrix with feature loadings matrix
#'  as attribute \code{"loadings"}.
#'  
#' For \code{addLDA}, a \code{\link[TreeSummarizedExperiment:TreeSummarizedExperiment-class]{TreeSummarizedExperiment}}
#'  object is returned containing the ordination matrix in reducedDims(..., name)
#'  with feature loadings matrix as attribute \code{"loadings"}.
#'  
#' @name addLDA
#' 
#' @examples
#' data(GlobalPatterns)
#' 
#' # Reduce the number of features 
#' tse <- agglomerateByPrevalence(GlobalPatterns, rank="Phylum", prevalence=0.99, update.tree = TRUE)
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
    function(x, k=2, assay.type = "counts", ...){
        .require_package("topicmodels")
        df <- as.data.frame(t(assay(x, assay.type)))
        # Estimate LDA model using VEM algorithm
        lda_model <- topicmodels::LDA(df, k)
        # Calculate scores and loadings
        posteriors <- topicmodels::posterior(lda_model, df)
        scores <- t(as.data.frame(posteriors$topics))
        loadings <- t(as.data.frame(posteriors$terms)) 
        # Add loadings as attribute of the scores matrix
        attr(scores, "loadings") <- loadings
        # Return scores with loadings as attribute
        return(scores)
    }
)

#' @export
#' @rdname addLDA
setMethod("addLDA", "SummarizedExperiment",
    function(x, k = 2, assay.type = "counts", name = "LDA", ...){
        scores <- t(getLDA(x, k, assay.type, ...))
        # Add scores matrix with loadings as attribute to reducedDims
        x <- .add_values_to_reducedDims(x, name, values = scores)
        return(x)
    }
)
