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
#' @param eval.metric \code{Character scalar}. Specifies evaluation metric that
#' will be used to select the model with the best fit. Must be in
#' \code{c("evar", "sparseness.basis", "sparseness.coef", "rss", "evar",
#'  "silhouette.coef", "silhouette.basis", "cophenetic", "dispersion")}.
#'  (Default: \code{"evar"})
#' 
#' @param ... optional arguments passed to \code{nmf::NMF}.
#' 
#' @return 
#' For \code{getNMF}, the ordination matrix with feature loadings matrix
#'  as attribute \code{"loadings"}.
#'  
#' For \code{addNMF}, a \code{\link[TreeSummarizedExperiment:TreeSummarizedExperiment-class]{TreeSummarizedExperiment}}
#'  object is returned containing the ordination matrix in reducedDims(..., name)
#'  with the following attributes: 
#'  "loadings" which is a matrix containing the feature loadings
#'  "NMF_output" which is the output of function \code{nmf::NMF} 
#'  "best_fit" which is the result of the best fit if k is a vector of integers
#'  
#' @details 
#' The functions \code{getNMF} and \code{addNMF} internally use \code{nmf::NMF} 
#'  compute the ordination matrix and 
#'  feature loadings.
#'  
#' If k is a vector of integers, NMF output is calculated for all the rank 
#'  values contained in k, and the best fit is selected based on 
#'  \code{eval.metric} value.
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
#' tse <- addNMF(tse, k = 2, name = "NMF")
#' 
#' # Extract feature loadings
#' loadings_NMF <- attr(reducedDim(tse, "NMF"), "loadings")
#' head(loadings_NMF)
#' 
#' # Estimate models with number of topics from 2 to 4
#' tse <- addNMF(tse, k = c(2,3,4), name = "NMF_4")
#' 
#' # Extract feature loadings
#' loadings_NMF_4 <- attr(reducedDim(tse, "NMF_4"), "loadings")
#' head(loadings_NMF_4)
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
    function(x, k = 2, assay.type = "counts", eval.metric = "evar", ...){
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
        
        mat <- t(assay(x))
        res <- NMF::nmf(mat, k)
        
        if( is(res, "NMF.rank") ){
          best_fit <- .get_best_nmf_fit(res, eval.metric, ...)
        } else{
          best_fit <- res
        }
        
        scores <- best_fit@fit@W
        loadings <- best_fit@fit@H
        
        attr(scores, "loadings") <- t(loadings)
        attr(scores, "NMF_output") <- res
        if( is(res, "NMF.rank") ){
          attr(scores, "best_fit") <- best_fit
        }
        # The NMF package is unloaded
        detach("package:NMF", unload = TRUE)
        # Return scores with loadings, metrics and model as attribute
        return(scores)
    }
)

#' @export
#' @rdname addNMF
setMethod("addNMF", "SummarizedExperiment",
    function(x, k = 2, assay.type = "counts", eval.metric = "evar",
              name = "NMF", ...){
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

################################ HELP FUNCTIONS ################################
# This function is for evaluating a fit of NMF models
.get_best_nmf_fit <- function(res, eval.metric){
  maximize <- c(
    "sparseness.basis" = TRUE,
    "sparseness.coef" = TRUE,
    "rss" = FALSE,
    "evar" = TRUE,
    "silhouette.coef" = TRUE,
    "silhouette.basis" = TRUE,
    "cophenetic" = TRUE,
    "dispersion" = FALSE)
  maximize <- maximize[ eval.metric ]
  measures <- res[["measures"]]
  values <- measures[[eval.metric]]
  if( maximize ){
    FUN <- which.max
  } else{
    FUN <- which.min
  }
  ind <- FUN(values)
  #
  model <- res[["fit"]][[ind]]
  return(model)
}
