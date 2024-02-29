#' These functions will be deprecated. Please use other functions instead.
#' 
#' @param x A
#'   \code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#'   object.
#'   
#' @param clust.col A single character value indicating the name of the 
#'   \code{rowData} (or \code{colData}) where the data will be stored.
#'   
#' @param ... Additional parameters to use altExps for example
#' @inheritParams bluster::clusterRows
#' @inheritParams runDMN
#' @inheritParams transformAssay
#' 
#' @name deprecate
NULL

#' @rdname deprecate
setGeneric("cluster", signature = c("x"),
            function(x, BLUSPARAM, assay.type = assay_name, 
                    assay_name = "counts", MARGIN = "features", full = FALSE, 
                    name = "clusters", clust.col = "clusters", ...)
                standardGeneric("cluster"))

#' @rdname deprecate
#' @export
#' @importFrom bluster clusterRows
setMethod("cluster", signature = c(x = "SummarizedExperiment"),
            function(x, BLUSPARAM, assay.type = assay_name, 
                    assay_name = "counts", MARGIN = "features", full = FALSE, 
                    name = "clusters", clust.col = "clusters", ...){
                .Deprecated(msg = paste0("'cluster' is deprecated.\n",
                                        "Use 'addCluster' instead."))
                addCluster(x)
            }
)

#' @rdname deprecate
#' @export
setMethod("calculateDMN", signature = c(x = "SummarizedExperiment"),
            function(x, assay.type = assay_name, assay_name = exprs_values, exprs_values = "counts", 
                    transposed = FALSE, ...){
                .Deprecated(old="calculateDMN", new="addCluster", 
                            "Now calculateDMN is deprecated. Use addCluster with DMMParam parameter instead.")
                mat <- assay(x, assay.type)
                if(!transposed){
                    mat <- t(mat)
                }
                calculateDMN(mat, ...)
            }
)

#' @rdname deprecate
#' @importFrom S4Vectors metadata<-
#' @export
runDMN <- function(x, name = "DMN", ...){
    .Deprecated(old="runDMN", new="addCluster", 
                "Now runDMN is deprecated. Use addCluster with DMMParam parameter instead.")
    if(!is(x,"SummarizedExperiment")){
        stop("'x' must be a SummarizedExperiment")
    }
    metadata(x)[[name]] <- calculateDMN(x, ...)
    x
}

#' @rdname deprecate
setGeneric("getDMN", signature = "x",
            function(x, name = "DMN", ...)
                standardGeneric("getDMN"))

#' @rdname deprecate
#' @importFrom DirichletMultinomial laplace AIC BIC
#' @export
setMethod("getDMN", signature = c(x = "SummarizedExperiment"),
            function(x, name = "DMN"){
                .Deprecated(old="getDMN", new="addCluster", 
                            "Now getDMN is deprecated. Use addCluster with DMMParam parameter and full parameter set as true instead.")
                .get_dmn(x, name)
            }
)

#' @rdname deprecate
#' @importFrom DirichletMultinomial laplace AIC BIC
#' @export
setMethod("bestDMNFit", signature = c(x = "SummarizedExperiment"),
            function(x, name = "DMN", type = c("laplace","AIC","BIC")){
                .Deprecated(old="bestDMNFit", new="addCluster", 
                            "Now bestDMNFit is deprecated. Use addCluster with DMMParam parameter and full parameter set as true instead.")
                dmn <- getDMN(x, name)
                fit_FUN <- .get_dmn_fit_FUN(type)
                .get_best_dmn_fit(dmn, fit_FUN)
            }
)

#' @rdname deprecate
setGeneric("getBestDMNFit", signature = "x",
            function(x, name = "DMN", type = c("laplace","AIC","BIC"), ...)
                standardGeneric("getBestDMNFit"))

#' @rdname deprecate
#' @importFrom DirichletMultinomial laplace AIC BIC
#' @export
setMethod("getBestDMNFit", signature = c(x = "SummarizedExperiment"),
            function(x, name = "DMN", type = c("laplace","AIC","BIC")){
                .Deprecated(old="getBestDMNFit", new="addCluster", 
                            "Now getBestDMNFit is deprecated. Use addCluster with DMMParam parameter and full parameter set as true instead.")
                dmn <- getDMN(x, name)
                fit_FUN <- .get_dmn_fit_FUN(type)
                dmn[[.get_best_dmn_fit(dmn, fit_FUN)]]
            }
)
