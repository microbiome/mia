#' These functions will be deprecated. Please use other functions instead.
#' 
#' @param x A
#'   \code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#'   object.
#'   
#' @param clust.col A single character value indicating the name of the 
#'   \code{rowData} (or \code{colData}) where the data will be stored.
#'   
#' @param exprs_values a single \code{character} value for specifying which
#'   assay to use for calculation.
#'   (Please use \code{assay.type} instead.)
#'   
#' @param transposed Logical scalar, is x transposed with cells in rows?
#' 
#' @param type the type of measure used for the goodness of fit. One of
#'   \sQuote{laplace}, \sQuote{AIC} or \sQuote{BIC}.
#'   
#' @inheritParams bluster::clusterRows
#' @inheritParams runDMN
#' @inheritParams transformAssay
#' 
#' @param x a \code{\link{SummarizedExperiment}} object -
#' 
#' @param ... Additional parameters to use altExps for example
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
                .Deprecated(msg = paste0("'cluster' is deprecated. ",
                                        "Use 'addCluster' instead."))
                addCluster(x)
            }
)

#' @rdname deprecate
setGeneric("addTaxonomyTree",
            signature = "x",
            function(x, ...)
                standardGeneric("addTaxonomyTree"))

#' @rdname deprecate
#' @export
setMethod("addTaxonomyTree", signature = c(x = "SummarizedExperiment"),
            function(x){
                .Deprecated(msg = paste0("'addTaxonomyTree' is deprecated.",
                                        "Use 'addHierarchyTree' instead."))
                addHierarchyTree(x)
            }
)

#' @rdname deprecate
#' @export
setMethod("calculateDMN", signature = c(x = "SummarizedExperiment"),
            function(x, assay.type = assay_name, assay_name = exprs_values, exprs_values = "counts", 
                    transposed = FALSE, ...){
                .Deprecated(msg = paste0("'calculateDMN' is deprecated. ",
                                        "Use 'addCluster' with DMMParam parameter instead."))
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
    .Deprecated(msg = paste0("'runDMN' is deprecated. ",
                            "Use 'addCluster' with DMMParam parameter instead."))
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
                .Deprecated(msg = paste0("'getDMN' is deprecated. ",
                                        "Use 'addCluster' with DMMParam parameter",
                                        "and full parameter set as true instead."))
                .get_dmn(x, name)
            }
)

#' @rdname deprecate
#' @importFrom DirichletMultinomial laplace AIC BIC
#' @export
setMethod("bestDMNFit", signature = c(x = "SummarizedExperiment"),
            function(x, name = "DMN", type = c("laplace","AIC","BIC")){
                .Deprecated(msg = paste0("'bestDMNFit' is deprecated. ",
                                        "Use 'addCluster' with DMMParam parameter",
                                        "and full parameter set as true instead."))
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
                .Deprecated(msg = paste0("'getBestDMNFit' is deprecated. ",
                                        "Use 'addCluster' with DMMParam parameter",
                                        "and full parameter set as true instead."))
                dmn <- getDMN(x, name)
                fit_FUN <- .get_dmn_fit_FUN(type)
                dmn[[.get_best_dmn_fit(dmn, fit_FUN)]]
            }
)

#' @rdname deprecate
setGeneric("taxonomyTree",
            signature = "x",
            function(x, ...)
                standardGeneric("taxonomyTree"))

#' @rdname deprecate
#' @export
setMethod("taxonomyTree", signature = c(x = "SummarizedExperiment"),
            function(x){
                .Deprecated(msg = paste0("'taxonomyTree' is deprecated. ",
                                        "Use 'getHierarchyTree' instead."))
                getHierarchyTree(x)
            }
)
