#' These functions will be deprecated. Please use other functions instead.
#' 
#' @param x A
#'   \code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#'   object.
#'   
#' @param ... Additional parameters. See dedicated function.
#' 
#' @name deprecate
NULL

#' @rdname deprecate
setGeneric("cluster", signature = c("x"),
            function(x,...)
                standardGeneric("cluster"))

#' @rdname deprecate
#' @export
#' @importFrom bluster clusterRows
setMethod("cluster", signature = c(x = "SummarizedExperiment"),
            function(x,...){
                .Deprecated(msg = paste0("'cluster' is deprecated. ",
                                        "Use 'addCluster' instead."))
                addCluster(x,...)
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
            function(x,...){
                .Deprecated(msg = paste0("'addTaxonomyTree' is deprecated.",
                                        "Use 'addHierarchyTree' instead."))
                addHierarchyTree(x,...)
            }
)

#' @rdname deprecate
#' @export
setMethod("calculateDMN", signature = c(x = "SummarizedExperiment"),
            function(x,...){
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
runDMN <- function(x,...){
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
            function(x,...)
                standardGeneric("getDMN"))

#' @rdname deprecate
#' @importFrom DirichletMultinomial laplace AIC BIC
#' @export
setMethod("getDMN", signature = c(x = "SummarizedExperiment"),
            function(x,...){
                .Deprecated(msg = paste0("'getDMN' is deprecated. ",
                                        "Use 'addCluster' with DMMParam parameter",
                                        "and full parameter set as true instead."))
                .get_dmn(x,...)
            }
)

#' @rdname deprecate
#' @importFrom DirichletMultinomial laplace AIC BIC
#' @export
setMethod("bestDMNFit", signature = c(x = "SummarizedExperiment"),
            function(x,...){
                .Deprecated(msg = paste0("'bestDMNFit' is deprecated. ",
                                        "Use 'addCluster' with DMMParam parameter",
                                        "and full parameter set as true instead."))
                dmn <- getDMN(x,...)
                fit_FUN <- .get_dmn_fit_FUN(...)
                .get_best_dmn_fit(dmn, fit_FUN)
            }
)

#' @rdname deprecate
setGeneric("getBestDMNFit", signature = "x",
            function(x,...)
                standardGeneric("getBestDMNFit"))

#' @rdname deprecate
#' @importFrom DirichletMultinomial laplace AIC BIC
#' @export
setMethod("getBestDMNFit", signature = c(x = "SummarizedExperiment"),
            function(x,...){
                .Deprecated(msg = paste0("'getBestDMNFit' is deprecated. ",
                                        "Use 'addCluster' with DMMParam parameter",
                                        "and full parameter set as true instead."))
                dmn <- getDMN(x,...)
                fit_FUN <- .get_dmn_fit_FUN(...)
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
            function(x,...){
                .Deprecated(msg = paste0("'taxonomyTree' is deprecated. ",
                                        "Use 'getHierarchyTree' instead."))
                getHierarchyTree(x,...)
            }
)
