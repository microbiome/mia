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
                .Deprecated(msg = paste0("'addTaxonomyTree' is deprecated. ",
                                        "Use 'addHierarchyTree' instead."))
                addHierarchyTree(x,...)
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

#' @rdname deprecate
#' @export
setGeneric("getExperimentCrossAssociation", signature = c("x"),
           function(x, ...)
               standardGeneric("getExperimentCrossAssociation"))

#' @rdname deprecate
#' @export
setMethod("getExperimentCrossAssociation", 
            signature = c(x = "MultiAssayExperiment"),
            function(x, ...){
                .Deprecated(msg = paste0("'getExperimentCrossAssociation' is ",
                                        "deprecated. Use ", 
                                        "'getCrossAssociation' instead."))
                getCrossAssociation(x, ...)
            }
)

#' @rdname deprecate
#' @export
setMethod("getExperimentCrossAssociation", signature = "SummarizedExperiment",
          function(x, ...){
              .Deprecated(msg = paste0("'getExperimentCrossAssociation' is ",
                                       "deprecated. Use ", 
                                       "'getCrossAssociation' instead."))
              getCrossAssociation(x, ...)
          }
)

#' @rdname deprecate
#' @export
setGeneric("testExperimentCrossAssociation", signature = c("x"),
           function(x, ...)
               standardGeneric("testExperimentCrossAssociation"))

#' @rdname deprecate
#' @export
setMethod("testExperimentCrossAssociation", signature = c(x = "ANY"),
          function(x, ...){
              .Deprecated(msg = paste0("'testExperimentCrossAssociation' is ",
                                       "deprecated. Use ", 
                                       "'getCrossAssociation' instead."))
              getCrossAssociation(x, test_significance = TRUE, ...)
          }
)

#' @rdname deprecate
#' @export
setGeneric("testExperimentCrossCorrelation", signature = c("x"),
           function(x, ...)
               standardGeneric("testExperimentCrossCorrelation"))

#' @rdname deprecate
#' @export
setMethod("testExperimentCrossCorrelation", signature = c(x = "ANY"),
          function(x, ...){
              .Deprecated(msg = paste0("'testExperimentCrossCorrelation' is ",
                                       "deprecated. Use ", 
                                       "'getCrossAssociation' instead."))
              getCrossAssociation(x, test_significance = TRUE, ...)
          }
)

#' @rdname deprecate
#' @export
setGeneric("getExperimentCrossCorrelation", signature = c("x"),
           function(x, ...)
               standardGeneric("getExperimentCrossCorrelation"))

#' @rdname deprecate
#' @export
setMethod("getExperimentCrossCorrelation", signature = c(x = "ANY"),
          function(x, ...){
              .Deprecated(msg = paste0("'getExperimentCrossCorrelation' is ",
                                       "deprecated. Use ", 
                                       "'getCrossAssociation' instead."))
              getCrossAssociation(x, ...)
          }
)
