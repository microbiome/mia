#' These functions will be deprecated. Please use other functions instead.
#' 
#' @param x A
#'   \code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#'   object.
#' @param value a matrix to store as the \sQuote{relabundance} assay
#' 
#' @param mat An abundance matrix
#' 
#' @param tree A phylogenetic tree
#'    
#' @param ... Additional parameters. See dedicated function.
#' 
#' @name deprecate
NULL

#' @rdname deprecate
#' @export
setGeneric("cluster", signature = c("x"),
            function(x,...)
                standardGeneric("cluster"))

#' @rdname deprecate
#' @export
setMethod("cluster", signature = c(x = "SummarizedExperiment"),
            function(x,...){
                .Deprecated(msg = paste0("'cluster' is deprecated. ",
                                        "Use 'addCluster' instead."))
                addCluster(x,...)
            }
)

#' @rdname deprecate
#' @export
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
#' @export
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
setGeneric("mergeRows",
            signature = "x",
            function(x, ...)
                standardGeneric("mergeRows"))

#' @rdname deprecate
#' @export
setMethod("mergeRows", signature = c(x = "SummarizedExperiment"),
            function(x, ...){
                .Deprecated(msg = paste0("'mergeRows' is deprecated. ",
                                        "Use 'agglomerateByVariable' with ",
                                        "parameter by = 'rows' instead."))
                agglomerateByVariable(x, by = "rows", ...)
            }
)

#' @rdname deprecate
#' @export
setMethod("mergeRows", signature = c(x = "TreeSummarizedExperiment"),
          function(x, ...){
              .Deprecated(msg = paste0("'mergeRows' is deprecated. ",
                                       "Use 'agglomerateByVariable' with ",
                                        "parameter by = 'rows' instead."))
              agglomerateByVariable(x, by = "rows", ...)
          }
)

#' @rdname deprecate
#' @export
setGeneric("mergeCols",
            signature = "x",
            function(x, ...)
                standardGeneric("mergeCols"))

#' @rdname deprecate
#' @export
setMethod("mergeCols", signature = c(x = "SummarizedExperiment"),
            function(x, ...){
                .Deprecated(msg = paste0("'mergeCols' is deprecated. ",
                                        "Use 'agglomerateByVariable' with ", 
                                        "parameter by = 'cols' instead."))
                agglomerateByVariable(x, by = "cols", ...)
            }
)

#' @rdname deprecate
#' @export
setMethod("mergeCols", signature = c(x = "TreeSummarizedExperiment"),
          function(x, ...){
              .Deprecated(msg = paste0("'mergeCols' is deprecated. ",
                                        "Use 'agglomerateByVariable' with ",
                                        "parameter by = 'cols' instead."))
              agglomerateByVariable(x, by = "cols", ...)
          }
)

#' @rdname deprecate
#' @export
setGeneric("mergeFeatures",
            signature = "x",
            function(x, ...)
                standardGeneric("mergeFeatures"))

#' @rdname deprecate
#' @export
setMethod("mergeFeatures", signature = c(x = "SummarizedExperiment"),
            function(x, ...){
                .Deprecated(msg = paste0("'mergeFeatures' is deprecated. ",
                                        "Use 'agglomerateByVariable' with ", 
                                        "parameter by = 'rows' instead."))
                agglomerateByVariable(x, by = "rows", ...)
            }
)

#' @rdname deprecate
#' @export
setMethod("mergeFeatures", signature = c(x = "TreeSummarizedExperiment"),
          function(x, ...){
              .Deprecated(msg = paste0("'mergeFeatures' is deprecated. ",
                                       "Use 'agglomerateByVariable' with ",
                                        "parameter by = 'rows' instead."))
              agglomerateByVariable(x, by = "rows", ...)
          }
)

#' @rdname deprecate
#' @export
setGeneric("mergeSamples",
            signature = "x",
            function(x, ...)
                standardGeneric("mergeSamples"))

#' @rdname deprecate
#' @export
setMethod("mergeSamples", signature = c(x = "SummarizedExperiment"),
            function(x, ...){
                .Deprecated(msg = paste0("'mergeSamples' is deprecated. ",
                                        "Use 'agglomerateByVariable' with ",
                                        "parameter by = 'cols' instead."))
                agglomerateByVariable(x, by = "cols", ...)
            }
)

#' @rdname deprecate
#' @export
setMethod("mergeSamples", signature = c(x = "TreeSummarizedExperiment"),
          function(x, ...){
              .Deprecated(msg = paste0("'mergeSamples' is deprecated. ",
                                       "Use 'agglomerateByVariable' with ", 
                                        "parameter by = 'cols' instead."))
              agglomerateByVariable(x, by = "cols", ...)
          }
)

#' @rdname deprecate
#' @export
setGeneric("mergeFeaturesByRank",
           signature = "x",
           function(x, ...)
               standardGeneric("mergeFeaturesByRank"))

#' @rdname deprecate
#' @export
setMethod("mergeFeaturesByRank", signature = c(x = "SummarizedExperiment"),
          function(x, ...){
              .Deprecated(msg = paste0("'mergeFeaturesByRank' is deprecated. ",
                                        "Use 'agglomerateByRank' instead."))
              x <- agglomerateByRank(x, ...)
              x
          }
)

#' @rdname deprecate
#' @export
setMethod("mergeFeaturesByRank", signature = c(x = "SingleCellExperiment"),
          function(x, ...){
              .Deprecated(msg = paste0("'mergeFeaturesByRank' is deprecated. ",
                                        "Use 'agglomerateByRank' instead."))
              x <- agglomerateByRank(x, ...)
              x
          }
)

#' @rdname deprecate
#' @export
setGeneric("mergeFeaturesByPrevalence", signature = "x",
           function(x, ...)
               standardGeneric("mergeFeaturesByPrevalence"))

#' @rdname deprecate
#' @export
setMethod("mergeFeaturesByPrevalence", signature = c(x = "SummarizedExperiment"),
          function(x, ...){
              .Deprecated(msg = paste0(
                  "'mergeFeaturesByPrevalence' is deprecated. ",
                  "Use agglomerateByPrevalence instead."))
              x <- agglomerateByPrevalence(x, ...)
              x 
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
setMethod("mergeFeaturesByRank", signature = c(x = "TreeSummarizedExperiment"),
          function(x, ...){
              .Deprecated(msg = paste0("'mergeFeaturesByRank' is deprecated. ",
                                        "Use 'agglomerateByRank' instead."))
              x <- agglomerateByRank(x, ...)
              x
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
              getCrossAssociation(x, test.signif = TRUE, ...)
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
              getCrossAssociation(x, test.signif = TRUE, ...)
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

#' @rdname deprecate
#' @export
loadFromBiom <- function(...) {
    .Deprecated(msg = paste0("'loadFromBiom' is deprecated.",
                            " Use 'importBIOM' instead."))
    importBIOM(...)
}

#' @rdname deprecate
#' @export
loadFromQIIME2 <- function(...) {
    .Deprecated(msg = paste0("'loadFromQIIME2' is deprecated.",
                            " Use 'importQIIME2' instead."))
    importQIIME2(...)
}

#' @rdname deprecate
#' @export
readQZA <- function(...) {
    .Deprecated(msg = paste0("'readQZA' is deprecated.",
                            " Use 'importQZA' instead."))
    importQZA(...)
}

#' @rdname deprecate
#' @export
loadFromMothur <- function(...) {
    .Deprecated(msg = paste0("'loadFromMothur' is deprecated.",
                            " Use 'importMothur' instead."))
    importMothur(...)
}

#' @rdname deprecate
#' @export   
loadFromMetaphlan <- function(...) {
    .Deprecated(msg = paste0("'loadFromMetaphlan' is deprecated.",
                            " Use 'importMetaPhlAn' instead."))
    importMetaPhlAn(...)
}

#' @rdname deprecate
#' @export    
loadFromHumann <- function(...) {
    .Deprecated(msg = paste0("'loadFromHumann' is deprecated.",
                            " Use 'importHUMAnN' instead."))
    importHUMAnN(...)
}

#' @rdname deprecate
#' @export
setGeneric("countDominantFeatures", signature = c("x"),
            function(x, ...)
                standardGeneric("countDominantFeatures"))

#' @rdname deprecate
#' @export
setMethod("countDominantFeatures", signature = c(x = "SummarizedExperiment"),
            function(x, ...){
                .Deprecated(msg = "'countDominantFeatures' function is deprecated. ",
                            "Use 'summarizeDominance' instead.")
                summarizeDominance(x, ...)
            }
)

#' @rdname deprecate
#' @export
setGeneric(
    "subsetByRareTaxa", signature = c("x"), function(x, ...) 
        standardGeneric("subsetByRareTaxa"))

#' @rdname deprecate
#' @export
setMethod(
    "subsetByRareTaxa", signature = c(x = "ANY"), function(x, ...){
        .Deprecated(
            msg = "'subsetByRareTaxa' is deprecated. ",
            "Use 'subsetByRare' instead.")
        subsetByRare(x, ...)
    }
)

#' @rdname deprecate
#' @export
setGeneric(
    "subsetByRareFeatures", signature = c("x"), function(x, ...) 
        standardGeneric("subsetByRareFeatures"))

#' @rdname deprecate
#' @export
setMethod(
    "subsetByRareFeatures", signature = c(x = "ANY"), function(x, ...){
        .Deprecated(
            msg = "'subsetByRareFeatures' is deprecated. ",
            "Use 'subsetByRare' instead.")
        subsetByRare(x, ...)
    }
)

#' @rdname deprecate
#' @export
setGeneric(
    "subsetByPrevalentTaxa", signature = c("x"), function(x, ...) 
        standardGeneric("subsetByPrevalentTaxa"))

#' @rdname deprecate
#' @export
setMethod(
    "subsetByPrevalentTaxa", signature = c(x = "ANY"), function(x, ...){
        .Deprecated(
            msg = "'subsetByPrevalentTaxa' is deprecated. Use ",
            "'subsetByPrevalent' instead.")
        subsetByPrevalent(x, ...)
    }
)

#' @rdname deprecate
#' @export
setGeneric(
    "subsetByPrevalentFeatures", signature = c("x"), function(x, ...) 
        standardGeneric("subsetByPrevalentFeatures"))

#' @rdname deprecate
#' @export
setMethod(
    "subsetByPrevalentFeatures", signature = c(x = "ANY"), function(x, ...){
        .Deprecated(
            msg = "'subsetByPrevalentFeatures' is deprecated. Use ",
            "'subsetByPrevalent' instead.")
        subsetByPrevalent(x, ...)
    }
)

#' @rdname deprecate
#' @export
setGeneric("countDominantTaxa", signature = c("x"),
           function(x, ...)
               standardGeneric("countDominantTaxa"))

#' @rdname deprecate
#' @export
setMethod("countDominantTaxa", signature = c(x = "SummarizedExperiment"),
          function(x, ...){
              .Deprecated(msg = "'countDominantTaxa' function is deprecated. ",
                          "Use 'summarizeDominance' instead.")
              summarizeDominance(x, ...)
          }
)

#' @rdname deprecate
#' @export
setGeneric("full_join", signature = c("x"),
           function(x, ...)
               standardGeneric("full_join"))

#' @rdname deprecate
#' @export
setMethod("full_join", signature = c(x = "ANY"),
          function(x, ...){
              .Deprecated(msg = paste0("'full_join' is deprecated. ",
                                       "Use 'mergeSEs' with 'join = full' ",
                                       "instead."))
              mergeSEs(x, join = "full", ...)
          }
)

#' @rdname deprecate
#' @export
setGeneric("inner_join", signature = c("x"),
           function(x, ...)
               standardGeneric("inner_join"))

#' @rdname deprecate
#' @export
setMethod("inner_join", signature = c(x = "ANY"),
            function(x, ...){
                .Deprecated(msg = paste0("'inner_join' is deprecated. ",
                                        "Use 'mergeSEs' with 'join = inner' ",
                                        "instead."))
                mergeSEs(x, join = "inner", ...)
            }
)

#' @rdname deprecate
#' @export
setGeneric("left_join", signature = c("x"),
           function(x, ...)
               standardGeneric("left_join"))

#' @rdname deprecate
#' @export
setMethod("left_join", signature = c(x = "ANY"),
          function(x, ...){
              .Deprecated(msg = paste0("'left_join' is deprecated. ",
                                        "Use 'mergeSEs' with 'join = left' ",
                                        "instead."))
              mergeSEs(x, join = "left", ...)
          }
)

#' @rdname deprecate
#' @export
setGeneric("right_join", signature = c("x"),
           function(x, ...)
               standardGeneric("right_join"))

#' @rdname deprecate
#' @export
setMethod("right_join", signature = c(x = "ANY"),
            function(x, ...){
                .Deprecated(msg = paste0("'right_join' is deprecated. ",
                                        "Use 'mergeSEs' with 'join = right' ",
                                        "instead."))
                mergeSEs(x, join = "right", ...)
            }
)

#' @rdname deprecate
#' @export    
plotNMDS <- function(x, ...){
    .Deprecated(msg = paste0("'plotNMDS' is deprecated. ",
                             "Use 'scater::plotReducedDim' with ",
                             "dimred = 'NMDS' instead."))
    plotReducedDim(x, ncomponents = 2, dimred = "NMDS",...)
}

#' @rdname deprecate
#' @export
setGeneric("estimateDivergence",signature = c("x"),
            function(x, ...)
                standardGeneric("estimateDivergence"))

#' @rdname deprecate
#' @export
setMethod("estimateDivergence", signature = c(x="SummarizedExperiment"),
            function(x, ...){
                .Deprecated(msg = paste0("'estimateDivergence' is deprecated. ",
                                        "Use 'addDivergence' instead."))
                addDivergence(x, ...)
            }
)

#' @rdname deprecate
#' @export
setGeneric("meltAssay",signature = "x",
            function(x, ...)
                standardGeneric("meltAssay"))

#' @rdname deprecate
#' @export
setMethod("meltAssay", signature = c(x="SummarizedExperiment"),
          function(x, ...){
              .Deprecated(msg = paste0("'meltAssay' is deprecated. ",
                                       "Use 'meltSE' instead."))
              meltSE(x, ...)
          }
)

#' @rdname deprecate
#' @export
setGeneric(
    "transformSamples", signature = c("x"), function(x,...)
        standardGeneric("transformSamples"))

#' @rdname deprecate
#' @export
setMethod(
    "transformSamples", signature = c(x = "SummarizedExperiment"),
    function(x, ...){
        .Deprecated(
            "'transformSamples' is deprecated. Use 'transformAssay' instead.")
        transformAssay(x, MARGIN = "samples", ...)
    }
)

#' @rdname deprecate
#' @export
setGeneric(
    "ZTransform", signature = c("x"), function(x, ...)
        standardGeneric("ZTransform"))

#' @rdname deprecate
#' @export
setMethod(
    "ZTransform", signature = c(x = "SummarizedExperiment"), function(x, ...){
        .Deprecated(
            "'Ztransform' is deprecated. Use 'transformAssay' instead.")
        transformAssay(x, method = "standardize", MARGIN = "features", ...)
    }
)

#' @rdname deprecate
#' @export
setGeneric(
    "relAbundanceCounts", signature = c("x"), function(x, ...)
        standardGeneric("relAbundanceCounts"))

#' @rdname deprecate
#' @export
setMethod(
    "relAbundanceCounts", signature = c(x = "SummarizedExperiment"),
    function(x, ...){
        .Deprecated(
            "'relAbundanceCounts' is deprecated. ",
            "Use 'transformAssay' instead.")
        transformAssay(x, method = "relabundance", MARGIN = "samples", ...)
    }
)

#' @rdname deprecate
#' @export
transformCounts <- function(x, ...){
    .Deprecated(
        "'transformCounts' is deprecated. Use 'transformAssay' instead.")
    return(transformAssay(x,...))
}

#' @rdname deprecate
#' @export
setGeneric(
    "transformFeatures", signature = c("x"), function(x, ...)
        standardGeneric("transformFeatures"))

#' @rdname deprecate
#' @export
setMethod(
    "transformFeatures", signature = c(x = "SummarizedExperiment"),
    function(x, ...){
        .Deprecated(
            "'transformFeatures' is deprecated. Use 'transformAssay' instead.")
        transformAssay(x, MARGIN = "features", ...)
    }
)

#' @rdname deprecate
#' @export
setGeneric(
    "getUniqueFeatures", signature = c("x"), function(x, ...)
        standardGeneric("getUniqueFeatures"))

#' @rdname deprecate
#' @export
setMethod(
    "getUniqueFeatures", signature = c(x = "SummarizedExperiment"),
    function(x,...){
        .Deprecated(
            msg = "'getUniqueFeatures' is deprecated. Use 'getUnique' instead.")
        getUnique(x,...)
    }
)

#' @rdname deprecate
#' @export
setGeneric(
    "getUniqueTaxa", signature = c("x"), function(x, ...)
        standardGeneric("getUniqueTaxa"))

#' @rdname deprecate
#' @export
setMethod(
    "getUniqueTaxa", signature = c(x = "SummarizedExperiment"), function(x,...){
        .Deprecated(
            msg = "'getUniqueTaxa' is deprecated. Use 'getUnique' instead.")
        getUnique(x,...)
    }
)

#' @rdname deprecate
#' @export
setGeneric(
    "getTopFeatures", signature = "x", function(x,...)
        standardGeneric("getTopFeatures"))

#' @rdname deprecate
#' @export
setMethod(
    "getTopFeatures", signature = c(x = "SummarizedExperiment"),
    function(x,...){
        .Deprecated(
            msg = "'getTopFeatures' is deprecated. Use 'getTop' instead.")
        getTop(x,...)
    }
)

#' @rdname deprecate
#' @export
setGeneric(
    "getTopTaxa", signature = "x", function(x,...)
        standardGeneric("getTopTaxa"))

#' @rdname deprecate
#' @export
setMethod(
    "getTopTaxa", signature = c(x = "SummarizedExperiment"), function(x,...){
        .Deprecated(msg = "'getTopTaxa' is deprecated. Use 'getTop' instead.")
        getTop(x,...)
    }
)

#' @rdname deprecate
#' @export
setGeneric(
    "getRareFeatures", signature = "x", function(x, ...)
        standardGeneric("getRareFeatures"))

#' @rdname deprecate
#' @export
setMethod(
    "getRareFeatures", signature = c(x = "ANY"), function(x,...){
        .Deprecated(
            msg = "'getRareFeatures' is deprecated. Use 'getRare' instead.")
        getRare(x,...)
    }
)

#' @rdname deprecate
#' @export
setMethod(
    "getRareFeatures", signature = c(x = "SummarizedExperiment"),
    function(x,...){
        .Deprecated(
            msg = "'getRareFeatures' is deprecated. Use 'getRare' instead.")
        getRare(x,...)
    }
)

#' @rdname deprecate
#' @export
setGeneric(
    "getRareTaxa", signature = "x", function(x, ...)
        standardGeneric("getRareTaxa"))

#' @rdname deprecate
#' @export
setMethod(
    "getRareTaxa", signature = c(x = "ANY"), function(x,...){
        .Deprecated(
            msg = "'getRareTaxa' is deprecated. Use 'getRare' instead.")
        getRare(x,...)
    }
)

#' @rdname deprecate
#' @export
setMethod(
    "getRareTaxa", signature = c(x = "SummarizedExperiment"), function(x,...){
        .Deprecated(msg = "'getRareTaxa' is deprecated. Use 'getRare' instead.")
        getRare(x,...)
    }
)

#' @rdname deprecate
#' @export
setGeneric(
    "getPrevalentFeatures", signature = "x", function(x, ...)
        standardGeneric("getPrevalentFeatures"))

#' @rdname deprecate
#' @export
setMethod(
    "getPrevalentFeatures", signature = c(x = "ANY"), function(x,...){
        .Deprecated(
            msg = "'getPrevalentFeatures' is deprecated. ",
            "Use 'getPrevalent' instead.")
        getPrevalent(x,...)
    }
)

#' @rdname deprecate
#' @export
setMethod(
    "getPrevalentFeatures", signature = c(x = "SummarizedExperiment"),
    function(x,...){
        .Deprecated(
            msg = "'getPrevalentFeatures' is deprecated. ",
            "Use 'getPrevalent' instead.")
        getPrevalent(x,...)
    }
)

#' @rdname deprecate
#' @export
setGeneric(
    "getPrevalentTaxa", signature = "x", function(x, ...)
        standardGeneric("getPrevalentTaxa"))

#' @rdname deprecate
#' @export
setMethod(
    "getPrevalentTaxa", signature = c(x = "ANY"), function(x,...){
        .Deprecated(
            msg = "'getPrevalentTaxa' is deprecated. ",
            "Use 'getPrevalent' instead.")
        getPrevalent(x,...)
    }
)

#' @rdname deprecate
#' @export
setMethod(
    "getPrevalentTaxa", signature = c(x = "SummarizedExperiment"),
    function(x,...){
        .Deprecated(
            msg = "'getPrevalentTaxa' is deprecated. ",
            "Use 'getPrevalent' instead.")
        getPrevalent(x,...)
    }
)

#' @rdname deprecate
#' @export
setGeneric(
    "subsampleCounts", signature = c("x"), function(x, ...)
        standardGeneric("subsampleCounts"))

#' @rdname deprecate
#' @export
setMethod(
    "subsampleCounts", signature = c(x = "SummarizedExperiment"),
    function(x, ...){
        .Deprecated(
            msg = "'subsampleCounts' is deprecated. ",
            "Use 'rarefyAssay' instead.")
        rarefyAssay(x, ...)
    }
)

#' @rdname deprecate
#' @export
setGeneric(
    "addPerSampleDominantFeatures", signature = c("x"), function(x, ...) 
        standardGeneric("addPerSampleDominantFeatures"))

#' @rdname deprecate
#' @export
setMethod(
    "addPerSampleDominantFeatures", signature = c(x = "SummarizedExperiment"),
    function(x, ...){
        .Deprecated(
            msg = "'addPerSampleDominantFeatures' is deprecated.",
            " Use 'addDominant' instead.")
        addDominant(x, ...)
    }
)

#' @rdname deprecate
#' @export
setGeneric(
    "addPerSampleDominantTaxa", signature = c("x"), function(x, ...) 
        standardGeneric("addPerSampleDominantTaxa"))

#' @rdname deprecate
#' @export
setMethod(
    "addPerSampleDominantTaxa", signature = c(x = "SummarizedExperiment"),
    function(x, ...){
        .Deprecated(
            msg = "'addPerSampleDominantTaxa' is deprecated. ",
            "Use 'addDominant' instead.")
        addDominant(x, ...)
    }
)

#' @rdname deprecate
#' @export
setGeneric(
    "perSampleDominantFeatures", signature = c("x"), function(x, ...) 
        standardGeneric("perSampleDominantFeatures"))

#' @rdname deprecate
#' @export
setMethod(
    "perSampleDominantFeatures", signature = c(x = "SummarizedExperiment"),
    function(x, ...){
        .Deprecated(
            msg = "'perSampleDominantFeatures' is deprecated. ",
            "Use 'getDominant' instead.")
        getDominant(x, ...)
    }
)

#' @rdname deprecate
#' @export
setGeneric(
    "perSampleDominantTaxa", signature = c("x"), function(x, ...) 
        standardGeneric("perSampleDominantTaxa"))

#' @rdname deprecate
#' @export
setMethod(
    "perSampleDominantTaxa", signature = c(x = "SummarizedExperiment"),
    function(x, ...){
        .Deprecated(
            msg = "'perSampleDominantTaxa' is deprecated. ",
            "Use 'getDominant' instead.")
        getDominant(x, ...)
    }
)

#' @rdname deprecate
#' @export
setGeneric("makePhyloseqFromTreeSE", signature = c("x"),
           function(x, ...)
             standardGeneric("makePhyloseqFromTreeSE"))

#' @rdname deprecate
#' @export
setMethod("makePhyloseqFromTreeSE", signature = c(x = "SummarizedExperiment"),
    function(x){
        .Deprecated(msg = paste0(
            "'makeTreeSummarizedExperimentFromPhyloseq' is deprecated.",
            " Use 'convertFromPhyloseq' instead."))
        convertFromPhyloseq(x)
    }
)

#' @rdname deprecate
#' @export
setMethod("makePhyloseqFromTreeSE", 
          signature = c(x = "TreeSummarizedExperiment"),
    function(x){
        .Deprecated(msg = paste0(
            "'makeTreeSummarizedExperimentFromPhyloseq' is deprecated.",
            " Use 'convertFromPhyloseq' instead."))
        convertFromPhyloseq(x)
    }
)

#' @rdname deprecate
#' @export
setGeneric("makePhyloseqFromTreeSummarizedExperiment", signature = c("x"),
           function(x)
             standardGeneric("makePhyloseqFromTreeSummarizedExperiment"))

#' @rdname deprecate
#' @export
setMethod("makePhyloseqFromTreeSummarizedExperiment", signature = c(x = "ANY"), 
          function(x){
            .Deprecated(msg = paste0(
              "'makePhyloseqFromTreeSummarizedExperiment' is deprecated.",
              " Use 'convertToPhyloseq' instead."))
            convertToPhyloseq(x)
          }
)

#' @rdname deprecate
#' @export
setGeneric("makeTreeSummarizedExperimentFromPhyloseq", signature = c("x"),
           function(x)
             standardGeneric("makeTreeSummarizedExperimentFromPhyloseq"))

#' @rdname deprecate
#' @export
setMethod("makeTreeSummarizedExperimentFromPhyloseq", signature = c(x = "ANY"), 
    function(x){
        .Deprecated(msg = paste0(
            "'makeTreeSummarizedExperimentFromPhyloseq' is deprecated.",
            " Use 'convertFromPhyloseq' instead."))
        convertFromPhyloseq(x)
    }
)

#' @rdname deprecate
#' @export
makeTreeSEFromBiom <- function(...){
  .Deprecated(msg = paste0(
    "'makeTreeSEFromBiom' is deprecated.",
    " Use 'convertFromBIOM' instead."))
  convertFromBIOM(...)
}

#' @rdname deprecate
#' @export
makeTreeSummarizedExperimentFromBiom <- function(...){
  .Deprecated(msg = paste0(
    "'makeTreeSummarizedExperimentFromBiom' is deprecated.",
    " Use 'convertFromBIOM' instead."))
  convertFromBIOM(...)
}

#' @rdname deprecate
#' @export
makeTreeSEFromDADA2 <- function(...) {
  .Deprecated(msg = paste0(
    "'makeTreeSEFromDADA2' is deprecated.",
    " Use 'convertFromDADA2' instead."))
  convertFromDADA2(...)
}

#' @rdname deprecate
#' @export
makeTreeSummarizedExperimentFromDADA2 <- function(...) {
  .Deprecated(msg = paste0(
    "'makeTreeSummarizedExperimentFromDADA2' is deprecated.",
    " Use 'convertFromDADA2' instead."))
  convertFromDADA2(...)
}

#' @rdname deprecate
#' @export
makeTreeSEFromPhyloseq <- function(x) {
  .Deprecated(msg = paste0(
    "'makeTreeSEFromPhyloseq' is deprecated.",
    " Use 'convertFromPhyloseq' instead."))
  convertFromPhyloseq(x)
}

#' @rdname deprecate
#' @export
setGeneric(
    "estimateEvenness", signature = c("x"), function(x, ...)
        standardGeneric("estimateEvenness"))

#' @rdname deprecate
#' @export
setMethod(
    "estimateEvenness", signature = c(x = "ANY"), function(x, ...){
        .Deprecated(msg = paste0(
            "'estimateEvenness' is deprecated. Use 'addAlpha' instead."))
        addAlpha(x, ...)
    }
)

#' @rdname deprecate
#' @export
setGeneric(
    "estimateRichness", signature = c("x"),
    function(x, ...) standardGeneric("estimateRichness"))

#' @rdname deprecate
#' @export
setMethod(
    "estimateRichness", signature = c(x="ANY"),
    function(x, ...){
        .Deprecated(msg = paste0(
          "'estimateRichness' is deprecated. Use 'addAlpha' instead."))
        addAlpha(x, ...)
    }
)

#' @rdname deprecate
#' @export
setGeneric(
    "estimateDiversity", signature = c("x"),
    function(x, ...) standardGeneric("estimateDiversity"))

#' @rdname deprecate
#' @export
setMethod(
    "estimateDiversity", signature = c(x = "ANY"), function(x, ...){
        .Deprecated(msg = paste0(
            "'estimateDiversity' is deprecated. Use 'addAlpha' instead."))
        addAlpha(x, ...)
    }
)

#' @rdname deprecate
#' @export
setGeneric(
    "estimateFaith", signature = c("x"), function(x, ...)
        standardGeneric("estimateFaith"))

#' @rdname deprecate
#' @export
setMethod(
    "estimateFaith", signature = c(x="ANY"),
    function(x, ...){
        .Deprecated(msg = paste0(
          "'estimateFaith' is deprecated. Use 'addAlpha' instead."))
        addAlpha(x, ...)
    }
)

#' @rdname deprecate
#' @export
setGeneric(
    "estimateDominance", signature = c("x"),
    function(x, ...) standardGeneric("estimateDominance"))

#' @rdname deprecate
#' @export
setMethod(
    "estimateDominance", signature = c(x="ANY"),
    function(x, ...){
        .Deprecated(msg = paste0(
          "'estimateDominance' is deprecated. Use 'addAlpha' instead."))
        addAlpha(x, ...)
    })
#' @rdname deprecate
#' @export
setGeneric("subsetSamples", signature = "x",
           function(x, ...)
               standardGeneric("subsetSamples"))
#' @rdname deprecate
#' @export
setGeneric("subsetFeatures", signature = "x",
           function(x, ...)
               standardGeneric("subsetFeatures"))
#' @rdname deprecate
#' @export
setGeneric("subsetTaxa", signature = "x",
           function(x, ...)
               standardGeneric("subsetTaxa"))
#' @rdname deprecate
#' @export
setMethod("subsetSamples", signature = "SummarizedExperiment",
          function(x, ...){
              .Deprecated(msg = paste0("subsetSamples is deprecated. Please ",
                                       "use '[]' for subsetting instead."))
              subset_args <- .get_subset_args(x, ...)
              x[subset_args$columns,subset_args$rows]
          }
)

#' @rdname deprecate
#' @export
setMethod("subsetFeatures", signature = "SummarizedExperiment",
          function(x, ...){
              .Deprecated(msg = paste0("subsetFeatures is deprecated. Please",
                                       " use '[]' for subsetting instead."))
              subset_args <- .get_subset_args(x, ...)
              x[subset_args$rows, subset_args$columns]
          }
)

#' @rdname deprecate
#' @export
setMethod("subsetTaxa", signature = "SummarizedExperiment",
          function(x, ...){
              .Deprecated(msg = paste0("subsetFeatures is deprecated. Please",
                                       " use '[]' for subsetting instead."))
              subsetFeatures(x, ...)
          }
)
#' @rdname deprecate
#' @export
setGeneric("relabundance", signature = c("x"),
           function(x, ...) standardGeneric("relabundance"))

#' @rdname deprecate
#' @export
setGeneric("relabundance<-", signature = c("x"),
           function(x, value) standardGeneric("relabundance<-"))

#' @rdname deprecate
#' @importFrom SummarizedExperiment assays
#' @export
setMethod("relabundance", signature = c(x = "SummarizedExperiment"),
    function(x){
        .Deprecated(msg = paste0("'relabundance' is deprecated\n",
                                 "Use 'assay(x, 'relabundance')' instead."))
        assays(x)[["relabundance"]]
    }
)

#' @rdname deprecate
#' @importFrom SummarizedExperiment assays<-
#' @export
setReplaceMethod("relabundance", signature = c(x = "SummarizedExperiment"),
    function(x, value){
        .Deprecated(msg = paste0("'relabundance' is deprecated\n",
                                 "Use 'assay(x, 'relabundance')' instead."))
        assays(x)[["relabundance"]] <- value
        x
    }
)

#' @rdname deprecate
#' @export
setGeneric("runOverlap", signature = c("x"),
          function(x, ...)
            standardGeneric("runOverlap"))

#' @rdname deprecate
#' @export
setMethod("runOverlap", signature = c(x = "SummarizedExperiment"),
    function(x, ...){
        .Deprecated(msg = paste0("'runOverlap' is deprecated\n",
                               "Use 'addDissimilarity' with parameter ",
                               "method = 'overlap' instead."))
        addDissimilarity(x, method = "overlap", ...)
    }
)

#' @rdname deprecate
#' @export
setGeneric("calculateOverlap", signature = c("x"),
           function(x, ...)
             standardGeneric("calculateOverlap"))

#' @rdname deprecate
#' @export
setMethod("calculateOverlap", signature = c(x = "SummarizedExperiment"),
          function(x, ...){
            .Deprecated(msg = paste0("'calculateOverlap' is deprecated\n",
                                     "Use 'getDissimilarity' with parameter ",
                                     "method = 'overlap' instead."))
            getDissimilarity(x, method = "overlap", ...)
          }
)

#' @rdname deprecate
#' @export
setGeneric("calculateJSD", signature = c("x"),
           function(x, ...)
             standardGeneric("calculateJSD"))

#' @rdname deprecate
#' @export
setMethod("calculateJSD", signature = c(x = "SummarizedExperiment"),
          function(x, ...){
            .Deprecated(msg = paste0("'calculateJSD' is deprecated\n",
                                     "Use 'getDissimilarity' with parameter ",
                                     "method = 'jsd' instead."))
            getDissimilarity(x, method = "jsd", ...)
          }
          
)

#' @rdname deprecate
#' @export
setGeneric("runJSD", signature = c("x"),
           function(x, ...)
             standardGeneric("runJSD"))

#' @rdname deprecate
#' @export
setMethod("runJSD", signature = c(x = "SummarizedExperiment"),
          function(x, ...){
            .Deprecated(msg = paste0("'runJSD' is deprecated\n",
                                     "Use 'getDissimilarity' with parameter ",
                                     "method = 'jsd' instead."))
            getDissimilarity(x, method = "jsd", ...)
          }
)

#' @rdname deprecate
#' @export
setGeneric("calculateUnifrac", signature = c("x"),
           function(x, ...)
             standardGeneric("calculateUnifrac"))

#' @rdname deprecate
#' @export
setMethod("calculateUnifrac", signature = c(x = "TreeSummarizedExperiment"),
          function(x, ...){
            .Deprecated(msg = paste0("'calculateUnifrac' is deprecated\n",
                                     "Use 'getDissimilarity' with parameter ",
                                     "method = 'unifrac' instead."))
            getDissimilarity(x, method = "unifrac", ...)
          }
)

#' @rdname deprecate
#' @export
runUnifrac <- function(mat, tree, ...){
  .Deprecated(msg = paste0("'runUnifrac' is deprecated\n",
                           "Use 'getDissimilarity' with parameter ",
                           "method = 'unifrac' instead."))
  getDissimilarity(t(mat), tree = tree, method = "unifrac", ...)
}