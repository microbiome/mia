#' Defunct functions in mia
#'
#' These functions were deprecated and are now defunct. Calling them gives an
#' error that names the replacement.
#'
#' @details
#' The functions and their replacements:
#' \itemize{
#'   \item \code{addPerSampleDominantFeatures}, \code{addPerSampleDominantTaxa}:
#'     \code{\link{addDominant}}
#'   \item \code{addTaxonomyTree}: \code{\link{addHierarchyTree}}
#'   \item \code{bestDMNFit}, \code{getBestDMNFit}, \code{getDMN},
#'     \code{runDMN}, and \code{calculateDMN} for a
#'     \code{SummarizedExperiment}: \code{\link{addCluster}} with
#'     \code{bluster::DmmParam()}
#'   \item \code{calculateJSD}, \code{runJSD}: \code{\link{getDissimilarity}}
#'     with \code{method = "jsd"}
#'   \item \code{calculateOverlap}: \code{\link{getDissimilarity}} with
#'     \code{method = "overlap"}
#'   \item \code{calculateUnifrac}, \code{runUnifrac}:
#'     \code{\link{getDissimilarity}} with \code{method = "unifrac"}
#'   \item \code{cluster}: \code{\link{addCluster}}
#'   \item \code{countDominantFeatures}, \code{countDominantTaxa}:
#'     \code{\link{summarizeDominance}}
#'   \item \code{estimateDiversity}, \code{estimateDominance},
#'     \code{estimateEvenness}, \code{estimateFaith}, \code{estimateRichness}:
#'     \code{\link{addAlpha}}
#'   \item \code{estimateDivergence}: \code{\link{addDivergence}}
#'   \item \code{full_join}, \code{inner_join}, \code{left_join},
#'     \code{right_join}: \code{\link{mergeSEs}} with \code{join}
#'   \item \code{getExperimentCrossAssociation},
#'     \code{getExperimentCrossCorrelation},
#'     \code{testExperimentCrossAssociation},
#'     \code{testExperimentCrossCorrelation}: \code{\link{getCrossAssociation}}
#'   \item \code{getPrevalentFeatures}, \code{getPrevalentTaxa}:
#'     \code{\link{getPrevalent}}
#'   \item \code{getRareFeatures}, \code{getRareTaxa}: \code{\link{getRare}}
#'   \item \code{getTopFeatures}, \code{getTopTaxa}: \code{\link{getTop}}
#'   \item \code{getUniqueFeatures}, \code{getUniqueTaxa}:
#'     \code{\link{getUnique}}
#'   \item \code{loadFromBiom}: \code{\link{importBIOM}}
#'   \item \code{loadFromHumann}: \code{\link{importHUMAnN}}
#'   \item \code{loadFromMetaphlan}: \code{\link{importMetaPhlAn}}
#'   \item \code{loadFromMothur}: \code{\link{importMothur}}
#'   \item \code{loadFromQIIME2}: \code{\link{importQIIME2}}
#'   \item \code{makePhyloseqFromTreeSE},
#'     \code{makePhyloseqFromTreeSummarizedExperiment}:
#'     \code{\link{convertToPhyloseq}}
#'   \item \code{makeTreeSEFromBiom},
#'     \code{makeTreeSummarizedExperimentFromBiom}:
#'     \code{\link{convertFromBIOM}}
#'   \item \code{makeTreeSEFromDADA2},
#'     \code{makeTreeSummarizedExperimentFromDADA2}:
#'     \code{\link{convertFromDADA2}}
#'   \item \code{makeTreeSEFromPhyloseq},
#'     \code{makeTreeSummarizedExperimentFromPhyloseq}:
#'     \code{\link{convertFromPhyloseq}}
#'   \item \code{meltAssay}: \code{\link{meltSE}}
#'   \item \code{mergeCols}, \code{mergeSamples}:
#'     \code{\link{agglomerateByVariable}} with \code{by = "cols"}
#'   \item \code{mergeFeatures}, \code{mergeRows}:
#'     \code{\link{agglomerateByVariable}} with \code{by = "rows"}
#'   \item \code{mergeFeaturesByPrevalence}: \code{\link{subsetBy}} with
#'     \code{collapse = TRUE}
#'   \item \code{mergeFeaturesByRank}: \code{\link{agglomerateByRank}}
#'   \item \code{perSampleDominantFeatures}, \code{perSampleDominantTaxa}:
#'     \code{\link{getDominant}}
#'   \item \code{plotNMDS}: \code{scater::plotReducedDim} with
#'     \code{dimred = "NMDS"}
#'   \item \code{readQZA}: \code{\link{importQZA}}
#'   \item \code{relabundance}, \code{relabundance<-}:
#'     \code{assay(x, "relabundance")}
#'   \item \code{runOverlap}: \code{\link{addDissimilarity}} with
#'     \code{method = "overlap"}
#'   \item \code{subsampleCounts}: \code{\link{rarefyAssay}}
#'   \item \code{subsetByPrevalentFeatures}, \code{subsetByPrevalentTaxa},
#'     \code{subsetByRareFeatures}, \code{subsetByRareTaxa}:
#'     \code{\link{subsetBy}}
#'   \item \code{subsetFeatures}, \code{subsetSamples}, \code{subsetTaxa}:
#'     \code{[} for subsetting
#'   \item \code{transformCounts}, \code{transformFeatures},
#'     \code{transformSamples}, \code{relAbundanceCounts}, \code{ZTransform}:
#'     \code{\link{transformAssay}}
#' }
#'
#' @param x Ignored.
#'
#' @param value Ignored.
#'
#' @param ... Ignored.
#'
#' @return These functions do not return; they give an error.
#'
#' @name mia-defunct
NULL

# Returns a defunct stub. 'new' is the replacement as shown in the message.
.defunct <- function(old, new){
    msg <- paste0("'", old, "' is defunct. Use ", new, " instead.")
    function(...) .Defunct(msg = msg)
}

#' @rdname mia-defunct
#' @export
addPerSampleDominantFeatures <- .defunct(
    "addPerSampleDominantFeatures", "'addDominant'")

#' @rdname mia-defunct
#' @export
addPerSampleDominantTaxa <- .defunct(
    "addPerSampleDominantTaxa", "'addDominant'")

#' @rdname mia-defunct
#' @export
addTaxonomyTree <- .defunct("addTaxonomyTree", "'addHierarchyTree'")

#' @rdname mia-defunct
#' @export
bestDMNFit <- .defunct(
    "bestDMNFit", "'addCluster' with 'bluster::DmmParam()' and full = TRUE")

#' @rdname mia-defunct
#' @export
calculateJSD <- .defunct(
    "calculateJSD", "'getDissimilarity' with method = 'jsd'")

#' @rdname mia-defunct
#' @export
calculateOverlap <- .defunct(
    "calculateOverlap", "'getDissimilarity' with method = 'overlap'")

#' @rdname mia-defunct
#' @export
calculateUnifrac <- .defunct(
    "calculateUnifrac", "'getDissimilarity' with method = 'unifrac'")

#' @rdname mia-defunct
#' @export
cluster <- .defunct("cluster", "'addCluster'")

#' @rdname mia-defunct
#' @export
countDominantFeatures <- .defunct(
    "countDominantFeatures", "'summarizeDominance'")

#' @rdname mia-defunct
#' @export
countDominantTaxa <- .defunct("countDominantTaxa", "'summarizeDominance'")

#' @rdname mia-defunct
#' @export
estimateDivergence <- .defunct("estimateDivergence", "'addDivergence'")

#' @rdname mia-defunct
#' @export
estimateDiversity <- .defunct("estimateDiversity", "'addAlpha'")

#' @rdname mia-defunct
#' @export
estimateDominance <- .defunct("estimateDominance", "'addAlpha'")

#' @rdname mia-defunct
#' @export
estimateEvenness <- .defunct("estimateEvenness", "'addAlpha'")

#' @rdname mia-defunct
#' @export
estimateFaith <- .defunct("estimateFaith", "'addAlpha'")

#' @rdname mia-defunct
#' @export
estimateRichness <- .defunct("estimateRichness", "'addAlpha'")

#' @rdname mia-defunct
#' @export
full_join <- .defunct("full_join", "'mergeSEs' with join = 'full'")

#' @rdname mia-defunct
#' @export
getBestDMNFit <- .defunct(
    "getBestDMNFit", "'addCluster' with 'bluster::DmmParam()' and full = TRUE")

#' @rdname mia-defunct
#' @export
getDMN <- .defunct(
    "getDMN", "'addCluster' with 'bluster::DmmParam()' and full = TRUE")

#' @rdname mia-defunct
#' @export
getExperimentCrossAssociation <- .defunct(
    "getExperimentCrossAssociation", "'getCrossAssociation'")

#' @rdname mia-defunct
#' @export
getExperimentCrossCorrelation <- .defunct(
    "getExperimentCrossCorrelation", "'getCrossAssociation'")

#' @rdname mia-defunct
#' @export
getPrevalentFeatures <- .defunct("getPrevalentFeatures", "'getPrevalent'")

#' @rdname mia-defunct
#' @export
getPrevalentTaxa <- .defunct("getPrevalentTaxa", "'getPrevalent'")

#' @rdname mia-defunct
#' @export
getRareFeatures <- .defunct("getRareFeatures", "'getRare'")

#' @rdname mia-defunct
#' @export
getRareTaxa <- .defunct("getRareTaxa", "'getRare'")

#' @rdname mia-defunct
#' @export
getTopFeatures <- .defunct("getTopFeatures", "'getTop'")

#' @rdname mia-defunct
#' @export
getTopTaxa <- .defunct("getTopTaxa", "'getTop'")

#' @rdname mia-defunct
#' @export
getUniqueFeatures <- .defunct("getUniqueFeatures", "'getUnique'")

#' @rdname mia-defunct
#' @export
getUniqueTaxa <- .defunct("getUniqueTaxa", "'getUnique'")

#' @rdname mia-defunct
#' @export
inner_join <- .defunct("inner_join", "'mergeSEs' with join = 'inner'")

#' @rdname mia-defunct
#' @export
left_join <- .defunct("left_join", "'mergeSEs' with join = 'left'")

#' @rdname mia-defunct
#' @export
loadFromBiom <- .defunct("loadFromBiom", "'importBIOM'")

#' @rdname mia-defunct
#' @export
loadFromHumann <- .defunct("loadFromHumann", "'importHUMAnN'")

#' @rdname mia-defunct
#' @export
loadFromMetaphlan <- .defunct("loadFromMetaphlan", "'importMetaPhlAn'")

#' @rdname mia-defunct
#' @export
loadFromMothur <- .defunct("loadFromMothur", "'importMothur'")

#' @rdname mia-defunct
#' @export
loadFromQIIME2 <- .defunct("loadFromQIIME2", "'importQIIME2'")

#' @rdname mia-defunct
#' @export
makePhyloseqFromTreeSE <- .defunct(
    "makePhyloseqFromTreeSE", "'convertToPhyloseq'")

#' @rdname mia-defunct
#' @export
makePhyloseqFromTreeSummarizedExperiment <- .defunct(
    "makePhyloseqFromTreeSummarizedExperiment", "'convertToPhyloseq'")

#' @rdname mia-defunct
#' @export
makeTreeSEFromBiom <- .defunct("makeTreeSEFromBiom", "'convertFromBIOM'")

#' @rdname mia-defunct
#' @export
makeTreeSEFromDADA2 <- .defunct("makeTreeSEFromDADA2", "'convertFromDADA2'")

#' @rdname mia-defunct
#' @export
makeTreeSEFromPhyloseq <- .defunct(
    "makeTreeSEFromPhyloseq", "'convertFromPhyloseq'")

#' @rdname mia-defunct
#' @export
makeTreeSummarizedExperimentFromBiom <- .defunct(
    "makeTreeSummarizedExperimentFromBiom", "'convertFromBIOM'")

#' @rdname mia-defunct
#' @export
makeTreeSummarizedExperimentFromDADA2 <- .defunct(
    "makeTreeSummarizedExperimentFromDADA2", "'convertFromDADA2'")

#' @rdname mia-defunct
#' @export
makeTreeSummarizedExperimentFromPhyloseq <- .defunct(
    "makeTreeSummarizedExperimentFromPhyloseq", "'convertFromPhyloseq'")

#' @rdname mia-defunct
#' @export
meltAssay <- .defunct("meltAssay", "'meltSE'")

#' @rdname mia-defunct
#' @export
mergeCols <- .defunct("mergeCols", "'agglomerateByVariable' with by = 'cols'")

#' @rdname mia-defunct
#' @export
mergeFeatures <- .defunct(
    "mergeFeatures", "'agglomerateByVariable' with by = 'rows'")

#' @rdname mia-defunct
#' @export
mergeFeaturesByPrevalence <- .defunct(
    "mergeFeaturesByPrevalence", "'subsetBy' with collapse = TRUE")

#' @rdname mia-defunct
#' @export
mergeFeaturesByRank <- .defunct("mergeFeaturesByRank", "'agglomerateByRank'")

#' @rdname mia-defunct
#' @export
mergeRows <- .defunct("mergeRows", "'agglomerateByVariable' with by = 'rows'")

#' @rdname mia-defunct
#' @export
mergeSamples <- .defunct(
    "mergeSamples", "'agglomerateByVariable' with by = 'cols'")

#' @rdname mia-defunct
#' @export
perSampleDominantFeatures <- .defunct(
    "perSampleDominantFeatures", "'getDominant'")

#' @rdname mia-defunct
#' @export
perSampleDominantTaxa <- .defunct("perSampleDominantTaxa", "'getDominant'")

#' @rdname mia-defunct
#' @export
plotNMDS <- .defunct(
    "plotNMDS", "'scater::plotReducedDim' with dimred = 'NMDS'")

#' @rdname mia-defunct
#' @export
readQZA <- .defunct("readQZA", "'importQZA'")

#' @rdname mia-defunct
#' @export
relabundance <- .defunct("relabundance", "assay(x, 'relabundance')")

#' @rdname mia-defunct
#' @export
`relabundance<-` <- function(x, value) relabundance()

#' @rdname mia-defunct
#' @export
relAbundanceCounts <- .defunct("relAbundanceCounts", "'transformAssay'")

#' @rdname mia-defunct
#' @export
right_join <- .defunct("right_join", "'mergeSEs' with join = 'right'")

#' @rdname mia-defunct
#' @export
runDMN <- .defunct("runDMN", "'addCluster' with 'bluster::DmmParam()'")

#' @rdname mia-defunct
#' @export
runJSD <- .defunct("runJSD", "'getDissimilarity' with method = 'jsd'")

#' @rdname mia-defunct
#' @export
runOverlap <- .defunct(
    "runOverlap", "'addDissimilarity' with method = 'overlap'")

#' @rdname mia-defunct
#' @export
runUnifrac <- .defunct(
    "runUnifrac", "'getDissimilarity' with method = 'unifrac'")

#' @rdname mia-defunct
#' @export
subsampleCounts <- .defunct("subsampleCounts", "'rarefyAssay'")

#' @rdname mia-defunct
#' @export
subsetByPrevalentFeatures <- .defunct(
    "subsetByPrevalentFeatures", "'subsetBy'")

#' @rdname mia-defunct
#' @export
subsetByPrevalentTaxa <- .defunct("subsetByPrevalentTaxa", "'subsetBy'")

#' @rdname mia-defunct
#' @export
subsetByRareFeatures <- .defunct("subsetByRareFeatures", "'subsetBy'")

#' @rdname mia-defunct
#' @export
subsetByRareTaxa <- .defunct("subsetByRareTaxa", "'subsetBy'")

#' @rdname mia-defunct
#' @export
subsetFeatures <- .defunct("subsetFeatures", "'[]' for subsetting")

#' @rdname mia-defunct
#' @export
subsetSamples <- .defunct("subsetSamples", "'[]' for subsetting")

#' @rdname mia-defunct
#' @export
subsetTaxa <- .defunct("subsetTaxa", "'[]' for subsetting")

#' @rdname mia-defunct
#' @export
taxonomyTree <- .defunct("taxonomyTree", "'getHierarchyTree'")

#' @rdname mia-defunct
#' @export
testExperimentCrossAssociation <- .defunct(
    "testExperimentCrossAssociation", "'getCrossAssociation'")

#' @rdname mia-defunct
#' @export
testExperimentCrossCorrelation <- .defunct(
    "testExperimentCrossCorrelation", "'getCrossAssociation'")

#' @rdname mia-defunct
#' @export
transformCounts <- .defunct("transformCounts", "'transformAssay'")

#' @rdname mia-defunct
#' @export
transformFeatures <- .defunct("transformFeatures", "'transformAssay'")

#' @rdname mia-defunct
#' @export
transformSamples <- .defunct("transformSamples", "'transformAssay'")

#' @rdname mia-defunct
#' @export
ZTransform <- .defunct("ZTransform", "'transformAssay'")
