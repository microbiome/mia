#' Get abundance values by \dQuote{SampleID} or \dQuote{FeatureID}
#'
#' These are functions for extracting abundances present in \code{assay(x)}.
#' These functions are convenience wrapper around subsetting columns or rows
#' from \code{assay(x,name)}.
#'
#' \code{getAbundanceSample} returns abundance values for all
#' \dQuote{FeatureIDs} in a user specified \dQuote{SampleID}.
#'
#' \code{getAbundanceFeature} returns abundance values in all \dQuote{SampleIDs}
#' for user specified \dQuote{FeatureID}.
#'
#' @param x A
#'  \code{\link[=SummarizedExperiment-class]{SummarizedExperiment}} object.
#'
#' @param sample_id A \dQuote{SampleID} from which user wants to extract the
#'   abundances of \dQuote{FeatureID}. This is essentially a column name in
#'   \code{assay(x)}.
#'
#' @param feature_id A \dQuote{FeatureID} for which user wants to extract the
#'   abundances from all of \dQuote{SampleID} in
#'   \code{\link[SummarizedExperiment:SummarizedExperiment-class]{assayNames}}.
#'   This is essentially a rowname in \code{assay(x)}.
#'
#' @param abund_values a \code{character} value to select an
#'   \code{\link[SummarizedExperiment:SummarizedExperiment-class]{assayNames}}
#'
#' @return \code{getAbundanceSample} and \code{getAbundanceFeature} return a
#'   numeric matrix of the abundance values for all
#'   \dQuote{SampleIDs}/\dQuote{FeatureIDs}
#'
#' @name getAbundance
#'
#' @author
#' Sudarshan A. Shetty
#'
#' @examples
#' # getAbundanceSample
#' data(GlobalPatterns)
#' getAbundanceSample(GlobalPatterns,
#'                    sample_id = 'CC1',
#'                    abund_values = 'counts')
#' # getAbundanceFeature
#' getAbundanceFeature(GlobalPatterns,
#'                     feature_id = '522457',
#'                     abund_values = 'counts')
NULL

#' @rdname getAbundance
#' @export
setGeneric("getAbundanceSample", signature = "x",
           function(x, sample_id, abund_values = "counts")
               standardGeneric("getAbundanceSample"))

#' @rdname getAbundance
#' @export
setMethod("getAbundanceSample", signature = c(x = "SummarizedExperiment"),
    function(x, sample_id = NULL, abund_values = "counts") {
        # check assay
        .check_abund_values(abund_values, x)
        # check if sampleid exists or matches
        .check_feature_sample_ids(sample_id, colnames(x))
        assay <- assay(x, abund_values)
        rowSums(assay[,colnames(assay) %in% sample_id,drop=FALSE])
    }
)

#' @importFrom SummarizedExperiment assay
.check_feature_sample_ids <- function(id, names,
                                      id_name = .get_name_in_parent(id)) {
    if (is.null(id) || !is.character(id) || length(id) > 1L) {
        stop("'",id_name,"' must be a single non-empty character value",
             call. = FALSE)
    } else if (!(id %in% names)) {
        stop("Please provide a valid '",id_name,"'", call. = FALSE)
    }
}

#' @rdname getAbundance
#' @export
setGeneric("getAbundanceFeature", signature = "x",
           function(x, feature_id, abund_values)
               standardGeneric("getAbundanceFeature"))

#' @rdname getAbundance
#' @export
setMethod("getAbundanceFeature", signature = c(x = "SummarizedExperiment"),
    function(x, feature_id = NULL, abund_values = "counts") {
        # check assay
        .check_abund_values(abund_values, x)
        # check if feature_id exists or matches
        .check_feature_sample_ids(feature_id, rownames(x))
        assay <- assay(x, abund_values)
        colSums(assay[rownames(assay) %in% feature_id,,drop=FALSE])
    }
)
