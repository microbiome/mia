#' @title Get abundance values by \dQuote{SampleID} or \dQuote{FeatureID}
#'
#' @description
#' These are basic functions for extracting abundances present in \code{assay(x)}.
#' These functions are convenience wrapper around subsetting columns or rows from
#' assay(x,name).
#'
#' \code{getAbundanceSample} returns abundance values for all \dQuote{FeatureIDs} in a
#' user specified \dQuote{SampleID}.
#'
#' \code{getAbundanceFeature} returns abundance values in all \dQuote{SampleIDs} for
#' user specified \dQuote{FeatureID}.
#'
#' @param x A
#'  \code{\link[=SummarizedExperiment-class]{SummarizedExperiment}} object.
#'
#' @param sample_id A \dQuote{SampleID} from which user wants to extract
#'  the abundances of \dQuote{FeatureID}. This is essentially a colname in assay(x).
#'
#' @param abund_values a \code{character} value to select an
#'   \code{\link[SummarizedExperiment:SummarizedExperiment-class]{assayNames}}
#'
#' @return \code{getAbundanceSample} An integer vector of the abundance values
#'         for all \dQuote{FeatureIDs} in user-specified \dQuote{SampleID}
#'
#' @name getAbundance
#'
#' @author
#' Sudarshan A. Shetty
#'
#' @examples
#'
#' # getAbundanceSample
#' data(GlobalPatterns, package = 'MicrobiomeExperiment')
#' getAbundanceSample(GlobalPatterns,
#'                    sample_id = 'CC1',
#'                    abund_values = 'counts')
NULL

#' @rdname getAbundance
#'
#' @export

setGeneric("getAbundanceSample", signature = "x", function(x, sample_id,
                                                           abund_values) standardGeneric("getAbundanceSample"))


#' @rdname getAbundance
#'
#' @aliases getAbundanceSample
#'
#' @export
#'
setMethod("getAbundanceSample", signature = c(x = "SummarizedExperiment"),
          function(x, sample_id = NULL, abund_values = "counts") {
              # check assay
              .check_abund_values(abund_values, x)
              # check if sampleid exists or matches
              .check_sample_ids_assays(x, sample_id, abund_values)
              assay(x, abund_values)[, sample_id]
          })


#' @importFrom SummarizedExperiment assay
#'

.check_sample_ids_assays <- function(x, sample_id, abund_values) {
    if (is.null(sample_id)) {
        stop("'sample_id' must be a single non-empty character value",
             call. = FALSE)
    } else if (isFALSE(any(sample_id %in% colnames(assay(x, abund_values))))) {
        stop("Please provide a valid 'sample_id'", call. = FALSE)
    }
}

#' @param feature_id A \dQuote{FeatureID} for which user wants to extract
#'  the abundances from all of \dQuote{SampleID} in
#'  \code{\link[SummarizedExperiment:SummarizedExperiment-class]{assayNames}}.
#'  This is essentially a rowname in assay(x).
#'
#' @return \code{getAbundanceFeature} An integer vector of the abundance values
#'         for \dQuote{FeatureID} in all \dQuote{SampleIDs}.
#'
#' @name getAbundance
#'
#' @examples
#'
#' # getAbundanceFeature
#' data(GlobalPatterns, package = 'MicrobiomeExperiment')
#' getAbundanceFeature(GlobalPatterns,
#'                     feature_id = '522457',
#'                     abund_values = 'counts')
NULL

#' @rdname getAbundance
#'
#' @export

setGeneric("getAbundanceFeature", signature = "x", function(x, feature_id,
                                                            abund_values) standardGeneric("getAbundanceFeature"))


#' @rdname getAbundance
#'
#' @aliases getAbundanceFeature
#'
#' @export
#'
setMethod("getAbundanceFeature", signature = c(x = "SummarizedExperiment"),
          function(x, feature_id = NULL, abund_values = "counts") {
              # check assay
              .check_abund_values(abund_values, x)
              # check if feature_id exists or matches
              .check_feature_ids_assays(x, feature_id, abund_values)
              assay(x, abund_values)[feature_id, ]
          })

#' @importFrom SummarizedExperiment assay
#'
.check_feature_ids_assays <- function(x, feature_id, abund_values) {
    if (is.null(feature_id)) {
        stop("'feature_id' must be a single non-empty character value",
             call. = FALSE)
    } else if (isFALSE(any(feature_id %in% rownames(assay(x, abund_values))))) {
        stop("'feature_id' must be in rownames(assay(x))", call. = FALSE)
    }
}
