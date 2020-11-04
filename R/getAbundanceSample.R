#' Returns abundance values for all \dQuote{FeatureID} in a \dQuote{SampleID}
#'
#' Basic extractor function for extracting \dQuote{FeatureID} abundances
#' from assay for a user specified \dQuote{SampleID} in a
#' \code{\link[=SummarizedExperiment-class]{SummarizedExperiment}}
#' object.
#'
#' @param x A
#'  \code{\link[=SummarizedExperiment-class]{SummarizedExperiment}} object.
#'
#' @param sample_id A \dQuote{SampleID} from which user wants to extract
#'  the abundances of \dQuote{FeatureID}
#'
#' @param abund_values a \code{character} value to select an
#'   \code{\link[SummarizedExperiment:SummarizedExperiment-class]{assayNames}}
#'
#' @return An integer vector of the abundance values for
#' all \dQuote{FeatureID} in
#' \code{\link[SummarizedExperiment:SummarizedExperiment-class]{assayNames}}
#' for specified \dQuote{SampleID}
#'
#' @name getAbundanceSample
#'
#' @author
#' Sudarshan A. Shetty
#'
#' @examples
#' data(GlobalPatterns, package = "MicrobiomeExperiment")
#' getAbundanceSample(GlobalPatterns,
#'                    sample_id = "CC1",
#'                    abund_values = "counts")
NULL

#' @rdname getAbundanceSample
#'
#' @export

setGeneric("getAbundanceSample", signature = "x",
           function(x, sample_id, abund_values)
               standardGeneric("getAbundanceSample"))



#' @rdname getAbundanceSample
#'
#' @export
#'
setMethod("getAbundanceSample", signature = c(x = "SummarizedExperiment"),
          function(x, sample_id = NULL, abund_values = "counts"){

              # check assay
              .check_abund_values(abund_values, x)
              # check if sampleid exists or matches
              .check_sample_ids_assays(x, sample_id, abund_values)

              assay(x, abund_values)[,sample_id]

}
)

#' @importFrom SummarizedExperiment assay
#'

.check_sample_ids_assays <- function(x, sample_id, abund_values){

    if(is.null(sample_id)) {
        stop("'sample_id' must be a single non-empty character value",
             call. = FALSE)
    } else if (isFALSE(any(sample_id %in% colnames(assay(x, abund_values))))) {

        stop("Please provide a valid 'sample_id'", call. = FALSE)
    }

}


