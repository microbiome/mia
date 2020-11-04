#' Returns abundance values for user specified \dQuote{FeatureID}
#' in all \dQuote{SampleID}
#'
#' Basic extractor function for extracting abundances for a user
#' specified \dQuote{FeatureID} from assay in a
#' \code{\link[=SummarizedExperiment-class]{SummarizedExperiment}}
#' object.
#'
#' @param x A
#'  \code{\link[=SummarizedExperiment-class]{SummarizedExperiment}} object.
#'
#' @param feature_id A \dQuote{FeatureID} for which user wants to extract
#'  the abundances from all of \dQuote{SampleID} in
#'  \code{\link[SummarizedExperiment:SummarizedExperiment-class]{assayNames}}
#'
#' @param abund_values a \code{character} value to select an
#'   \code{\link[SummarizedExperiment:SummarizedExperiment-class]{assayNames}}
#'
#' @return An integer vector of the abundance values for
#'  \dQuote{FeatureID} in all \dQuote{SampleID} present in
#'  \code{\link[SummarizedExperiment:SummarizedExperiment-class]{assayNames}}
#'
#' @name getAbundanceFeature
#'
#' @author
#' Sudarshan A. Shetty
#'
#' @examples
#' data(GlobalPatterns, package = "MicrobiomeExperiment")
#' getAbundanceFeature(GlobalPatterns,
#'                     feature_id = "522457",
#'                     abund_values = "counts")
NULL

#' @rdname getAbundanceFeature
#'
#' @export

setGeneric("getAbundanceFeature", signature = "x",
           function(x, feature_id, abund_values)
               standardGeneric("getAbundanceFeature"))


#' @rdname getAbundanceFeature
#'
#' @export
#'
setMethod("getAbundanceFeature", signature = c(x = "SummarizedExperiment"),
          function(x, feature_id = NULL, abund_values = "counts"){

              # check assay
              .check_abund_values(abund_values, x)
              # check if feature_id exists or matches
              .check_feature_ids_assays(x, feature_id, abund_values)

              assay(x, abund_values)[feature_id,]

          }
)


#' @importFrom SummarizedExperiment assay
#'
.check_feature_ids_assays <- function(x, feature_id, abund_values){


    if(is.null(feature_id)) {

        stop("'feature_id' must be a single non-empty character value",
             call. = FALSE)

    } else if (isFALSE(any(feature_id %in% rownames(assay(x, abund_values))))) {

        stop("Please provide a valid 'feature_id'", call. = FALSE)

    }

}
