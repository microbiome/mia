#' Return abundance values for all \dQuote{FeatureID} in a \dQuote{SampleID}
#'
#' This is a basic function for extracting abundances of all \dQuote{FeatureID}
#' for a user specified \dQuote{SampleID} from assay in a
#' \code{\link[=SummarizedExperiment-class]{SummarizedExperiment}} object.
#'
#' This function is aimed at simplifying basic/repetitive tasks and is a convenience
#' wrapper around subsetting columns from assay(x,name).
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


#' Return abundance values for user specified \dQuote{FeatureID}
#' in all \dQuote{SampleID}
#'
#' This is a basic function for extracting abundances of user specified \dQuote{FeatureID}
#' from assay in a \code{\link[=SummarizedExperiment-class]{SummarizedExperiment}} object.
#'
#' This function is aimed at simplifying basic/repetitive tasks and is a convenience
#' wrapper around subsetting rows from assay(x,name).
#'
#' @param x A
#'  \code{\link[=SummarizedExperiment-class]{SummarizedExperiment}} object.
#'
#' @param feature_id A \dQuote{FeatureID} for which user wants to extract
#'  the abundances from all of \dQuote{SampleID} in
#'  \code{\link[SummarizedExperiment:SummarizedExperiment-class]{assayNames}}.
#'  This is essentially a rowname in assay(x).
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
