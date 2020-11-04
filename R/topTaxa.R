#' Returns top N \dQuote{FeatureID}
#'
#' Basic extractor function for extracting top N \dQuote{FeatureID}
#' from assay for a user specified \dQuote{SampleID} in a
#' \code{\link[=SummarizedExperiment-class]{SummarizedExperiment}}
#' object.
#'
#' @param x A
#'  \code{\link[=SummarizedExperiment-class]{SummarizedExperiment}} object.
#'
#' @param top Numeric value, how many top taxa to return. Default return top five taxa.
#'
#' @param method Specify the method to determine top taxa. Either sum, mean or median.
#'               Default is 'mean'
#'
#' @param abund_values a \code{character} value to select an
#'   \code{\link[SummarizedExperiment:SummarizedExperiment-class]{assayNames}}
#'
#' @importFrom DelayedMatrixStats rowSums2 rowMeans2 rowMedians
#'
#' @return A vector of the top N \dQuote{FeatureID}
#' @name topTaxa
#'
#' @author
#' Sudarshan A. Shetty
#'
#' @examples
#' data(GlobalPatterns, package = "MicrobiomeExperiment")
#'
#' top_taxa <- topTaxa(GlobalPatterns,
#'                     method="mean",
#'                     top=5,
#'                     abund_values="counts")
#'
NULL

#' @rdname topTaxa
#'
#' @export

setGeneric("topTaxa", signature = "x",
           function(x, top, method, abund_values)
               standardGeneric("topTaxa"))

.check_max_taxa <- function(x, top, abund_values){

    if(top > nrow(assay(x,abund_values))){
        stop("top value is greater than maximum number of taxa in assay", call. = FALSE)
    }

}

#' @rdname topTaxa
#'
#' @export
#'
setMethod("topTaxa", signature = c(x = "SummarizedExperiment"),
          function(x, top=5, method="mean", abund_values = "counts"){

              # check max taxa
              .check_max_taxa(x, top, abund_values)

              # check assay
              .check_abund_values(abund_values, x)
              #assay(x, abund_values)[,sample_id]

              if(method=="sum"){

                  taxs <- rowSums2(assay(x, abund_values))

              } else if (method=="mean"){

                  taxs <- rowMeans2(assay(x, abund_values))

              } else if (method=="median") {

                  taxs <- rowMedians(assay(x, abund_values))

              }

              names(taxs) <- rownames(assay(x))
              taxs <- sort(taxs,decreasing = TRUE)[1:top]
              return(names(taxs))
          }
)





