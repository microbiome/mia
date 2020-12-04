#' Prevalent abundance
#'
#' This function calculates the joint abundance of the prevalent taxa.
#'
#' @param x a
#'   \code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#'   object
#'
#' @param abund_values A single character value for selecting the
#'   \code{\link[SummarizedExperiment:SummarizedExperiment-class]{assay}}
#'   to use for prevalence calculation.
#'
#' @inheritParams getPrevalentTaxa
#'
#' @details The core abundance index gives the relative proportion of the core
#' species (in between 0 and 1). The core taxa are defined as those that exceed the
#' given population prevalence threshold at the given detection level.
#'
#' @return
#' a named \code{numeric} vector. It includes joint abundance of prevalent taxa,
#' that is calculated for individual samples.
#'
#' @examples
#' data(esophagus, package = "MicrobiomeExperiment")
#' esophagus <- as(esophagus, "MicrobiomeExperiment")
#' prevalentAbundance(esophagus)
#'
#' @name prevalentAbundance
#' @export
#'
#' @author Leo Lahti and Tuomas Borman. Contact: \url{microbiome.github.io}
#'
#' @export
#'
setGeneric("prevalentAbundance", signature = "x",
           function(x, abund_values = "relabundance", ...)
               standardGeneric("prevalentAbundance"))

#' @rdname getPrevalence
#' @export
setMethod("prevalentAbundance", signature = c(x = "MicrobiomeExperiment"),
          function(x, abund_values = "relabundance", detection = 0.1/100, prevalence = 50/100,
                   include_lowest = FALSE, as_relative=FALSE, ...){

              #Adds relative abundance table to the data
              x <- relAbundanceCounts(x)

              #Saves the relative abundances (or counts, if wanted) to the variable "values"
              values <- assays(x)[[abund_values]]

              # Core members
              cm <- getPrevalentTaxa(x, detection=detection, prevalence=prevalence, include_lowest=include_lowest,
                                     as_relative=as_relative, rank=NULL)

              if (length(cm) == 0) {
                  warning("With the given abundance and prevalence
            thresholds, no taxa were found. Returning NA for prevalentAbundance. Try to
            change detection and prevalence parameters.");
                  ret <- NA
              }

              # Pick the core and calculate abundance
              if (is.vector(values)) {
                  ret <- values # 1 OTU x N samples
              } else if (length(cm) == 1) {
                  ret <- values[cm,]
              } else if (ncol(values) > 1 && length(cm) > 1) {
                  ret <- colSums(values[cm, ])
              } else if (ncol(values) == 1) {
                  ret <- sum(values[cm, ])
              }
              return(ret)


          }
)
