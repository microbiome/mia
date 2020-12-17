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
#' @param ... additional parameters passed to \code{getPrevalentAbundance}
#'
#' @details The core abundance index gives the relative proportion of the core
#' species (in between 0 and 1). The core taxa are defined as those that exceed the
#' given population prevalence threshold at the given detection level.
#'
#' @return
#' a named \code{numeric} vector. It includes joint abundance of prevalent taxa,
#' that is calculated for individual samples.
#'
#' @name getPrevalentAbundance
#'
#' @author Leo Lahti and Tuomas Borman. Contact: \url{microbiome.github.io}
#'
#' @examples
#' data(esophagus)
#' getPrevalentAbundance(esophagus, abund_values = "counts")
NULL

#' @rdname getPrevalentAbundance
#' @export
setGeneric("getPrevalentAbundance", signature = "x",
           function(x, abund_values = "relabundance", ...)
               standardGeneric("getPrevalentAbundance"))

#' @rdname getPrevalentAbundance
#' @export
setMethod("getPrevalentAbundance", signature = c(x = "ANY"),
    function(x, ...){
        x <- .calc_rel_abund(x)
        cm <- getPrevalentTaxa(x, ...)
        if (length(cm) == 0) {
          stop("With the given abundance and prevalence thresholds, no taxa ",
               "were found. Try to change detection and prevalence parameters.",
               call. = FALSE)
        }
        colSums(x[cm, ,drop=FALSE])
    }
)

#' @rdname getPrevalence
#' @export
setMethod("getPrevalentAbundance", signature = c(x = "SummarizedExperiment"),
    function(x, abund_values = "counts", ...){
        # check assay
        .check_abund_values(abund_values, x)
        #
        getPrevalentAbundance(assay(x,abund_values))
    }
)
