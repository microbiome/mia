#' mia Package.
#'
#' \code{mia} implements tools for microbiome analysis
#'
#' @name mia-package
#' @docType package
#' @seealso \link[MicrobiomeExperiment:MicrobiomeExperiment-class]{MicrobiomeExperiment} class
NULL

#' @import methods
#' @import MicrobiomeExperiment
#' @import TreeSummarizedExperiment
#' @import scater
#' @importFrom dplyr %>%
#' @importFrom rlang sym :=
NULL

#' @title mia datasets
#'
#' @description
#' These datasets are conversions of the phyloseq datasets \code{GlobalPatterns}
#' \code{enterotype}, \code{esophagus} and \code{soilrep}.
#'
#' \code{taxa} contains example taxonomic data derived from the
#'   \code{metagenomeFeatures} package.
#'
#' @name mia-datasets
#' @docType data
#' @keywords datasets
#' @usage data(GlobalPatterns)
"GlobalPatterns"
#' @name mia-datasets
#' @usage data(enterotype)
"enterotype"
#' @name mia-datasets
#' @usage data(esophagus)
"esophagus"
#' @name mia-datasets
#' @usage data(soilrep)
"soilrep"
