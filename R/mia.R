#' \code{mia} Package.
#'
#' \code{mia} implements tools for microbiome analysis based on the
#' SummarizedExperiment, SingleCellExperiment and TreeSummarizedExperiment
#' infrastructure. Data wrangling and analysis are the main scope of this
#' package.
#'
#' @name mia-package
#' @docType package
#' @seealso \link[TreeSummarizedExperiment:TreeSummarizedExperiment-class]{TreeSummarizedExperiment} class
NULL

#' @import methods
#' @import TreeSummarizedExperiment
#' @import DelayedArray
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
#' \code{dmn_se} contains an example \code{SummarizedExperiment} derived
#' from data in the \code{DirichletMultinomal} package. See \code{?calculateDMN}
#' for more details.
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
#' @name mia-datasets
#' @usage data(dmn_se)
"dmn_se"
