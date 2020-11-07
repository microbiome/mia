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
#' @importFrom dplyr %>%
#' @importFrom rlang sym :=
NULL

#' @title MicrobiomeExperiment datasets
#'
#' @description
#' These datasets are conversions of the phyloseq datasets \code{GlobalPatterns}
#' \code{enterotype}, \code{esophagus}, \code{soilrep} and \code{MovingPictures}.
#'
#' The phyloseq version for MovingPictures data set was obtained from
#' \url{github.com/twbattaglia/MicrobeDS} before conversion.
#' 
#' \code{taxa} contains example taxonomic data derived from the
#'   \code{metagenomeFeatures} package.
#'
#' @name MicrobiomeExperiment-datasets
#' @docType data
#' @keywords datasets
#' @usage data(GlobalPatterns)
"GlobalPatterns"
#' @name MicrobiomeExperiment-datasets
#' @usage data(enterotype)
"enterotype"
#' @name MicrobiomeExperiment-datasets
#' @usage data(esophagus)
"esophagus"
#' @name MicrobiomeExperiment-datasets
#' @usage data(soilrep)
"soilrep"
#' @name MicrobiomeExperiment-datasets
#' @usage data(MovingPictures)
"MovingPictures"
