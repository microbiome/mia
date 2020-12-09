#' Loading a biom file
#'
#' For convenciance a few functions are available to convert data from a
#' \sQuote{biom} file or object into a
#' \code{\link[MicrobiomeExperiment:MicrobiomeExperiment-class]{MicrobiomeExperiment}}
#'
#' @param file biom file location
#'
#' @return An object of class
#'   \code{\link[MicrobiomeExperiment:MicrobiomeExperiment-class]{MicrobiomeExperiment}}
#'
#' @name makeMicrobiomeExperimentFromBiom
#'
#' @examples
#' if(requireNamespace("biomformat")) {
#'   library(biomformat)
#'   # load from file
#'   rich_dense_file  = system.file("extdata", "rich_dense_otu_table.biom",
#'                                  package = "biomformat")
#'   me <- loadFromBiom(rich_dense_file)
#'
#'   # load from object
#'   x1 <- biomformat::read_biom(rich_dense_file)
#'   me <- makeMicrobiomeExperimentFromBiom(x1)
#'   me
#' }
NULL

#' @rdname makeMicrobiomeExperimentFromBiom
#'
#' @export
loadFromBiom <- function(file) {
    .require_package("biomformat")
    biom <- biomformat::read_biom(file)
    makeMicrobiomeExperimentFromBiom(biom)
}

#' @rdname makeMicrobiomeExperimentFromBiom
#'
#' @param obj object of type \code{\link[biomformat:read_biom]{biom}}
#'
#' @export
makeMicrobiomeExperimentFromBiom <- function(obj){
    # input check
    .require_package("biomformat")
    if(!is(obj,"biom")){
        stop("'obj' must be a 'biom' object")
    }
    #
    counts <- as(biomformat::biom_data(obj), "matrix")
    sample_data <- biomformat::sample_metadata(obj)
    feature_data <- biomformat::observation_metadata(obj)

    MicrobiomeExperiment(assays = list(counts = counts),
                         colData = sample_data,
                         rowData = feature_data)
}
