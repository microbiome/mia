#' Loading a biom file
#'
#' For convenience a few functions are available to convert data from a
#' \sQuote{biom} file or object into a
#' \code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#'
#' @param file biom file location
#'
#' @return An object of class
#'   \code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#'
#' @name makeSummarizedExperimentFromBiom
#' @seealso
#' \code{\link[=makeTreeSummarizedExperimentFromphyloseq]{makeTreeSummarizedExperimentFromphyloseq}}
#' \code{\link[=makeTreeSummarizedExperimentFromDADA2]{makeTreeSummarizedExperimentFromDADA2}}
#' \code{\link[=loadFromQIIME2]{loadFromQIIME2}}
#' \code{\link[=loadFromMothur]{loadFromMothur}}
#'
#' @examples
#' if(requireNamespace("biomformat")) {
#'   library(biomformat)
#'   # load from file
#'   rich_dense_file  = system.file("extdata", "rich_dense_otu_table.biom",
#'                                  package = "biomformat")
#'   tse <- loadFromBiom(rich_dense_file)
#'
#'   # load from object
#'   x1 <- biomformat::read_biom(rich_dense_file)
#'   tse <- makeSummarizedExperimentFromBiom(x1)
#'   tse
#' }
NULL

#' @rdname makeSummarizedExperimentFromBiom
#'
#' @export
loadFromBiom <- function(file) {
    .require_package("biomformat")
    biom <- biomformat::read_biom(file)
    makeSummarizedExperimentFromBiom(biom)
}

#' @rdname makeSummarizedExperimentFromBiom
#'
#' @param obj object of type \code{\link[biomformat:read_biom]{biom}}
#'
#' @export
makeSummarizedExperimentFromBiom <- function(obj){
    # input check
    .require_package("biomformat")
    if(!is(obj,"biom")){
        stop("'obj' must be a 'biom' object")
    }
    #
    counts <- as(biomformat::biom_data(obj), "matrix")
    sample_data <- biomformat::sample_metadata(obj)
    feature_data <- biomformat::observation_metadata(obj)
    
    # If sample_data is not included in the file, it's NULL, which leads to error
    # when object is created. --> NULL is replaced with empty data frame.
    if( is.null(sample_data) ){
        sample_data <- S4Vectors:::make_zero_col_DataFrame(ncol(counts))
    }

    SummarizedExperiment(assays = list(counts = counts),
                         colData = sample_data,
                         rowData = feature_data)
}
