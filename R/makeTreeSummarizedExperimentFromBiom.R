#' Loading a biom file
#'
#' For convenience a few functions are available to convert data from a
#' \sQuote{biom} file or object into a
#' \code{\link[TreeSummarizedExperiment:TreeSummarizedExperiment-class]{TreeSummarizedExperiment}}
#'
#' @param file biom file location
#'
#' @return An object of class
#'   \code{\link[TreeSummarizedExperiment:TreeSummarizedExperiment-class]{TreeSummarizedExperiment}}
#'
#' @name makeTreeSEFromBiom
#' @seealso
#' \code{\link[=makeTreeSummarizedExperimentFromPhyloseq]{makeTreeSummarizedExperimentFromPhyloseq}}
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
#'   se <- loadFromBiom(rich_dense_file)
#'
#'   # load from object
#'   x1 <- biomformat::read_biom(rich_dense_file)
#'   se <- makeTreeSEFromBiom(x1)
#'   # Convert SE to TreeSE
#'   tse <- as(se, "TreeSummarizedExperiment")
#'   tse
#' }
NULL

#' @rdname makeTreeSEFromBiom
#'
#' @export
loadFromBiom <- function(file) {
    .require_package("biomformat")
    biom <- biomformat::read_biom(file)
    makeTreeSEFromBiom(biom)
}

#' @rdname makeTreeSEFromBiom
#'
#' @param obj object of type \code{\link[biomformat:read_biom]{biom}}
#'
#' @export
#' @importFrom S4Vectors make_zero_col_DFrame
makeTreeSEFromBiom <- function(obj){
    # input check
    .require_package("biomformat")
    if(!is(obj,"biom")){
        stop("'obj' must be a 'biom' object")
    }
    #
    counts <- as(biomformat::biom_data(obj), "matrix")
    sample_data <- biomformat::sample_metadata(obj)
    feature_data <- biomformat::observation_metadata(obj)
    
    # colData and rowData are initialized with empty tables with rownames if they are NULL
    if( is.null(sample_data) ){
        sample_data <- S4Vectors::make_zero_col_DFrame(ncol(counts))
        rownames(sample_data) <- colnames(counts)
    }
    if( is.null(feature_data) ){
        feature_data <- S4Vectors::make_zero_col_DFrame(nrow(counts))
        rownames(feature_data) <- rownames(counts)
    }
    
    TreeSummarizedExperiment(assays = list(counts = counts),
                            colData = sample_data,
                            rowData = feature_data)
}
