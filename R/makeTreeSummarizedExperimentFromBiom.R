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
#' \code{\link[=makeTreeSEFromPhyloseq]{makeTreeSEFromPhyloseq}}
#' \code{\link[=makeTreeSEFromDADA2]{makeTreeSEFromDADA2}}
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
    
    # colData is initialized with empty tables with rownames if it is NULL
    if( is.null(sample_data) ){
        sample_data <- S4Vectors::make_zero_col_DFrame(ncol(counts))
        rownames(sample_data) <- colnames(counts)
    # Otherwise convert it into correct format
    } else{
        # Get the maximum length of list
        max_length <- max( lengths(sample_data) )
        # Get the column names from the taxa info that has all the columns that occurs
        # in the data
        colnames <- names( head( sample_data[ lengths(sample_data) == 
                                                  max_length ], 1)[[1]] )
        # Append the data with NAs if some samples do not have all the info
        sample_data <- lapply(sample_data, function(x){
            length(x) <- max_length 
            return(x)
        })
        # Create a data.frame from the list
        sample_data <- do.call(rbind, sample_data)
        # Add correct colnames
        colnames(sample_data) <- colnames
        # Convert into DataFrame
        sample_data <- DataFrame(sample_data)
    }
    # rowData is initialized with empty tables with rownames if it is NULL
    if( is.null(feature_data) ){
        feature_data <- S4Vectors::make_zero_col_DFrame(nrow(counts))
        rownames(feature_data) <- rownames(counts)
    # Otherwise convert it into correct format
    } else{
        # Feature data is a list of taxa info
        # Get the maximum length of list
        max_length <- max( lengths(feature_data) )
        # Get the column names from the taxa info that has all the levels that occurs
        # in the data
        colnames <- names( head( feature_data[ lengths(feature_data) == 
                                                   max_length ], 1)[[1]] )
        # Convert the list so that all invidual taxa info have the max length
        # of the list objects. All vectors are appended with NAs, if they do not
        # have all the levels. E.g., if only Kingdom level is found, all lower
        # ranks are now NA
        feature_data <- lapply(feature_data, function(x){
            length(x) <- max_length 
            return(x)
        })
        # Create a data.frame from the list
        feature_data <- do.call(rbind, feature_data)
        # Add correct colnames
        colnames(feature_data) <- colnames
        # Convert into DataFrame
        feature_data <- DataFrame(feature_data)
    }
    
    TreeSummarizedExperiment(assays = list(counts = counts),
                            colData = sample_data,
                            rowData = feature_data)
}

####################### makeTreeSummarizedExperimentFromBiom #######################
#' @param obj object of type \code{\link[biomformat:read_biom]{biom}}
#' @rdname makeTreeSEFromBiom
#' @export
makeTreeSummarizedExperimentFromBiom <- function(obj){
    makeTreeSEFromBiom(obj)
}
