#' Loading a biom file
#'
#' For convenience a few functions are available to convert data from a
#' \sQuote{biom} file or object into a
#' \code{\link[TreeSummarizedExperiment:TreeSummarizedExperiment-class]{TreeSummarizedExperiment}}
#'
#' @param file biom file location
#' 
#' @param removeTaxaPrefixes \code{TRUE} or \code{FALSE}: Should
#' taxonomic prefixes be removed? (default \code{removeTaxaPrefixes = FALSE})
#' 
#' @param rankFromPrefix \code{TRUE} or \code{FALSE}: If file does not have
#' taxonomic ranks on feature table, should they be scraped from prefixes?
#' (default \code{rankFromPrefix = FALSE})
#' 
#' @param ... optional arguments (not used).
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
#'   se <- loadFromBiom(rich_dense_file, removeTaxaPrefixes = TRUE, rankFromPrefix = TRUE)
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
loadFromBiom <- function(file, ...) {
    .require_package("biomformat")
    biom <- biomformat::read_biom(file)
    makeTreeSEFromBiom(biom, ...)
}

#' @rdname makeTreeSEFromBiom
#'
#' @param obj object of type \code{\link[biomformat:read_biom]{biom}}
#'
#' @export
#' @importFrom S4Vectors make_zero_col_DFrame
makeTreeSEFromBiom <- function(
        obj, removeTaxaPrefixes = FALSE, rankFromPrefix = FALSE, ...){
    # input check
    .require_package("biomformat")
    if(!is(obj,"biom")){
        stop("'obj' must be a 'biom' object", call. = FALSE)
    }
    if( !.is_a_bool(removeTaxaPrefixes) ){
        stop("'removeTaxaPrefixes' must be TRUE or FALSE.", call. = FALSE)
    }
    if( !.is_a_bool(rankFromPrefix) ){
        stop("'rankFromPrefix' must be TRUE or FALSE.", call. = FALSE)
    }
    #
    counts <- as(biomformat::biom_data(obj), "matrix")
    sample_data <- biomformat::sample_metadata(obj)
    feature_data <- biomformat::observation_metadata(obj)
    
    # colData is initialized with empty tables with rownames if it is NULL
    if( is.null(sample_data) ){
        sample_data <- S4Vectors::make_zero_col_DFrame(ncol(counts))
        rownames(sample_data) <- colnames(counts)
    # Otherwise convert it into correct format if it is a list
    } else if( is(sample_data, "list") ){
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
    }
    # rowData is initialized with empty tables with rownames if it is NULL
    if( is.null(feature_data) ){
        feature_data <- S4Vectors::make_zero_col_DFrame(nrow(counts))
        rownames(feature_data) <- rownames(counts)
    # Otherwise convert it into correct format if it is a list
    } else if( is(feature_data, "list") ){
        # Feature data is a list of taxa info
        # Get the maximum length of list
        max_length <- max( lengths(feature_data) )
        # Get the column names from the taxa info that has all the levels that occurs
        # in the data
        colnames <- names( head( feature_data[ lengths(feature_data) == 
                                                   max_length ], 1)[[1]] )
        # Convert the list so that all individual taxa info have the max length
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
    }
    
    # Replace taxonomy ranks with ranks found based on prefixes
    if( rankFromPrefix && all(
        unlist(lapply(colnames(feature_data),
                      function(x) !x %in% TAXONOMY_RANKS)))){
        # Find ranks
        ranks <- lapply(colnames(feature_data),
                        .replace_colnames_based_on_prefix, x=feature_data)
        # Replace old ranks with found ranks
        colnames(feature_data) <- ranks
    }
    
    # Remove prefixes if specified and rowData includes info
    if(removeTaxaPrefixes && ncol(feature_data) > 0){
        # Patterns for superkingdom, domain, kingdom, phylum, class, order, family,
        # genus, species
        patterns <- "sk__|([dkpcofgs]+)__"
        feature_data <- lapply(
            feature_data,
            gsub, pattern = patterns, replacement = "")
        feature_data <- as.data.frame(feature_data)
    }
    
    # Adjust row and colnames
    rownames(counts) <- rownames(feature_data) <- biomformat::rownames(obj)
    colnames(counts) <- rownames(sample_data) <- biomformat::colnames(obj)
    
    # Convert into DataFrame
    sample_data <- DataFrame(sample_data)
    feature_data <- DataFrame(feature_data)
    # Convert into list
    assays <- SimpleList(counts = counts)
    
    # Create TreeSE
    tse <- TreeSummarizedExperiment(
        assays = assays,
        colData = sample_data,
        rowData = feature_data)
    return(tse)
}

####################### makeTreeSummarizedExperimentFromBiom #######################
#' @param obj object of type \code{\link[biomformat:read_biom]{biom}}
#' @rdname makeTreeSEFromBiom
#' @export
makeTreeSummarizedExperimentFromBiom <- function(obj, ...){
    makeTreeSEFromBiom(obj, ...)
}

################################ HELP FUNCTIONS ################################
# Find taxonomy rank based on prefixes. If found, return
# corresponding rank. Otherwise, return the original
# rank that is fed to function.
.replace_colnames_based_on_prefix <- function(colname, x){
    # Get column
    col = x[ , colname]
    # List prefixes
    prefixes <- c(
        "^d__",
        "^k__",
        "^p__",
        "^c__",
        "^o__",
        "^f__",
        "^g__",
        "^s__"
    )
    # Find which prefix is found from each column value, if none.
    found_rank <- lapply(
        prefixes, FUN = function(pref){all(grepl(pattern=pref, col))})
    found_rank <- unlist(found_rank)
    # If only one prefix was found (like it should be), get the corresponding
    # rank name.
    if( sum(found_rank) ){
        colname <- TAXONOMY_RANKS[found_rank]
        # Make it capitalized
        colname <- paste0(toupper(substr(colname, 1, 1)),
                          substr(colname, 2, nchar(colname)))
    }
    return(colname)    
}
