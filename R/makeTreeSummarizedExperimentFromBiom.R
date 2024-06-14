#' Loading a biom file
#'
#' For convenience a few functions are available to convert data from a
#' \sQuote{biom} file or object into a
#' \code{\link[TreeSummarizedExperiment:TreeSummarizedExperiment-class]{TreeSummarizedExperiment}}
#'
#' @param file biom file location
#' 
#' @param prefix.rm \code{TRUE} or \code{FALSE}: Should
#' taxonomic prefixes be removed? The prefixes is removed only from detected
#' taxa columns meaning that \code{rank.from.prefix} should be enabled in the most cases.
#' (default \code{prefix.rm = FALSE})
#' 
#' @param removeTaxaPrefixes Deprecated. Use \code{prefix.rm} instead.
#' 
#' @param rank.from.prefix \code{TRUE} or \code{FALSE}: If file does not have
#' taxonomic ranks on feature table, should they be scraped from prefixes?
#' (default \code{rank.from.prefix = FALSE})
#' 
#' @param rankFromPrefix Deprecated.Use \code{rank.from.prefix} instead.
#' 
#' @param artifact.rm \code{TRUE} or \code{FALSE}: If file have
#' some taxonomic character naming artifacts, should they be removed.
#' (default \code{artifact.rm = FALSE})
#' 
#' @param ... additional arguments 
#'   \itemize{
#'        \item \code{patter}: \code{character} value specifying artifacts
#'        to be removed. If \code{patterns = "auto"}, special characters
#'        are removed. (default: \code{pattern = "auto"})
#'    }
#' 
#' @return An object of class
#'   \code{\link[TreeSummarizedExperiment:TreeSummarizedExperiment-class]{TreeSummarizedExperiment}}
#'
#' @name makeTreeSEFromBiom
#' @seealso
#' \code{\link[=makeTreeSEFromPhyloseq]{makeTreeSEFromPhyloseq}}
#' \code{\link[=makeTreeSEFromDADA2]{makeTreeSEFromDADA2}}
#' \code{\link[=importQIIME2]{importQIIME2}}
#' \code{\link[=importMothur]{importMothur}}
#'
#' @examples
#' # Load biom file
#' library(biomformat)
#' biom_file <- system.file("extdata", "rich_dense_otu_table.biom",
#'                          package = "biomformat")
#' 
#' # Make TreeSE from biom file
#' tse <- importBIOM(biom_file)
#' 
#' # Make TreeSE from biom object
#' biom_object <- biomformat::read_biom(biom_file)
#' tse <- makeTreeSEFromBiom(biom_object)
#' 
#' # Get taxonomyRanks from prefixes and remove prefixes
#' tse <- importBIOM(biom_file,
#'                     rank.from.prefix = TRUE,
#'                     prefix.rm = TRUE)
#' 
#' # Load another biom file
#' biom_file <- system.file("extdata/testdata", "Aggregated_humanization2.biom",
#'                          package = "mia")
#' 
#' # Clean artifacts from taxonomic data
#' tse <- importBIOM(biom_file,
#'                     artifact.rm = TRUE)
NULL

#' @rdname makeTreeSEFromBiom
#'
#' @export
importBIOM <- function(file, ...) {
    .require_package("biomformat")
    biom <- biomformat::read_biom(file)
    makeTreeSEFromBiom(biom, ...)
}

#' @rdname makeTreeSEFromBiom
#'
#' @param x object of type \code{\link[biomformat:read_biom]{biom}}
#'
#' @export
#' @importFrom S4Vectors make_zero_col_DFrame DataFrame
#' @importFrom dplyr %>% bind_rows
makeTreeSEFromBiom <- function(
        x, prefix.rm = removeTaxaPrefixes, 
        removeTaxaPrefixes = FALSE, rank.from.prefix = rankFromPrefix, 
        rankFromPrefix = FALSE,
        artifact.rm = FALSE, ...){
    # input check
    .require_package("biomformat")
    if(!is(x,"biom")){
        stop("'x' must be a 'biom' object", call. = FALSE)
    }
    if( !.is_a_bool(prefix.rm) ){
        stop("'prefix.rm' must be TRUE or FALSE.", call. = FALSE)
    }
    if( !.is_a_bool(rank.from.prefix) ){
        stop("'rank.from.prefix' must be TRUE or FALSE.", call. = FALSE)
    }
    if( !.is_a_bool(artifact.rm) ){
        stop("'artifact.rm' must be TRUE or FALSE.", call. = FALSE)
    }
    #
    counts <- as(biomformat::biom_data(x), "matrix")
    sample_data <- biomformat::sample_metadata(x)
    feature_data <- biomformat::observation_metadata(x)
    
    # colData is initialized with empty tables with rownames if it is NULL
    if( is.null(sample_data) ){
        sample_data <- S4Vectors::make_zero_col_DFrame(ncol(counts))
        rownames(sample_data) <- colnames(counts)
    # Otherwise convert it into correct format if it is a list
    } else if( is(sample_data, "list") ){
        # Merge list of data.frames into one
        sample_data <- bind_rows(sample_data)
        sample_data < as.data.frame(sample_data)
    }
    # rowData is initialized with empty tables with rownames if it is NULL
    if( is.null(feature_data) ){
        feature_data <- S4Vectors::make_zero_col_DFrame(nrow(counts))
        rownames(feature_data) <- rownames(counts)
    # Otherwise convert it into correct format if it is a list
    } else if( is(feature_data, "list") ){
        # Feature data is a list of taxa info. Dfs are merged together 
        # differently than sample metadata since the column names are only 
        # "Taxonomy". If there is only one taxonomy level, the column name does 
        # not get a suffix.
        # --> bind rows based on the index of column.
        
        # Get the maximum length of list
        max_length <- max( lengths(feature_data) )
        # Get the column names from the taxa info that has all the levels that 
        # occurs in the data
        colnames <- names( head(
            feature_data[ lengths(feature_data) == max_length ], 1)[[1]])
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
        # Transposing feature_data and make it df object
        feature_data <- as.data.frame(feature_data)
        # Add column that includes all the data
        feature_data[["taxonomy_unparsed"]] <- apply(feature_data, 1, paste0, collapse = ";")
        # Add correct colnames
        colnames(feature_data) <- c(colnames, "taxonomy_unparsed")
    }
    # If there is only one column in the feature data, it is the most probable
    # that the taxonomy is not parsed. Try to parse it.
    if( ncol(feature_data) == 1 ){
        colnames(feature_data) <- "taxonomy_unparsed"
        tax_tab <- .parse_taxonomy(feature_data, column_name = colnames(feature_data))
        feature_data <- cbind(tax_tab, feature_data)
        feature_data <- as.data.frame(feature_data)
    }
    
    # Clean feature_data from possible character artifacts if specified
    if( artifact.rm ){
        feature_data <- .detect_taxa_artifacts_and_clean(feature_data, ...)
    }
    
    # Replace taxonomy ranks with ranks found based on prefixes
    if( rank.from.prefix && all(
        unlist(lapply(colnames(feature_data),
                        function(x) !x %in% TAXONOMY_RANKS)))){
        # Find ranks
        ranks <- lapply(colnames(feature_data),
                        .replace_colnames_based_on_prefix, x=feature_data)
        # Replace old ranks with found ranks
        colnames(feature_data) <- unlist(ranks)
    }
    
    # Remove prefixes if specified and rowData includes info
    if(prefix.rm && ncol(feature_data) > 0){
        feature_data <- .remove_prefixes_from_taxa(feature_data, ...)
    }
    
    # Adjust row and colnames
    rownames(counts) <- rownames(feature_data) <- biomformat::rownames(x)
    colnames(counts) <- rownames(sample_data) <- biomformat::colnames(x)
    
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

####################### makeTreeSummarizedExperimentFromBiom ###################
#' @param x object of type \code{\link[biomformat:read_biom]{biom}}
#' @rdname makeTreeSEFromBiom
#' @export
makeTreeSummarizedExperimentFromBiom <- function(x, ...){
    makeTreeSEFromBiom(x, ...)
}

################################ HELP FUNCTIONS ################################
# This function removes prefixes from taxonomy names
.remove_prefixes_from_taxa <- function(
        feature_tab, prefixes = "sk__|([dkpcofgs]+)__",
        only.taxa.col = TRUE, ...){
    if( !.is_a_bool(only.taxa.col) ){
        stop("'only.taxa.col' must be TRUE or FALSE.", call. = FALSE)
    }
    #
    # Subset by taking only taxonomy info if user want to remove the pattern 
    # only from those. (Might be too restricting, e.g., if taxonomy columns are 
    # not detected in previous steps. That is way the default is FALSE)
    if( only.taxa.col ){
        ind <- tolower(colnames(feature_tab)) %in% TAXONOMY_RANKS
        temp <- feature_tab[, ind, drop = FALSE]
    } else{
        ind <- rep(TRUE, ncol(feature_tab))
        temp <- feature_tab
    }
    
    # If there are columns left for removing the pattern
    if( ncol(temp) > 0 ){
        # Remove patterns
        temp <- lapply(
            temp, gsub, pattern = prefixes, replacement = "")
        temp <- as.data.frame(temp)
        # If cell had only prefix, it is now empty string. Convert to NA
        temp[ temp == "" ] <- NA
        # Combine table
        feature_tab[, ind] <- temp
    }
    return(feature_tab)
}

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
    found_rank <- lapply(prefixes, FUN = function(pref){
        all(grepl(pattern = pref, col) | is.na(col)) && !all(is.na(col))
        })
    found_rank <- unlist(found_rank)
    # If only one prefix was found (like it should be), get the corresponding
    # rank name.
    if( sum(found_rank) == 1 ){
        colname <- TAXONOMY_RANKS[found_rank]
        # Make it capitalized
        colname <- paste0(toupper(substr(colname, 1, 1)),
                            substr(colname, 2, nchar(colname)))
    }
    return(colname)    
}

# Detect and clean non wanted characters from Taxonomy data if needed.
.detect_taxa_artifacts_and_clean <- function(x, pattern = "auto", ...) {
    #
    if( !.is_non_empty_character(pattern) ){
        stop("'pattern' must be a single character value.", call. = FALSE)
    }
    #
    row_names <- rownames(x)
    # Remove artifacts
    if( pattern == "auto" ){
        .require_package("stringr")
        # Remove all but these characters
        pattern <- "[[:alnum:]]|-|_|\\[|\\]|,|;\\||[[:space:]]"
        x <- lapply(x, function(col){
            # Take all specified characters as a matrix where each column is a 
            # character
            temp <- stringr::str_extract_all(col, pattern = pattern, simplify = TRUE)
            # Collapse matrix to strings
            temp <- apply(temp, 1, paste, collapse = "")
            # Now NAs are converted into characters. Convert them back
            temp[ temp == "NA" ] <- NA
            # Convert also empty strings to NA
            temp[ temp == "" ] <- NA
            return(temp)
        })
    } else{
        # Remove pattern specified by user
        x <- lapply(x, gsub, pattern = pattern, replacement = "")
    }
    x <- as.data.frame(x)
    # Add rownames because they are dropped while removing artifacts
    rownames(x) <- row_names
    return(x)
}
