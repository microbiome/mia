#' Import taxpasta-specific BIOM results to
#' \code{\link[TreeSummarizedExperiment:TreeSummarizedExperiment-class]{TreeSummarizedExperiment}}
#' 
#' @name importTaxpasta
#' 
#' @details
#' \code{importTaxpasta} imports data that is returned from TAXonomic Profile
#' Aggregation and STAndardisation (taxpasta) pipeline. See more information on
#' taxpasta from
#' \href{https://taxpasta.readthedocs.io/en/latest/}{taxpasta documentation}.
#'
#' @param file A single \code{character} value specifying the file path to a
#' BIOM file.
#'
#' @return A \code{\link[TreeSummarizedExperiment:TreeSummarizedExperiment-class]{TreeSummarizedExperiment}} object.
#'
#' @examples
#' 
#' # File path to BIOM file
#' file_path <- system.file("extdata/testdata/complete.biom", package = "mia")
#' # Import BIOM as TreeSE
#' tse <- importTaxpasta(file_path)
#'
#' @seealso
#' \code{\link[=importBIOM]{importBIOM}}
#' \code{\link[=convertFromBIOM]{convertFromBIOM}}
#'
NULL

#'
#' @rdname importTaxpasta
#' @export
#'
#' @importFrom SummarizedExperiment rowData
importTaxpasta <- function(file) {
    # Check dependencies.
    .require_package("rhdf5")
    .require_package("biomformat")
    
    # Validate the input.
    if(!.is_non_empty_string(file) ){
        stop("'filename' must be a single character value.", call. = FALSE)
    }
    if( !file.exists(file) ){
        stop("'", file, "' not found.", call. = FALSE)
    }
    
    # We read our own HDF5 array to later be able to read observation group
    # metadata, which [biomformat::read_biom()] currently doesn't do.
    raw <- rhdf5::h5read(file, "/", read.attributes = TRUE)
    biom <- .create_biom(raw)
    # Convert BIOM to TreeSE
    tse <- convertFromBIOM(biom)
    
    # IF we have taxonomic information, we add a hierarchy to the TreeSE.
    if( !is.null(raw$observation$`group-metadata`$ranks) ){
        # Set ranks of mia based on found ranks
        ranks <- .get_ranks(raw)
        setTaxonomyRanks(ranks)
        # Create rowData and rowTree
        rowData(tse) <- .create_row_data(biom, ranks)
        tse <- addHierarchyTree(tse)
        # Agglomerate to all existing ranks
        tse <- agglomerateByRanks(tse, agglomerate.tree = TRUE)
    } else{
        # Without taxonomic information, we return a simple TreeSE.
        warning(
            "The BIOM file does not contain taxonomy information;",
            "unable to generate a taxonomic tree.", call. = FALSE)
    }
    return(tse)
}

#' Create BIOM object from raw HDF5 array
#'
#' \code{.create_biom} is an internal function used by
#' \code{\link[=importTaxpasta]{importTaxpasta}}, that more or less replicates
#' \code{\link[biomformat:read_hdf5_biom]{read_hdf5_biom}}.
#'
#' @param h5array A raw HDF5 array read from a BIOM file.
#'
#' @return A \code{\link[biomformat:biom]{biom}} object.
#'
#' @keywords internal
#' @noRd
.create_biom <- function(h5array) {
    # Get abundance data, row data and column data along with shape
    data <- biomformat:::generate_matrix(h5array)
    rows <- biomformat:::generate_metadata(h5array$observation)
    columns <- biomformat:::generate_metadata(h5array$sample)
    shape <- c(length(data), length(data[[1]]))
    # Get metadata on the data
    id <- attr(h5array, "id")
    vs <- attr(h5array, "format-version")
    format <- sprintf("Biological Observation Matrix %s.%s", vs[1], vs[2])
    format_url <- attr(h5array, "format-url")
    type <- attr(h5array, "type")
    generated_by <- attr(h5array, "generated-by")
    date <- attr(h5array, "creation-date")
    matrix_type <- "dense"
    matrix_element_type <- "int"
    # Create BIOM object
    result <- biomformat:::biom(biomformat:::namedList(
        id, format, format_url, type, generated_by, date, matrix_type,
        matrix_element_type, rows, columns, shape, data
    ))
    return(result)
}

#' Get taxonomic ranks from observation group metadata
#'
#' \code{.get_ranks} is an internal function used by
#' \code{\link[=importTaxpasta]{importTaxpasta}}, that reads the string of
#' semi-colon separated ranks from the group observation metadata.
#'
#' @param h5array A raw HDF5 array read from a BIOM file.
#'
#' @return A character vector of taxonomic ranks in order.
#'
#' @keywords internal
#' @noRd
.get_ranks <- function(h5array) {
    ranks <- strsplit(h5array$observation$`group-metadata`$ranks,
        ";",
        fixed = TRUE
    )[[1]]
    return(ranks)
}

#' Recreate observation metadata
#'
#' \code{.create_row_data} is an internal function used by
#' \code{\link[=importTaxpasta]{importTaxpasta}}, that returns a copy of the
#' BIOM object's observation metadata, where the headers of the taxonomy
#' columns have been replaced with the correct taxonomic rank names.
#'
#' @param biom A BIOM object.
#' @param ranks A character vector of taxonomic ranks in order.
#'
#' @return A \code{\link[base:data.frame]{data.frame}} with observation
#' metadata.
#'
#' @keywords internal
#' @noRd
.create_row_data <- function(biom, ranks) {
    # Get taxonomy table
    meta <- biomformat::observation_metadata(biom)
    # Get column names of taxonomy table
    column.names <- colnames(meta)
    # Taxonomy ranks should start with "taxonomy", take all of them
    indeces <- startsWith(column.names, "taxonomy")
    # Check that they match with ranks got from the data
    if( sum(indeces) != length(ranks) ){
        stop(
            "The number of generic taxonomy* columns differs",
            "from the number of ranks.", call. = FALSE)
    }
    # And replace column names with the set of ranks
    colnames(meta) <- replace(column.names, indeces, ranks)
    return(meta)
}
