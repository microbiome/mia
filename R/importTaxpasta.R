#' Import taxpasta-specific BIOM results to `TreeSummarizedExperiment`
#'
#' @param filename A character vector with one element that denotes the file
#'   path to a BIOM file.
#'
#' @return A [TreeSummarizedExperiment::TreeSummarizedExperiment()] object.
#' @export
#'
#' @examples
#' tse <- importTaxpasta("profiles.biom")
#'
#' @importFrom rhdf5 h5read
#' @importFrom SummarizedExperiment rowData
#' @importFrom SingleCellExperiment altExps
importTaxpasta <- function(filename) {
    # Validate the input.
    if (!.is_non_empty_string(filename)) {
        stop("'filename' must be a single character value.", call. = FALSE)
    }
    if (!file.exists(filename)) {
        stop(filename, " not found.", call. = FALSE)
    }

    # We read our own HDF5 array to later be able to read observation group
    # metadata, which [biomformat::read_biom()] currently doesn't do.
    raw <- h5read(filename, "/", read.attributes = TRUE)
    biom <- .createBiom(raw)

    # Without taxonomic information, we return a simple TreeSE.
    if (is.null(raw$observation$`group-metadata`$ranks)) {
        warning("The BIOM file does not contain taxonomy information; unable to generate a taxonomic tree.")
        return(makeTreeSEFromBiom(biom))
    }

    # With taxonomic information, we add a hierarchy to the TreeSE.
    ranks <- .getRanks(raw)
    tse <- makeTreeSEFromBiom(biom)
    rowData(tse) <- .createRowData(biom, ranks)
    setTaxonomyRanks(ranks)
    tse <- addHierarchyTree(tse)
    altExps(tse) <- splitByRanks(tse, agglomerate.tree = TRUE)

    return(tse)
}

#' Create BIOM object from raw HDF5 array
#'
#' `.createBiom()` is an internal function used by [importTaxpasta()], that more or
#' less replicates [biomformat::read_hdf5_biom()].
#'
#' @param h5array A raw HDF5 array read from a BIOM file.
#'
#' @return A [biomformat::biom()] object.
#'
#' @keywords internal
#' @noRd
.createBiom <- function(h5array) {
    data <- biomformat:::generate_matrix(h5array)
    rows <- biomformat:::generate_metadata(h5array$observation)
    columns <- biomformat:::generate_metadata(h5array$sample)
    shape <- c(length(data), length(data[[1]]))

    id <- attr(h5array, "id")
    vs <- attr(h5array, "format-version")
    format <- sprintf("Biological Observation Matrix %s.%s", vs[1], vs[2])
    format_url <- attr(h5array, "format-url")
    type <- attr(h5array, "type")
    generated_by <- attr(h5array, "generated-by")
    date <- attr(h5array, "creation-date")
    matrix_type <- "dense"
    matrix_element_type <- "int"

    return(
        biomformat:::namedList(
            id, format, format_url, type, generated_by, date, matrix_type,
            matrix_element_type, rows, columns, shape, data
        ) |> biomformat:::biom()
    )
}

#' Get taxonomic ranks from observation group metadata
#'
#' `.getRanks()` is an internal function used by [importTaxpasta()], that reads the
#' string of semi-colon separated ranks from the group observation metadata.
#'
#' @param h5array A raw HDF5 array read from a BIOM file.
#'
#' @return A character vector of taxonomic ranks in order.
#'
#' @keywords internal
#' @noRd
.getRanks <- function(h5array) {
    ranks <- strsplit(h5array$observation$`group-metadata`$ranks, ";", fixed = TRUE)[[1]]
    return(ranks)
}

#' Recreate observation metadata
#'
#' `.createRowData()` is an internal function used by [importTaxpasta()], that
#' returns a copy of the BIOM object's observation metadata, where the headers
#' of the taxonomy columns have been replaced with the correct taxonomic rank
#' names.
#'
#' @param biom A BIOM object.
#' @param ranks A character vector of taxonomic ranks in order.
#'
#' @return A [data.frame()] with observation metadata.
#'
#' @importFrom biomformat observation_metadata
#'
#' @keywords internal
#' @noRd
.createRowData <- function(biom, ranks) {
    meta <- observation_metadata(biom)
    column.names <- colnames(meta)
    indeces <- startsWith(column.names, "taxonomy")

    if (sum(indeces) != length(ranks)) {
        stop("The number of generic taxonomy* columns differs from the number of ranks.")
    }

    colnames(meta) <- replace(column.names, indeces, ranks)
    return(meta)
}
