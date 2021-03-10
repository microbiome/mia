#' Converting a \code{\link[=SummarizedExperiment-class]{SummarizedExperiment}}
#' object into a long data.frame
#'
#' \code{metlAssaay} Converts a
#' \code{\link[=SummarizedExperiment-class]{SummarizedExperiment}} object into a
#' long data.frame which can be used for \code{tidyverse}-tools.
#'
#' @details
#' If the \code{colData} contains a column \dQuote{SampleID} or the
#' \code{rowData} contains a column \dQuote{FeatureID}, they will be renamed to
#' \dQuote{SampleID_col} and \dQuote{FeatureID_row}, if row names or column
#' names are set.
#'
#' @param x A numeric matrix or a
#'   \code{\link[=SummarizedExperiment-class]{SummarizedExperiment}}
#'
#' @param add_col_data \code{NULL}, \code{TRUE} or a \code{character} vector to
#'   select information from the \code{colData} to add to the molten assay data.
#'   If \code{add_col_data = NULL} no data will be added, if
#'   \code{add_col_data = TRUE} all data will be added and if
#'   \code{add_col_data} is a \code{character} vector, it will be used to subset
#'   to given column names in \code{colData}. (default:
#'   \code{add_col_data = NULL})
#'
#' @param add_row_data \code{NULL}, \code{TRUE} or a \code{character} vector to
#'   select information from the \code{rowData} to add to the molten assay data.
#'   If \code{add_row_data = NULL} no data will be added, if
#'   \code{add_row_data = TRUE} all data will be added and if
#'   \code{add_row_data} is a \code{character} vector, it will be used to subset
#'   to given column names in \code{rowData}. (default:
#'   \code{add_row_data = NULL})
#'
#' @param abund_values a \code{character} value to select an
#'   \code{\link[SummarizedExperiment:SummarizedExperiment-class]{assayNames}}
#'
#' @param ... optional arguments currently not used.
#'
#' @return A \code{tibble} with the molten data. The assay values are given in a
#' column named like the selected assay \code{abund_values}. In addition, a
#' column \dQuote{FeatureID} will contain the rownames, if set, and analogously
#' a column \dQuote{SampleID} with the colnames, if set
#'
#' @name meltAssay

#' @author
#' Sudarshan A. Shetty
#'
#' @examples
#' data(GlobalPatterns)
#' molten_se <- meltAssay(GlobalPatterns,
#'                        add_row_data = TRUE,
#'                        add_col_data = TRUE,
#'                        abund_values = "counts")
#' molten_se
NULL

#' @rdname meltAssay
#' @export
setGeneric("meltAssay",
           signature = "x",
           function(x,
                    add_row_data = NULL,
                    add_col_data = NULL,
                    abund_values = "counts", ...)
               standardGeneric("meltAssay")
)

MIA_MELT_ROW_NAME <- "FeatureID"
MIA_MELT_COL_NAME <- "SampleID"

.norm_add_row_data <- function(add_row_data, x){
    if(is.null(add_row_data)){
        return(NULL)
    }

    if(anyNA(add_row_data)){
        stop("'add_row_data' contains NA.", call. = FALSE)
    }

    cn <- colnames(rowData(x))
    if(is.logical(add_row_data) && length(add_row_data) == 1L && add_row_data){
        add_row_data <- cn
    } else if (isFALSE(all(add_row_data %in% cn))) {
        stop("Please provide valid column names with 'add_row_data' matching ",
             "those in 'rowData(x)'", call. = FALSE)
    }

    if(!is.null(rownames(x)) && MIA_MELT_ROW_NAME %in% add_row_data){
        warning("'x' contains a column '",MIA_MELT_ROW_NAME,"' in its ",
                "rowData(), which will ",
                "be renamed to '",MIA_MELT_ROW_NAME,"_row'", call. = FALSE)
    }

    add_row_data
}

.norm_add_col_data <- function(add_col_data, x){
    if(is.null(add_col_data)){
        return(NULL)
    }

    if(anyNA(add_col_data)){
        stop("'add_col_data' contains NA.", call. = FALSE)
    }

    cn <- colnames(colData(x))
    if(is.logical(add_col_data) && length(add_col_data) == 1L && add_col_data){
        add_col_data <- cn
    } else if (isFALSE(all(add_col_data %in% cn))) {
        stop("Please provide valid column names with 'add_col_data' matching ",
             "those in 'colData(x)'", call. = FALSE)
    }

    if(!is.null(colnames(x)) && MIA_MELT_COL_NAME %in% add_col_data){
        warning("'x' contains a column '",MIA_MELT_COL_NAME,"' in its ",
                "colData(), which will ",
                "be renamed to '",MIA_MELT_COL_NAME,"_col'", call. = FALSE)
    }

    add_col_data
}

.col_switch_name <- function(name){
    paste0(name,"_col")
}

.row_switch_name <- function(name){
    paste0(name,"_row")
}

#' @importFrom dplyr mutate select
.format_molten_assay <- function(molten_assay, x){
    if(is.null(rownames(x)) &&
       !is.null(molten_assay[,.row_switch_name(MIA_MELT_ROW_NAME)]) &&
       !anyNA(molten_assay[,.row_switch_name(MIA_MELT_ROW_NAME)]) &&
       !anyDuplicated(rowData(x)[,MIA_MELT_ROW_NAME])){
        molten_assay <- molten_assay %>%
            select(!sym(MIA_MELT_ROW_NAME)) %>%
            dplyr::rename(!!sym(MIA_MELT_ROW_NAME) := !!sym(.row_switch_name(MIA_MELT_ROW_NAME)))
    }
    if(is.null(colnames(x)) &&
       !is.null(molten_assay[,.col_switch_name(MIA_MELT_COL_NAME)]) &&
       !anyNA(molten_assay[,.col_switch_name(MIA_MELT_COL_NAME)]) &&
       !anyDuplicated(colData(x)[,MIA_MELT_COL_NAME])){
        molten_assay %>%
            select(!sym(MIA_MELT_COL_NAME)) %>%
            dplyr::rename(sym(MIA_MELT_COL_NAME) := !!sym(.col_switch_name(MIA_MELT_COL_NAME)))
    }

    molten_assay %>%
        mutate(!!sym(MIA_MELT_ROW_NAME) := factor(!!sym(MIA_MELT_ROW_NAME)),
               !!sym(MIA_MELT_COL_NAME) := factor(!!sym(MIA_MELT_COL_NAME)))
}


#' @rdname meltAssay
#'
#' @export
setMethod("meltAssay", signature = c(x = "SummarizedExperiment"),
    function(x,
             add_row_data = NULL,
             add_col_data = NULL,
             abund_values = "counts", ...) {
        # input check abund_values
        .check_abund_values(abund_values, x)
        # check selected colnames
        add_row_data <- .norm_add_row_data(add_row_data, x)
        add_col_data <- .norm_add_col_data(add_col_data, x)

        molten_assay <- .melt_assay(x, abund_values)

        if(!is.null(add_row_data)){
            molten_assay <-
                .add_row_data_to_molten_assay(molten_assay, x, add_row_data)
        }

        if(!is.null(add_col_data)){
            molten_assay <-
                .add_col_data_to_molten_assay(molten_assay, x, add_col_data)
        }

        .format_molten_assay(molten_assay, x)
    }
)


# Melts assay in SummarizedExperiment
#' @importFrom SummarizedExperiment assay
#' @importFrom tibble rownames_to_column
#' @importFrom tidyr pivot_longer
#' @importFrom rlang sym
.melt_assay <- function(x, abund_values) {
    molten_assay <- assay(x, abund_values) %>%
        data.frame() %>%
        rownames_to_column(MIA_MELT_ROW_NAME) %>%
        # SampleID is unique sample id
        pivot_longer(!sym(MIA_MELT_ROW_NAME),
                     values_to = abund_values,
                     names_to = MIA_MELT_COL_NAME)

    return(molten_assay)
}


# Combines molten assay with rowData i.e. taxonomy table.
#' @importFrom dplyr left_join
.add_row_data_to_molten_assay <- function(molten_assay, x, add_row_data) {
    rd <- .get_row_data_frame(x, add_row_data)

    molten_assay <- molten_assay %>%
        left_join(rd, by = MIA_MELT_ROW_NAME)

    return(molten_assay)
}


# Converts rowData to data.frame to avoid issues with subsetting and rownames.
#' @importFrom SummarizedExperiment rowData
#' @importFrom rlang sym
#' @importFrom tibble rownames_to_column
#' @importFrom dplyr rename
.get_row_data_frame <- function(x, add_row_data){
    rd <- SummarizedExperiment::rowData(x)[,add_row_data] %>%
        data.frame()
    if(MIA_MELT_ROW_NAME %in% add_row_data){
        rd <- rd %>%
            dplyr::rename(!!sym(.row_switch_name(MIA_MELT_ROW_NAME)) := !!sym(MIA_MELT_ROW_NAME))
    }
    rd %>%
        rownames_to_column(MIA_MELT_ROW_NAME)
}


# Combines molten assay and rowData i.e. taxonomy table with
.add_col_data_to_molten_assay <- function(molten_assay, x, add_col_data) {
    cd <- .get_col_data_frame(x, add_col_data)

    molten_assay <- molten_assay %>%
        left_join(cd, by = MIA_MELT_COL_NAME)

    return(molten_assay)
}


# Converts colData to data.frame to avoid issues with
#' @importFrom SummarizedExperiment colData
#' @importFrom rlang sym
#' @importFrom tibble rownames_to_column
#' @importFrom dplyr rename
.get_col_data_frame <- function(x, add_col_data){
    cd <- SummarizedExperiment::colData(x)[,add_col_data] %>%
        data.frame()
    if(MIA_MELT_COL_NAME %in% add_col_data){
        cd <- cd %>%
            dplyr::rename(!!sym(.col_switch_name(MIA_MELT_COL_NAME)) := !!sym(MIA_MELT_COL_NAME))
    }
    cd %>%
        rownames_to_column(MIA_MELT_COL_NAME)
}
