#' Converting a \code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#' object into a long data.frame
#'
#' \code{meltSE} Converts a
#' \code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}} object into a
#' long data.frame which can be used for \code{tidyverse}-tools.
#'
#' @details
#' If the \code{colData} contains a column \dQuote{SampleID} or the
#' \code{rowData} contains a column \dQuote{FeatureID}, they will be renamed to
#' \dQuote{SampleID_col} and \dQuote{FeatureID_row}, if row names or column
#' names are set.
#'
#' @param x A numeric matrix or a
#'   \code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#'   
#' @param assay.type a \code{character} value to select an
#'   \code{\link[SummarizedExperiment:SummarizedExperiment-class]{assayNames}}
#'   
#' @param assay_name a single \code{character} value for specifying which
#'   assay to use for calculation.
#'   (Please use \code{assay.type} instead. At some point \code{assay_name}
#'   will be disabled.)
#'
#' @param add.col \code{NULL}, \code{TRUE} or a \code{character} vector to
#'   select information from the \code{colData} to add to the molten assay data.
#'   If \code{add.col = NULL} no data will be added, if
#'   \code{add.col = TRUE} all data will be added and if
#'   \code{add.col} is a \code{character} vector, it will be used to subset
#'   to given column names in \code{colData}. (default:
#'   \code{add.col = NULL})
#' 
#' @param add_col_data Deprecated. Use \code{add.col} instead.
#'
#' @param add.row \code{NULL}, \code{TRUE} or a \code{character} vector to
#'   select information from the \code{rowData} to add to the molten assay data.
#'   If \code{add.row = NULL} no data will be added, if
#'   \code{add.row = TRUE} all data will be added and if
#'   \code{add.row} is a \code{character} vector, it will be used to subset
#'   to given column names in \code{rowData}. (default:
#'   \code{add.row = NULL})
#' 
#' @param add_row_data Deprecated. Use \code{add.row} instead.
#'
#' @param row.name a \code{character} scalar to use as the output's name
#'   for the feature identifier. (default: \code{row.name = "FeatureID"})
#' 
#' @param feature_name Deprecated. Use \code{row.name} instead.
#'
#' @param col.name a \code{character} scalar to use as the output's name
#'   for the sample identifier. (default: \code{col.name = "SampleID"})
#' 
#' @param sample_name Deprecated. Use \code{col.name} instead.
#'
#' @param ... optional arguments:
#' \itemize{
#'   \item check_names: A boolean value passed to data.frame function's check.name
#'   argument. Determines if sample names are checked that they are syntactically 
#'   valid variable names and are not duplicated. If they are not, sample names 
#'   are modified. (default: \code{check_names = TRUE})
#' }
#'
#' @return A \code{tibble} with the molten data. The assay values are given in a
#' column named like the selected assay \code{assay.type}. In addition, a
#' column \dQuote{FeatureID} will contain the rownames, if set, and analogously
#' a column \dQuote{SampleID} with the colnames, if set
#'
#' @name meltSE

#' @author
#' Sudarshan A. Shetty
#'
#' @examples
#' data(GlobalPatterns)
#' molten_tse <- meltSE(GlobalPatterns,
#'                         assay.type = "counts",
#'                         add.row = TRUE,
#'                         add.col = TRUE
#'                         )
#' molten_tse
NULL

#' @rdname meltSE
#' @export
setGeneric("meltSE",
           signature = "x",
           function(x,
                    assay.type = assay_name, assay_name = "counts",
                    add.row = add_row_data,
                    add_row_data = NULL,
                    add.col = add_col_data,
                    add_col_data = NULL,
                    row.name = feature_name,
                    feature_name = "FeatureID",
                    col.name = sample_name,
                    sample_name = "SampleID",
                    ...)
               standardGeneric("meltSE")
)

.norm_add_row_data <- function(add.row, x, row.name){
    if(is.null(add.row)){
        return(NULL)
    }
    if(anyNA(add.row)){
        stop("'add.row' contains NA.", call. = FALSE)
    }
    cn <- colnames(rowData(x))
    if(is.logical(add.row) && length(add.row) == 1L && add.row){
        add.row <- cn
    } else if (isFALSE(all(add.row %in% cn))) {
        stop("Please provide valid column names with 'add.row' matching ",
             "those in 'rowData(x)'", call. = FALSE)
    }
    if(!is.null(rownames(x)) && row.name %in% add.row){
        warning("'x' contains a column '",row.name,"' in its ",
                "rowData(), which will ",
                "be renamed to '",row.name,"_row'", call. = FALSE)
    }
    add.row
}

.norm_add_col_data <- function(add.col, x, col.name){
    if(is.null(add.col)){
        return(NULL)
    }
    if(anyNA(add.col)){
        stop("'add.col' contains NA.", call. = FALSE)
    }
    cn <- colnames(colData(x))
    if(is.logical(add.col) && length(add.col) == 1L && add.col){
        add.col <- cn
    } else if (isFALSE(all(add.col %in% cn))) {
        stop("Please provide valid column names with 'add.col' matching ",
             "those in 'colData(x)'", call. = FALSE)
    }
    if(!is.null(colnames(x)) && col.name %in% add.col){
        warning("'x' contains a column '",col.name,"' in its ",
                "colData(), which will ",
                "be renamed to '",col.name,"_col'", call. = FALSE)
    }
    add.col
}

.col_switch_name <- function(name){
    paste0(name,"_col")
}

.row_switch_name <- function(name){
    paste0(name,"_row")
}

#' @importFrom dplyr mutate select
.format_molten_assay <- function(molten_assay, x,
                                 row.name,
                                 col.name){
    if(is.null(rownames(x)) &&
       .row_switch_name(row.name) %in% colnames(molten_assay) &&
       !anyNA(molten_assay[,.row_switch_name(row.name)]) &&
       !anyDuplicated(rowData(x)[,row.name])){
        molten_assay <- molten_assay %>%
            select(!sym(row.name)) %>%
            dplyr::rename(!!sym(row.name) := !!sym(.row_switch_name(row.name)))
    }
    if(is.null(colnames(x)) &&
       .col_switch_name(col.name) %in% colnames(molten_assay) &&
       !anyNA(molten_assay[,.col_switch_name(col.name)]) &&
       !anyDuplicated(colData(x)[,col.name])){
        molten_assay %>%
            select(!sym(col.name)) %>%
            dplyr::rename(!!sym(col.name) := !!sym(.col_switch_name(col.name)))
    }
    molten_assay %>%
        mutate(!!sym(row.name) := factor(!!sym(row.name)),
               !!sym(col.name) := factor(!!sym(col.name)))
}


#' @rdname meltSE
#'
#' @export
setMethod("meltSE", signature = c(x = "SummarizedExperiment"),
    function(x,
            assay.type = assay_name, assay_name = "counts", 
            add.row = NULL,
            add.col = NULL,
            row.name = feature_name,
            feature_name = "FeatureID",
            col.name = sample_name,
            sample_name = "SampleID",
            ...) {
        # input check
        .check_assay_present(assay.type, x)
        if(!.is_a_string(row.name)){
            stop("'row.name' must be a single non-empty character value.",
                 call. = FALSE)
        }
        if(!.is_a_string(col.name)){
            stop("'col.name' must be a single non-empty character value.",
                 call. = FALSE)
        }
        # check if rownames are duplicated, and if they are, modify
        if( any(duplicated(rownames(x))) ){
            rownames(x) <- make.unique(rownames(x))
            warning("rownames(x) included duplicates.",
                    " rownames(x) are made unique. ",
                    call. = FALSE)
        }
        # check selected colnames
        add.row <- .norm_add_row_data(add.row, x, row.name)
        add.col <- .norm_add_col_data(add.col, x, col.name)
        molten_assay <- .melt_assay(x, assay.type, row.name, col.name, ...)
        if(!is.null(add.row)){
            molten_assay <-
                .add_row_data_to_molten_assay(molten_assay, x, add.row,
                                              row.name)
        }
        if(!is.null(add.col)){
            molten_assay <-
                .add_col_data_to_molten_assay(molten_assay, x, add.col,
                                              col.name, ...)
        }
        .format_molten_assay(molten_assay, x, row.name, col.name)
    }
)

# Melts assay in SummarizedExperiment
#' @importFrom SummarizedExperiment assay
#' @importFrom tibble rownames_to_column
#' @importFrom tidyr pivot_longer
#' @importFrom rlang sym
.melt_assay <- function(x, assay.type, row.name, col.name, check_names = FALSE,...) {
    mat <- assay(x, assay.type) %>%
        as.matrix() 
    rownames(mat) <- rownames(x)
    colnames(mat) <- colnames(x)
    mat %>%
        data.frame(check.names = check_names) %>%
        rownames_to_column(row.name) %>%
        # SampleID is unique sample id
        pivot_longer(!sym(row.name),
                     values_to = assay.type,
                     names_to = col.name)
}

# Combines molten assay with rowData i.e. taxonomy table.
#' @importFrom SummarizedExperiment rowData
#' @importFrom rlang sym
#' @importFrom tibble rownames_to_column
#' @importFrom dplyr rename
.add_row_data_to_molten_assay <- function(molten_assay, x, add.row,
                                          row.name) {
    rd <- SummarizedExperiment::rowData(x)[,add.row,drop=FALSE] %>%
        data.frame()
    if(row.name %in% add.row){
        rd <- rd %>%
            dplyr::rename(!!sym(.row_switch_name(row.name)) := !!sym(row.name))
    }
    rd <- rd %>%
        rownames_to_column(row.name)
    molten_assay %>%
        dplyr::left_join(rd, by = row.name)
}

# Combines molten assay and rowData i.e. taxonomy table with
#' @importFrom SummarizedExperiment colData
#' @importFrom rlang sym
#' @importFrom tibble rownames_to_column
#' @importFrom dplyr rename
.add_col_data_to_molten_assay <- function(molten_assay, x, add.col,
                                          col.name, check_names = FALSE,...) {
    cd <- SummarizedExperiment::colData(x)[,add.col,drop=FALSE] %>%
        data.frame()
    # This makes sure that sample names match
    if(check_names == TRUE){
        rownames(cd) <- make.names(rownames(cd))
    }
    if(col.name %in% add.col){
        cd <- cd %>%
            dplyr::rename(!!sym(.col_switch_name(col.name)) := !!sym(col.name))
    }
    cd <- cd %>%
        rownames_to_column(col.name)
    molten_assay %>%
        dplyr::left_join(cd, by = col.name)
}
