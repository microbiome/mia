#' Converting a
#' \code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#' object into a long data.frame
#'
#' \code{meltSE} Converts a
#' \code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#' object into a long data.frame which can be used for \code{tidyverse}-tools.
#'
#' @details
#' If the \code{colData} contains a column \dQuote{SampleID} or the
#' \code{rowData} contains a column \dQuote{FeatureID}, they will be renamed to
#' \dQuote{SampleID_col} and \dQuote{FeatureID_row}, if row names or column
#' names are set.
#'
#' @inheritParams getDominant
#' @inheritParams getDissimilarity
#'
#' @param add.col \code{Logical scalar}. \code{NULL}, or
#' \code{character vector}. Used to select information from the \code{colData}
#' to add to the molten assay data. If \code{add.col = NULL} no data will
#' be added, if \code{add.col = TRUE} all data will be added and if
#' \code{add.col} is a \code{character} vector, it will be used to subset
#' to given column names in \code{colData}. (Default: \code{NULL})
#'
#' @param add_col_data Deprecated. Use \code{add.col} instead.
#'
#' @param add.row \code{Logical scalar} or \code{Character vector}. To
#' select information from the \code{rowData} to add to the molten assay data.
#' If \code{add.row = NULL} no data will be added, if
#' \code{add.row = TRUE} all data will be added and if
#' \code{add.row} is a \code{character} vector, it will be used to subset
#' to given column names in \code{rowData}. (Default: \code{NULL})
#'
#' @param add_row_data Deprecated. Use \code{add.row} instead.
#'
#' @param add.dimred \code{Logical scalar} or \code{Character vector}. To
#' select information from the \code{reducedDim} to add to the molten assay
#' data. If \code{add.dimred = NULL} no data will be added, if
#' \code{add.dimred = TRUE} all data will be added and if
#' \code{add.dimred} is a \code{character} vector, it will be used to subset
#' to given names in \code{reducedDimNames(x)}. (Default: \code{NULL})
#'
#' @param row.name \code{Character scalar}. To use as the output's name
#' for the feature identifier. (Default: \code{"FeatureID"})
#'
#' @param feature_name Deprecated. Use \code{row.name} instead.
#'
#' @param col.name \code{Character scalar}. To use as the output's name
#' for the sample identifier. (Default: \code{"SampleID"})
#'
#' @param sample_name Deprecated. Use \code{col.name} instead.
#'
#' @param ... optional arguments:
#' \itemize{
#'   \item check.names: \code{Logical scalar}. Passed to data.frame
#'   function's check.name argument. Determines if sample names are checked
#'   that they are syntactically valid variable names and are not duplicated.
#'   If they are not, sample names are modified. (Default: \code{TRUE})
#' }
#'
#' @return A \code{tibble} with the molten data. The assay values are given in a
#' column named like the selected assay \code{assay.type}. In addition, a
#' column \dQuote{FeatureID} will contain the rownames, if set, and analogously
#' a column \dQuote{SampleID} with the colnames, if set.
#'
#' @name meltSE

#'
#' @examples
#' data(GlobalPatterns)
#' molten_tse <- meltSE(
#'     GlobalPatterns,
#'     assay.type = "counts",
#'     add.row = TRUE,
#'     add.col = TRUE
#'     )
#' molten_tse
NULL

#' @rdname meltSE
#'
#' @export
setMethod("meltSE", signature = c(x = "SummarizedExperiment"),
    function(
        x,
        assay.type = assay_name, assay_name = "counts",
        add.row = add_row_data,
        add_row_data = NULL,
        add.col = add_col_data,
        add_col_data = NULL,
        row.name = feature_name,
        feature_name = "FeatureID",
        col.name = sample_name,
        sample_name = "SampleID",
        ...){
        # input check
        .check_assay_present(assay.type, x)
        if( !.is_a_string(row.name) ){
            stop("'row.name' must be a single non-empty character value.",
                call. = FALSE)
        }
        if( !.is_a_string(col.name) ){
            stop("'col.name' must be a single non-empty character value.",
                call. = FALSE)
        }
        # Check row and column names
        x <- .check_dimnames(x, ...)
        # Check selected colnames
        add.row <- .norm_add_row_data(
            add.row, x, row.name, .internal_MARGIN = "row")
        add.col <- .norm_add_row_data(
            add.col, x, col.name, .internal_MARGIN = "col")
        # Input check end
        # Melt the abundance table
        molten_assay <- .melt_assay(x, assay.type, row.name, col.name)
        # Add rowData to melted abundance table
        if( !is.null(add.row) ){
            molten_assay <- .add_row_data_to_molten_assay(
                molten_assay, x, add.row, row.name, .internal_MARGIN = "row")
        }
        # Add colData to melted abundance table
        if( !is.null(add.col) ){
            molten_assay <- .add_row_data_to_molten_assay(
                molten_assay, x, add.col, col.name, .internal_MARGIN = "col")
        }
        # Factorize feature and sample names
        molten_assay <- .format_molten_assay(
            molten_assay, x, row.name, col.name, ...)
        return(molten_assay)
    }
)

#' @rdname meltSE
#'
#' @export
setMethod("meltSE", signature = c(x = "SingleCellExperiment"),
    function(x, add.dimred = NULL, ...){
        add.dimred <- .check_dimred_for_melting(x, add.dimred)
        df <- callNextMethod(x, ...)
        if( !is.null(add.dimred) ){
            df <- .add_dimred_to_melted_data(df, x, add.dimred, ...)
        }
        return(df)
    }
)

################################ HELP FUNCTIONS ################################

# This function converts sample names so that they include only characters that
# are commonly supported.
.check_dimnames <- function(
        x, check.names = check_names, check_names = FALSE, ...) {
    #
    if( !.is_a_bool(check.names) ){
        stop("'check.names' must be TRUE or FALSE.", call. = FALSE)
    }
    #
    # There must be rownames
    if( is.null(rownames(x)) ){
        rownames(x) <- paste0("row", seq_len(nrow(x)))
        warning("'x' did not have rownames. They are added.", call. = FALSE)
    }
    # There must be colnames
    if( is.null(colnames(x)) ){
        colnames(x) <- paste0("col", seq_len(ncol(x)))
        warning("'x' did not have colnames. They are added.", call. = FALSE)
    }
    # Check if rownames are duplicated, and if they are, modify
    if( any(duplicated(rownames(x))) ){
        rownames(x) <- make.unique(rownames(x))
        warning("rownames(x) included duplicates. ",
                "rownames(x) are made unique. ",  call. = FALSE)
    }
    # Check if colnames are duplicated, and if they are, modify
    if( any(duplicated(colnames(x))) ){
        colnames(x) <- make.unique(colnames(x))
        warning("colnames(x) included duplicates. ",
                "colnames(x) are made unique. ",  call. = FALSE)
    }
    if( check.names ){
        colnames(x) <- make.names(colnames(x))
    }
    return(x)
}

# This function checks whether the specified parameters are correct for adding
# rowData or colData.
.norm_add_row_data <- function(add.row, x, row.name, .internal_MARGIN = "row"){
    # Check .internal_MARGIN
    if( !(.is_a_string(.internal_MARGIN) &&
            .internal_MARGIN %in% c("row", "col")) ){
        stop("'.internal_MARGIN' must be 'row' or 'col'.", call. = FALSE)
    }
    MARGIN_FUN <- switch(
        .internal_MARGIN,
        "row" = rowData,
        "col" = colData
        )
    #
    # If user has specified add.row to be FALSE, convert it to NULL to unify
    # input
    if( .is_a_bool(add.row) && !add.row ){
        add.row <- NULL
    }
    # IF user want to add rowData, check parameters
    if( !is.null(add.row) ){
        # The value should contain only string or single boolean; no NA.
        if( anyNA(add.row) ){
            stop("'add.", .internal_MARGIN, "' contains NA.", call. = FALSE)
        }
        # If user specified boolean, all the columns are added
        cn <- colnames(MARGIN_FUN(x))
        if( .is_a_bool(add.row) && add.row ){
            add.row <- cn
        }
        # Ensure that each variable is only once (user might have specified
        # the same value multiple times)
        add.row <- unique(add.row)
        # At this point, the variable includes character vector.
        # Check that the input includes only those values that can be found.
        if( any(!add.row %in% cn ) ){
            stop("All columns specified by 'add.", .internal_MARGIN, "' ",
                "must match with columns in ",.internal_MARGIN, "Data()",
                call. = FALSE)
        }
    }
    return(add.row)
}

# This function melts assay in SummarizedExperiment
#' @importFrom SummarizedExperiment assay
#' @importFrom tibble rownames_to_column
#' @importFrom tidyr pivot_longer
#' @importFrom rlang sym
.melt_assay <- function(x, assay.type, row.name, col.name, ...){
    # Get assay, ensure that it is a matrix. Ensure that it has row and colnames
    # from the TreeSE
    mat <- assay(x, assay.type) |>
        as.matrix()
    rownames(mat) <- rownames(x)
    colnames(mat) <- colnames(x)
    # Convert it to long format
    mat <- mat |>
        # Convert to data.frame, add feature names to column
        data.frame(check.names = FALSE) |>
        rownames_to_column(row.name) |>
        # Melt data to long format. There are now 3 columns: for feature names
        # sample names, and for abundance values.
        pivot_longer(
            !sym(row.name), values_to = assay.type, names_to = col.name)
    return(mat)
}

# Combines molten assay with rowData i.e. taxonomy table or colData i.e.
# sample metadata.
#' @importFrom SummarizedExperiment rowData colData
#' @importFrom rlang sym
#' @importFrom tibble rownames_to_column
#' @importFrom dplyr rename
.add_row_data_to_molten_assay <- function(
        molten_assay, x, add.row, row.name, .internal_MARGIN = "row", ...){
    # Check .internal_MARGIN
    if( !(.is_a_string(.internal_MARGIN) &&
            .internal_MARGIN %in% c("row", "col")) ){
        stop("'.internal_MARGIN' must be 'row' or 'col'.", call. = FALSE)
    }
    MARGIN_FUN <- switch(
        .internal_MARGIN,
        "row" = rowData,
        "col" = colData
    )
    #
    # Get rowData/colData, and only specified columns
    rd <- MARGIN_FUN(x)
    rd <- rd[ , colnames(rd) %in% add.row, drop = FALSE] |>
        data.frame(check.names = FALSE)
    # Now get those variables from rowData that match with the
    # user-specified columns. For instance, it might be that there are
    # columns with duplicated names. They are made unique..
    if( any(duplicated(colnames(rd))) ){
        warning(.internal_MARGIN, "Data() includes duplicated column ",
                "names that are converted to unique.", call. = FALSE)
        colnames(rd) <- make.unique(colnames(rd))
    }
    # If the name of feature column in melted data is the same as one of the
    # columns in rowData, rename the column in rowData.
    if( row.name %in% colnames(rd) ){
        warning("'x' contains a column '", row.name, "' in its ",
                "'", .internal_MARGIN, "Data(), which will be renamed to '",
                row.name, "_", .internal_MARGIN, "'", call. = FALSE)
        rd <- rd |>
            dplyr::rename(
                !!sym(paste0(row.name, "_", .internal_MARGIN))
                := !!sym(row.name))
    }
    # Add feature names to column in rowData
    rd <- rd |>
        rownames_to_column(row.name)
    # Add rowData to melted data
    molten_assay <- molten_assay |>
        dplyr::left_join(rd, by = row.name)
    return(molten_assay)
}

# This function factorize feature and sample names and modifies colnames if
# specified.
#' @importFrom dplyr mutate sym
.format_molten_assay <- function(
        molten_assay, x, row.name, col.name, check.names = check_names,
        check_names = FALSE, ...) {
    #
    if( !.is_a_bool(check.names) ){
        stop("'check.names' must be TRUE or FALSE.", call. = FALSE)
    }
    #
    # Factorize feature names and sample names
    molten_assay <- molten_assay |>
        mutate(
            !!sym(row.name) := factor(!!sym(row.name)),
            !!sym(col.name) := factor(!!sym(col.name)))
    # If user wants to ensure that column names include only supported
    # characters. E.g., spaces are replaces with dot.
    if( check.names ){
        colnames(molten_assay) <- make.names( colnames(molten_assay) )
    }
    return(molten_assay)
}

# This function validates the add.dimred parameter. The output is standardized.
.check_dimred_for_melting <- function(x, dimred){
    # If FALSE, convert to NULL
    if( length(dimred) == 1L && is.logical(dimred) && !dimred ){
        dimred <- NULL
    }
    # The input can be NULL, i.e., dimred is not added
    if( !is.null(dimred) ){
        # Check that there are dimred in the TreeSE
        dimnames <- reducedDimNames(x)
        if( length(dimnames) == 0L ){
            stop("'add.dimred' is specified but reducedDim(x) is empty.",
                call. = FALSE)
        }
        # If the value is TRUE, user want to add all the dimreds
        if( length(dimred) == 1L && is.logical(dimred) && dimred ){
            dimred <- dimnames
        }
        # If user has specified certain values, check that they exist...
        dimnames <- setNames(
            c(dimnames, seq_len(length(dimnames))), c(dimnames, dimnames))
        if( !all(dimred %in% dimnames) ){
            stop("All reduced dimensionalities specified by 'dimred' must be ",
                "present in reducedDim(x).", call. = FALSE)
        }
        # And match with the values. This works also for indices, that is why
        # we do it like this.
        dimred <- names(dimnames)[ match(dimred, dimnames) ] |> unique()
    }
    return(dimred)
}

# This function adds dimred to data.frame. The inptu is already validated to
# work.
#' @importFrom dplyr rename
#' @importFrom tibble rownames_to_column
.add_dimred_to_melted_data <- function(
        df, x, dimred, col.name = sample_name, sample_name = "SampleID",...){
    # Loop through all dimreds and get the results
    dim_data <- lapply(dimred, function(dim) reducedDim(x, dim))
    names(dim_data) <- dimred
    # Get all the column names and check if there are duplicated names. If there
    # are, add suffix from dimred
    column_names <- lapply(dim_data, colnames) |> unlist()
    if( anyDuplicated(column_names) ){
        warning("reducedDim() includes duplicated column ",
                "names that are converted to unique.", call. = FALSE)
        dim_data <- lapply(dimred, function(dim){
            temp <- dim_data[[dim]]
            colnames(temp) <- paste0(colnames(temp), "_", dim)
            return(temp)
        })
    }
    # Convert to data.frame
    dim_data <- do.call(cbind, dim_data) |> as.data.frame()
    # Check if reducedDim includes column named similarly to column including
    # sample names. If there are, modify the column from reducedDim.
    if( col.name %in% colnames(dim_data) ){
        warning("'x' contains a column '", col.name, "' in its ",
                "reducedDim which will be renamed to '",
                col.name, "_dimred'", call. = FALSE)
        dim_data <- dim_data |> rename(
            !!sym(paste0(col.name, "_dimred")) := !!sym(col.name))
    }
    # Add reducedDm data to the table
    dim_data <- dim_data |>
        rownames_to_column(col.name)
    df <- df |>
        dplyr::left_join(dim_data, by = col.name)
    return(df)
}