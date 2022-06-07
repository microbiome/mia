#' Import data to \code{TreeSummarizedExperiment}
#
#' 
#' @param ... additional arguments:
#' 
#'
#' @details
#' 
#' @return  A
#' \code{\link[TreeSummarizedExperiment:TreeSummarizedExperiment-class]{TreeSummarizedExperiment}}
#' object
#'
#' @name TreeSE
#' @seealso
#' \code{\link[=makeTreeSummarizedExperimentFromPhyloseq]{makeTreeSummarizedExperimentFromPhyloseq}}
#' \code{\link[=makeSummarizedExperimentFromBiom]{makeSummarizedExperimentFromBiom}}
#' \code{\link[=makeTreeSummarizedExperimentFromDADA2]{makeTreeSummarizedExperimentFromDADA2}}
#' \code{\link[=loadFromQIIME2]{loadFromQIIME2}}
#' \code{\link[=loadFromMothur]{loadFromMothur}}
#'
#' @export
#' @author Leo Lahti and Tuomas Borman. Contact: \url{microbiome.github.io}
#' 
#' @references
#' 
#' @examples
#' 
#' 
NULL

TreeSE <- function(counts, rowData = NULL, colData = NULL,  rowTree = NULL, colTree = NULL, ...){
    ################################ Input check ###############################
    if( missing(counts) ){
        stop("'counts' must be provided.", call. = FALSE)
    }
    ############################## Input check end #############################
    # Initialize a list of arguments
    args <- list()
    
    # Check counts
    counts <- .check_counts(counts)
    
    # If rowSata is provided
    if( !is.null(rowData) ){
        # Check rowData
        rowData <- .check_rowdata_and_coldata(rowData)
        # Map rows
        temp_list <- .map_rows(counts, rowData, rowTree, ...)
        counts <- temp_list$counts
        temp_list$counts <- NULL
        args <- append(args, temp_list)
    # If the rowData is not provided
    } else{
        rowData <- S4Vectors::make_zero_col_DFrame(nrow(counts))
        rownames(rowData) <- rownames(counts)
        temp_list <- list(rowData = rowData)
        args <- append(args, temp_list)
    }
    
    if( !is.null(colData) ){
        # Check colData
        colData <- .check_rowdata_and_coldata(colData)
        # Map cols
        temp_list <- .map_cols(counts, colData, colTree, ...)
        counts <- temp_list$counts
        temp_list$counts <- NULL
        args <- append(args, temp_list)
    } else{
        colData <- S4Vectors::make_zero_col_DFrame(ncol(counts))
        rownames(colData) <- colnames(counts)
        temp_list <- list(colData = colData)
        args <- append(args, temp_list)
    }
    
    # Create assays list
    assays <- SimpleList(counts = counts)
    temp_list <- list(assays = assays)
    args <- append(args, temp_list)
    #args <- append(args, ...)
    
    # Create a TreeSE
    # tse <- TreeSummarizedExperiment(assays = assays,
    #                                 rowData = rowData,
    #                                 colData = colData,
    #                                 ...
    #                                 )
    tse <- do.call(TreeSummarizedExperiment, args)
    tse
    #tse <- .addRowTree(tse, rowTree)
}

################################ HELP FUNCTIONS ################################

# This function checks counts
.check_counts <- function(counts){
    # Check that it is in supported format
    if( !.is_matrix_like(counts) ){
        stop("'counts' is not matrix-like object!",
             "Please check it.", call. = FALSE)
    }
    # Ensure that rownames and colnames are stored
    rownames <- rownames(counts)
    colnames <- colnames(counts)
    # Convert into matrix
    counts <- as.matrix(counts)
    # Assign rownames and colnames back
    rownames(counts) <- rownames
    colnames(counts) <- colnames
    # If rownames are not provided
    if( is.null(rownames(counts)) ){
        stop("'counts' does not have rownames!",
             "Please add them.", call. = FALSE)
    }
    # If colnames are not provided
    if( is.null(colnames(counts)) ){
        stop("'counts' does not have colnames! ",
             "Please add them.", call. = FALSE)
    }
    return(counts)
}

# This function checks rowData and colData
.check_rowdata_and_coldata <- function(table){
    # Check that it is in supported format
    if( !.is_matrix_like(table) ){
        stop("'", deparse(substitute(table)), "' is not matrix-like object!",
             "Please check it.", call. = FALSE)
    }
    # Store rownames and colnames
    rownames <- rownames(table)
    colnames <- colnames(table)
    # Convert into DataFrame
    table <- DataFrame(table)
    # Assign rownames and colnames
    rownames(table) <- rownames
    colnames(table) <- colnames
    # If rownames are not provided
    if( is.null(rownames(table)) ){
        stop("'", deparse(substitute(table)), "' does not have rownames!",
             "Please add them.", call. = FALSE)
    }
    return(table)
}

# Checks if the object is in supported format
.is_matrix_like <- function(table){
    result <- is.matrix(table) || is.data.frame(table) || class(table) == "DFrame"
    return(result)
}

# Map rows between assay, rowData, and rowTree
.map_rows <- function(counts, rowData, rowTree, rowNodeLab = NULL, ...){
    # Check that the rownames match
    if( any(suppressWarnings( rownames(counts) != rownames(rowData))) ){
        
        # Is there names that do not match
        mapping_counts <- match(rownames(counts), rownames(rowData))
        # Store names that are missing. If none are, it is NULL
        missing_names_counts <- NULL
        if( any(is.na(mapping_counts))){
            missing_names_counts <- rownames(counts)[is.na(mapping_counts)]
        }
        # Give warning if there are features that are in counts but that are missing
        # from rowData
        if( !is.null(missing_names_counts) ){
            warning("'counts' included ", length(missing_names_counts), 
                           " features that are missing from rowData.", "
                           \nPlease check these features for errors:\n",
                           paste0(paste0("'", missing_names_counts, "'"), collapse = " and "), call. = FALSE)
            
        }
        
        # Find features that are in rowData but that are not in counts
        mapping_rowdata <- match(rownames(rowData), rownames(counts))
        # Find names, if were found.
        missing_names_rowdata <- NULL
        if( any(is.na(mapping_rowdata))){
            missing_names_rowdata <- rownames(rowData)[is.na(mapping_rowdata)]
        }
        # Give warning if there were additional feature in rowData
        if( !is.null(missing_names_rowdata) ){
            warning("'rowData' included ", length(missing_names_rowdata), 
                    " additional festures that are removed. \nPlease check these features for errors:\n",
                    paste0(paste0("'", missing_names_rowdata, "'"), collapse = ", "), call. = FALSE)
        }
        
        # If there were not additional features
        if( is.null(missing_names_counts) && is.null(missing_names_rowdata) ){
            warning("'rowData' is ordered based on 'counts'!", call. = FALSE)
        }
        # Order rowData based on counts
        rowData <- rowData[ mapping_counts, ]
        rownames(rowData) <- rownames(counts)
    }
    # If rowTree is not NULL, check that
    if( !is.null(rowTree) ){
        if( !(class(rowTree) == "tree" || class(rowTree) == "phylo") ){
            stop("'rowTree' is not in supported format. Please check it.", call. = FALSE)
        }
        # Check if rowNodeLab match
        if( !is.null(rowNodeLab) ){
            if( is.vector(rowNodeLab) && length(rowNodeLab) == rownames(rowData) && 
                length(rowNodeLab) == length(rowTree$tip.label )){
                warning("'rowNodeLabs' do not match with the number of features ",
                        "in 'counts' or it is not a vector. ", 
                        "'rowTree' is not added.", call. = FALSE)
                rowTree <- NULL
                rowNodeLab <- NULL
            }
        # Check node labels if they match
        } else if( !all(rownames(rowData) %in% rowTree$tip.label) ){
            warning("Not all features in 'rowTree' are included in 'counts'. ", 
                    "'rowTree' is not added.", call. = FALSE)
            # Do not add rowTree
            rowTree <- NULL
            rowNodeLab <- NULL
        } else{
            # Ensure that the order is correct
            rowNodeLab <- rowTree$tip.label[match(rownames(rowData) == rowTree$tip.label) ]
        }
    }
    result_list <- list(counts = counts, rowData = rowData, 
                        rowTree = rowTree, rowNodeLab = rowNodeLab)
    return(result_list)
}



.map_cols <- function(counts, colData, colTree, colNodeLab = NULL, ...){
    if( any(suppressWarnings( colnames(counts) != rownames(colData))) ){
        
        # Is there names that do not match
        mapping_counts <- match(colnames(counts), rownames(colData))
        missing_names_counts <- NULL
        
        if( any(is.na(mapping_counts))){
            missing_names_counts <- colnames(counts)[is.na(mapping_counts)]
        }
        
        
        if( !is.null(missing_names_counts) ){
            warning("'counts' included ", length(missing_names_counts), 
                           " samples that are missing from colData.", 
                    "\nPlease check these samples for errors:\n",
                    paste0(paste0("'", missing_names_counts, "'"), collapse = " and "), call. = FALSE)
            
        }
        
        # Check colData
        mapping_coldata <- match(rownames(colData), colnames(counts))
        missing_names_coldata <- NULL
        
        if( any(is.na(mapping_coldata))){
            missing_names_coldata <- rownames(colData)[is.na(mapping_coldata)]
        }
        
        if( !is.null(missing_names_coldata) ){
            warning("'colData' included ", length(missing_names_coldata), 
                           " additional samples that are removed. \nPlease check these samples for errors:\n",
                           paste0(paste0("'", missing_names_coldata, "'"), collapse = ", "), call. = FALSE)
            
        }
        
        if( is.null(missing_names_counts) && is.null(missing_names_coldata) ){
            warning("'colData' is ordered based on 'counts'!", call. = FALSE)
        }
        
        # Order colData based on counts
        colData <- colData[ mapping_counts, ]
        rownames(colData) <- colnames(counts)
    }
    # If colTree is not NULL, check that
    if( !is.null(colTree) ){
        if( !(class(colTree) == "tree" || class(colTree) == "phylo") ){
            stop("'colTree' is not in supported format. Please check it.", call. = FALSE)
        }
        # Check if colNodeLab match
        if( !is.null(colNodeLab) ){
            if( is.vector(colNodeLab) && length(colNodeLab) == rownames(colData) && 
                length(colNodeLab) == length(colTree$tip.label )){
                warning("'colNodeLabs' do not match with the number of samples ",
                        "in 'counts' or it is not a vector. ", 
                        "'colTree' is not added.", call. = FALSE)
                colTree <- NULL
                colNodeLab <- NULL
            }
        }
        # Check node labels if they match
        if( !all(rownames(colData) %in% colTree$tip.label) ){
            warning("Not all samples in 'colTree' are included in 'counts'. ", 
                    "'colTree' is not added.", call. = FALSE)
            # Do not add colTree
            colTree <- NULL
            colNodeLab <- NULL
        } else{
            # Ensure that the order is correct
            colNodeLab <- colTree$tip.label[ match(rownames(colData) == colTree$tip.label) ]
        }
    }
    result_list <- list(counts = counts, colData = colData, 
                        colTree = colTree, colNodeLab = colNodeLab)
    return(result_list)
}
