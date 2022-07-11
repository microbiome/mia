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

TreeSE <- function(counts, rowData = NULL, colData = NULL,  
                   rowTree = NULL, rowNodeLab = NULL, 
                   colTree = NULL, colNodeLab = NULL, ...){
    ################################ Input check ###############################
    if( missing(counts) ){
        stop("'counts' must be provided.", call. = FALSE)
    }
    ############################## Input check end #############################
    # Initialize a list of arguments
    args <- list()
    
    # Check counts
    counts <- .check_counts(counts)
    
    # If rowData is provided
    if( !is.null(rowData) ){
        # Check rowData
        rowData <- .check_rowdata_and_coldata(rowData)
        # Map rows
        .map_rows(counts, rowData)
    # If the rowData is not provided
    } else{
        rowData <- S4Vectors::make_zero_col_DFrame(nrow(counts))
        rownames(rowData) <- rownames(counts)
    }
    temp_list <- list(rowData = rowData)
    args <- append(args, temp_list)
    
    if( !is.null(rowTree) && !.is_a_bool(rowTree) ){
        .check_rowTree_and_colTree(rowTree, rowNodeLab)
        .add_row_tree(args[["rowData"]], rowTree, rowNodeLab)
    }
    temp_list <- list(rowTree = rowTree, rowNodeLab = rowNodeLab)
    args <- append(args, temp_list)
    
    if( !is.null(colData) ){
        # Check colData
        colData <- .check_rowdata_and_coldata(colData)
        # Map cols
        .map_cols(counts, colData, colTree, ...)
    } else{
        colData <- S4Vectors::make_zero_col_DFrame(ncol(counts))
        rownames(colData) <- colnames(counts)
    }
    temp_list <- list(colData = colData)
    args <- append(args, temp_list)
    
    if( !is.null(colTree) ){
        .check_rowTree_and_colTree(colTree, colNodeLab)
        .add_col_tree(args[["colData"]], colTree, colNodeLab)
    }
    temp_list <- list(colTree = colTree, colNodeLab = colNodeLab)
    args <- append(args, temp_list)
    
    # Create assays list
    assays <- SimpleList(counts = counts)
    temp_list <- list(assays = assays)
    args <- append(args, temp_list)
    
    # Create a TreeSE)
    tse <- do.call(TreeSummarizedExperiment, args)
    
    if( .is_a_bool(rowTree) && rowTree ){
        tse <- addTaxonomyTree(tse)
    }
    return(tse)
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

.check_rowTree_and_colTree <- function(tree, nodeLab){
    if( !(class(tree) == "tree" || class(tree) == "phylo") ){
        stop("'", deparse(substitute(tree)), 
             "' is not in supported format. Please check it.", call. = FALSE)
    }
    if( !(is.vector(nodeLab) || is.null(nodeLab)) ){
        stop("'", deparse(substitute(nodeLab)), "' must be a vector or NULL.", call. = FALSE)
    }
}

# Map rows between assay, rowData
.map_rows <- function(counts, rowData){
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
            message("'counts' included ", length(missing_names_counts), 
                           " features that are missing from rowData.", "
                           \nPlease check these features for errors:\n",
                           paste0(paste0("'", missing_names_counts, "'"), collapse = " and "))
            
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
            message("'rowData' included ", length(missing_names_rowdata), 
                    " additional features. \nPlease check these features for errors:\n",
                    paste0(paste0("'", missing_names_rowdata, "'"), collapse = ", "))
        }
        
        if( is.null(missing_names_counts) && is.null(missing_names_rowdata) ){
            message("Feature names between 'counts' and 'rowData' are in different order.")
        }
        # Give error
        stop("'rownames' do no not match!", call. = FALSE)
    }
    return(TRUE)
}



.add_row_tree <- function(rowData, rowTree, rowNodeLab){
    # Check if rowNodeLab match
    if( !is.null(rowNodeLab) ){
        if( ! (length(rowNodeLab) == nrow(rowData) && 
            all(rowNodeLab %in% rowTree$tip.label) )){
            message("'rowNodeLabs' do not match with the number of features ",
                    "or labels cannot be found from the 'rowTree'")
        } else{
            return(TRUE)
        }
    # Check node labels if they match
    } else if( !( all(rownames(rowData) %in% rowTree$tip.label) ) ){
        message("'rowTree' do not include all the features of 'rowData'.")
    } else{
        return(TRUE)
    }
    stop("'rowTree' does not match with features.", call. = FALSE)
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
            message("'counts' included ", length(missing_names_counts), 
                           " samples that are missing from colData.", 
                    "\nPlease check these samples for errors:\n",
                    paste0(paste0("'", missing_names_counts, "'"), collapse = " and "))
            
        }
        
        # Check colData
        mapping_coldata <- match(rownames(colData), colnames(counts))
        missing_names_coldata <- NULL
        
        if( any(is.na(mapping_coldata))){
            missing_names_coldata <- rownames(colData)[is.na(mapping_coldata)]
        }
        
        if( !is.null(missing_names_coldata) ){
            message("'colData' included ", length(missing_names_coldata), 
                           " additional samples. \nPlease check these samples for errors:\n",
                           paste0(paste0("'", missing_names_coldata, "'"), collapse = ", "))
            
        }
        if( is.null(missing_names_counts) && is.null(missing_names_colData) ){
            message("sample names between 'counts' and 'colData' are in different order.")
        }
        # Give error
        stop("'colnames' do no not match!", call. = FALSE)
    }
    return(TRUE)
}

.add_col_tree <- function(colData, colTree, colNodeLab){
    # Check if colNodeLab match
    if( !is.null(colNodeLab) ){
        if( length(colNodeLab) == rownames(colData) && 
            all(colNodeLab %in% colTree$tip.label )){
            message("'colNodeLabs' do not match with the number of samples ",
                    "or labels cannot be found from the 'colTree'.")
        }
    }
    # Check node labels if they match
    if( !( all(rownames(colData) %in% colTree$tip.label) ) ){
        message("'colTree' does not include all the samples of 'colData'.")
    } else{
        return(TRUE)
    }
    stop("'colTree' does not match with samples.", call. = FALSE)
}

