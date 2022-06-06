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

TreeSE <- function(counts, rowData = NULL, rowTree = NULL, colData = NULL, colTree = NULL, ...){
    ################################ Input check ###############################
    if( missing(counts) ){
        stop("'counts' should be provided.", call. = FALSE)
    }
    ############################## Input check end #############################
    # Check counts
    counts <- .check_counts(counts)
    # Check rowData
    rowData <- .check_rowdata(rowData)
    # Check colData
    colData <- .check_coldata(colData)
    
    if( !is.null(rowData) ){
        # Map rows
        counts_rowData <- .map_rows(counts, rowData)
        counts <- counts_rowData$counts
        rowData <- counts_rowData$rowData
    } else{
        rowData <- S4Vectors::make_zero_col_DFrame(nrow(counts))
        rownames(rowData) <- rownames(counts)
    }
    
    if( !is.null(colData) ){
        # Map cols
        counts_colData <- .map_cols(counts, colData)
        counts <- counts_colData$counts
        colData <- counts_colData$colData
    } else{
        colData <- S4Vectors::make_zero_col_DFrame(ncol(counts))
        rownames(colData) <- colnames(counts)
    }
    
    # Create assays list
    assays <- SimpleList(counts = counts)
    # Create a TreeSE
    tse <- TreeSummarizedExperiment(assays = assays,
                                    rowData = rowData,
                                    colData = colData,
                                    ...
                                    )
    #tse <- .addRowTree(tse, rowTree)
}

################################ HELP FUNCTIONS ################################

# This function checks counts
.check_counts <- function(counts){
    if( is.null(rownames(counts)) ){
        stop("'counts' does not have rownames!", call. = FALSE)
    }
    # Convert into matrix
    counts <- as.matrix(counts)
    return(counts)
}

# This function checks rowData
.check_rowdata <- function(rowData){
    if( !is.null(rowData) ){
        if( is.null(rownames(rowData)) ){
            stop("'rowData' does not have rownames!", call. = FALSE)
        }
        # Convert into DataFrame
        rowData <- DataFrame(rowData)
    }
    return(rowData)
}

# This function checks colData
.check_coldata <- function(colData){
    if( !is.null(colData) ){
        if( is.null(rownames(colData)) ){
            stop("'colData' do not have rownames!", call. = FALSE)
        }
        # Convert into DataFrame
        colData <- DataFrame(colData)
    }
    return(colData)
}

.map_rows <- function(counts, rowData){
    if( !is.null(rowData) ){
        if( any(suppressWarnings( rownames(counts) != rownames(rowData))) ){
            # Is there names that do not match
            mapping_counts <- match(rownames(counts), rownames(rowData))
            missing_names_counts <- NULL
            if( any(is.na(mapping_counts))){
                missing_names_counts <- rownames(counts)[is.na(mapping_counts)]
            }
            
            
            if( !is.null(missing_names_counts) ){
                warning("'counts' included ", length(missing_names_counts), 
                               " features that are missing from rowData.", "
                               \nPlease check these features for errors:\n",
                               paste0(paste0("'", missing_names_counts, "'"), collapse = " and "), call. = FALSE)
                
            }
            
            # Check rowData
            mapping_rowdata <- match(rownames(rowData), rownames(counts))
            missing_names_rowdata <- NULL
            
            
            if( any(is.na(mapping_rowdata))){
                missing_names_rowdata <- rownames(rowData)[is.na(mapping_rowdata)]
            }
            
            if( !is.null(missing_names_rowdata) ){
                warning("'rowData' included ", length(missing_names_rowdata), 
                        " additional festures that are removed. \nPlease check these features for errors:\n",
                        paste0(paste0("'", missing_names_rowdata, "'"), collapse = ", "), call. = FALSE)
            }
            
            if( is.null(missing_names_counts) && is.null(missing_names_rowdata) ){
                warning("'rowData' is ordered based on 'counts'!", call. = FALSE)
            }
            
            # Order rowData based on counts
            rowData <- rowData[ mapping_counts, ]
            rownames(rowData) <- rownames(counts)
            }
    result_list <- list(counts = counts, rowData = rowData)
    return(result_list)
    }
}

.map_cols <- function(counts, colData){
    if( !is.null(colData) ){
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
    }
    result_list <- list(counts = counts, colData = colData)
    return(result_list)
}

.addRowTree <- function(tse, rowTree){
    if( !is.null(rowTree) ){
        if( rowTree ){
            taxa <- rowData(tse)
            taxa_res <- resolveLoop(as.data.frame(taxa))
            taxa_tree <- toTree(data = taxa_res)
            taxa_tree$tip.label <- getTaxonomyLabels(tse)
            rowNodeLab <- getTaxonomyLabels(tse, make_unique = FALSE)
            tse <- changeTree(tse, rowTree = taxa_tree, rowNodeLab = rowNodeLab)
        }
        else if( any(rownames(tse) != rowTree$tip.label) ){
            warning("'rowTree', does not match with rownames of 'counts'. ",
                    "Use 'rowTree' if you want to add generated tree based on counts...", call. = FALSE)
        } else{
            rowTree(tse) <- rowTree
        }
    }
    return(tse)
}

.addColTree <- function(tse, colTree){
    if( !is.null(colTree)  ){
        if( any(rownames(tse) != colTree$tip.label) ){
            warning("'colTree', does not match with colnames of 'counts'. ",
                    "'colTree is not added.", call. = FALSE)
        } else{
            rowTree(tse) <- rowTree
        }
    }
    return(tse)
}
