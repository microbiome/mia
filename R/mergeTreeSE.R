#' Merge TreeSE objects into single TreeSE object.
#' 
#' @param x a \code{\link{SummarizedExperiment}} object or a list of 
#' \code{\link{SummarizedExperiment}} objects.
#' 
#' @param y a \code{\link{SummarizedExperiment}} object 
#' 
#' @param abund_values A single character value for selecting the
#' \code{\link[SummarizedExperiment:SummarizedExperiment-class]{assay}}
#' to be merged. (By default: \code{abund_values = "counts"})
#' 
#' @param join A single character value for selecting the joining method.
#' Must be 'full', 'inner', 'left', or 'right'. 'left' and 'right' are disabled
#' when more than two objects are being merged.  (By default: \code{join = "full"})
#' 
#' @param missing_values NA, 0, or a single character values specifying the notation
#' of missing values. (By default: \code{missing_values = 0})
#' 
#' @param verbose A single boolean value to choose whether to show messages. 
#' (By default: \code{verbose = TRUE})
#'
#' @param ... optional arguments (not used).
#'
#' @return A single \code{TreeSummarizedExperiment} object.
#'
#' @details
#' This function merges multiple \code{SummarizedExperiment} objects. It combines
#' \code{rowData}, \code{assays}, and \code{colData} so that the output includes
#' each unique row and column ones. If, for example, all rows are not shared with
#' individual objects, there are missing values in \code{assays}. The notation of missing
#' can be specified with the \code{missing_values} argument. 
#' 
#' Compared to \code{cbind} and \code{rbind} \code{mergeTreeSE} 
#' allows more freely merging since \code{cbind} and \code{rbind} expect 
#' that rows and columns are matching, respectively.
#' 
#' You can choose joining methods from 'full', 'inner',
#' 'left', and 'right'. In all the methods, all the samples are included in the 
#' result object. However, with different methods, it is possible to choose which 
#' rows are included.
#' 
#' \itemize{
#'   \item{\code{full} -- all unique features}
#'   \item{\code{inner} -- all shared features}
#'   \item{\code{left} -- all the features of the first object}
#'   \item{\code{right} -- all the features of the second object}
#' }
#'
#' @seealso
#' \itemize{
#'   \item{\code{TreeSummarizedExperiment::cbind}}
#'   \item{\code{TreeSummarizedExperiment::rbind}}
#' }
#'
#' @name mergeTreeSE
#' @export
#'
#' @author Leo Lahti and Tuomas Borman. Contact: \url{microbiome.github.io}
#' 
#' @examples
#' data(GlobalPatterns)
#' data(esophagus)
#' data(enterotype)
#' 
#' # Take only subsets so that it wont take so long
#' tse1 <- GlobalPatterns[1:100, ]
#' tse2 <- esophagus
#' tse3 <- enterotype[1:100, ]
#' 
#' # Merge two TreeSEs
#' tse <- mergeTreeSE(tse1, tse2)
#' 
#' # Merge a list of TreeSEs
#' list <- SimpleList(tse1, tse2, tse3)
#' tse <- mergeTreeSE(list, abund_values = "counts", missing_values = 0)
#' tse
#' 
#' # With 'join', it is possible to specify the merging method. Subsets are used
#' # here just to show the functionality
#' tse <- mergeTreeSE(tse[1:10, 1:10],  tse[5:100, 11:20], join = "left")
#' tse
#' 
NULL

################################### Generic ####################################

#' @rdname mergeTreeSE
#' @export
setGeneric("mergeTreeSE", signature = c("x"),
        function(x, ... )
            standardGeneric("mergeTreeSE"))

###################### Function for SimpleList of TreeSEs ######################

#' @rdname mergeTreeSE
#' @export
setMethod("mergeTreeSE", signature = c(x = "SimpleList"),
        function(x, abund_values = "counts", join = "full", 
                 missing_values = 0, verbose = TRUE, ... ){
            ################## Input check ##################
            # Can the abund_value the found form all the objects
            abund_values_bool <- lapply(x, .assay_cannot_be_found, abund_values = abund_values)
            abund_values_bool <- unlist(abund_values_bool)
            if( any(abund_values_bool) ){
                stop("'abund_values' must specify an assay from assays. 'abund_values' ",
                     "cannot be found at least in one TreeSE.",
                     call. = FALSE)
            }
            # Check join
            if( !(.is_a_string(join) &&
                join %in% c("full", "inner", "left", "right") ) ){
                stop("'join' must be 'full', 'inner', 'left', or 'right'.",
                     call. = FALSE)
            }
            # Check if join is not available
            if( length(x) > 2 && length(join) != 1L && 
                !join %in% c("full", "inner") ){
                stop("Joining method 'left' and 'right' are not available ",
                     "when more than two objects are being merged.",
                     call. = FALSE)
            }
            # Is missing_values one of the allowed ones
            missing_values_bool <- length(missing_values) == 1L &&
                (is.numeric(missing_values) && missing_values == 0) ||
                .is_a_string(missing_values) || is.na(missing_values)
            # If not then give error
            if(  !missing_values_bool ){
                stop("'missing_values' must be 0, NA, or a single character value.",
                     call. = FALSE)
            }
            # Check verbose
            if( !.is_a_bool(verbose) ){
                stop("'verbose' must be TRUE or FALSE.",
                     call. = FALSE)
            }
            ################ Input check end ################
            # Give message if TRUE
            if( verbose ){
                message("Merging with ", join, " join...\n1/", length(x))
            }
            # Take first element and remove it from the list
            tse <- x[[1]]
            x[[1]] <- NULL
            # Get rowTree to include if match with result data
            tse <- as(tse, "TreeSummarizedExperiment")
            row_tree <- rowTree(tse)
            # Remove all information but rowData, colData, metadata and assay
            row_data <- rowData(tse)
            col_data <- colData(tse)
            assay <- assay(tse, abund_values)
            assays <- SimpleList(name = assay)
            names(assays) <- abund_values
            metadata <- metadata(tse)
            
            tse <- TreeSummarizedExperiment(assays = assays,
                                            rowData = row_data,
                                            colData = col_data,
                                            metadata = metadata,
                                            rowTree = row_tree
                                            )
            
            # Lopp through individual TreeSEs and add them to tse
            if( length(x) > 0 ){
                for( i in 1:length(x) ){
                    # Give message if TRUE
                    if( verbose ){
                        message(i+1, "/", length(x)+1)
                    }
                    temp <- x[[i]]
                    tse <- .merge_TreeSE(temp, tse_original = tse, abund_values = abund_values,
                                         join = join,
                                         missing_values = missing_values)
                }
            }

            return(tse)
        }
)

########################### Function for two TreeSEs ###########################

#' @rdname mergeTreeSE
#' @export
setMethod("mergeTreeSE", signature = c(x = "SummarizedExperiment"),
        function(x, y = NULL, ...){
            ################## Input check ##################
            # Check y
            class <- class(y)
            if( !(class == "SummarizedExperiment" || 
                  class == "TreeSummarizedExperiment") ){
                stop("'y' must be a 'SummarizedExperiment' object.",
                     call. = FALSE)
            } 
            ################ Input check end ################
            # Create a list based on TreeSEs
            list <- SimpleList(x, y)
            # Call the function for list
            mergeTreeSE(list, ...)
        }
)

########################### Function for list TreeSEs ##########################

#' @rdname mergeTreeSE
#' @export
setMethod("mergeTreeSE", signature = c(x = "list"),
          function(x, ...){
              # Convert into a list
              x <- SimpleList(x)
              # Call the function for list
              mergeTreeSE(x, ...)
          }
)

################################ HELP FUNCTIONS ################################
############################ .assay_cannot_be_found #############################
# This function checks that the assay can be found from TreeSE. If it cannot be found
# --> TRUE, if it can be found --> FALSE

# Input: the name of the assay and TreSE object
# Output: TRUE or FALSE
.assay_cannot_be_found <- function(abund_values, tse){
    # Check if the abund_values can be found. If yes, then FALSE. If not, then TRUE
    tryCatch(
        {
            .check_assay_present(abund_values, tse)
            return(FALSE)
            
        },
        error = function(cond) {
            return(TRUE)
        }
    )
}

################################ .merge_TreeSE #################################
# This function merges two TreeSE objects into one.

# Input: Two TreeSEs, the name of the assay, joining method, and the value to
# denote missing values that might occur when object do not share same features, e.g.
# Output: A single TreeSE
.merge_TreeSE <- function(tse_original, tse, abund_values, join, missing_values){
    # Merge rowData
    rowdata <- .merge_rowdata(tse_original, tse, join)
    # Merge colData
    coldata <- .merge_coldata(tse_original, tse, join)
    # Merge assay
    assay <- .merge_assay(tse_original, tse, abund_values, join, missing_values, rowdata, coldata)
    assays <- SimpleList(name = assay)
    names(assays) <- abund_values
    # Combine metadata
    metadata <- c( metadata(tse_original), metadata(tse) )
    # Get row trees
    row_tree1 <- rowTree(tse_original)
    row_tree2 <- rowTree(tse)
    
    # Create TreeSE from the data
    tse <- TreeSummarizedExperiment(assays = assays,
                                    rowData = rowdata,
                                    colData = coldata,
                                    metadata = metadata)
    
    # If labels of the 1st tree match with data, add tree1
    if( !is.null( rownames(tse) ) && !is.null( row_tree1 ) &&
        all( rownames(tse) %in% row_tree1$tip.label ) ){
        rowTree(tse) <- row_tree1
    # If labels of the 2nd tree match with data, add tree2
    } else if( !is.null( rownames(tse) ) && !is.null( row_tree2 ) &&
               all( rownames(tse) %in% row_tree2$tip.label ) ){
        rowTree(tse) <- row_tree2
    }
    
    return(tse)
}

################################ .merge_assay ##################################
# This function merges assays.

# Input: Two TreeSEs, the name of the assay, joining method, value to denote
# missing values, merged rowData, and merged colData
# Output: Merged assay
.merge_assay <- function(tse_original, tse, abund_values, join,
                         missing_values, rd, cd){
    # Take assays
    assay1 <- assay(tse_original, abund_values)
    assay2 <- assay(tse, abund_values)
    
    # Merge two assays into one
    assay <- .join_two_tables(assay1, assay2, join)
    
    # Fill missing values
    assay[ is.na(assay) ] <- missing_values
    # Convert into matrix
    assay <- as.matrix(assay)
    
    # Order the assay based on rowData and colData
    assay <- assay[ match(rownames(rd), rownames(assay)), , drop = FALSE ]
    assay <- assay[ , match(rownames(cd), colnames(assay)), drop = FALSE]
    
    return(assay)
}

############################### .merge_rowdata #################################
# This function merges rowDatas,

# Input: Two TreeSEs and joining method
# Output: Merged rowData
.merge_rowdata <- function(tse_original, tse, join){
    # Take rowDatas
    rd1 <- rowData(tse_original)
    rd2 <- rowData(tse)
    
    # Convert column names to lower
    if( length(colnames(rd1)) > 0 ){
        colnames(rd1) <- tolower(colnames(rd1))
    }
    if( length(colnames(rd2)) > 0 ){
        colnames(rd2) <- tolower(colnames(rd2))
    }
    
    # Merge rowdata
    rd <- .join_two_tables(rd1, rd2, join)
    
    # Get column indices that match with taxonomy ranks
    ranks_ind <- match( TAXONOMY_RANKS, colnames(rd) )
    # Remove NAs
    ranks_ind <- ranks_ind[ !is.na(ranks_ind) ]
    # Get the data in correct order, take only column that have ranks
    rd_rank <- rd[ , ranks_ind, drop = FALSE]
    # Take other columns
    rd_other <- rd[ , !ranks_ind, drop = FALSE]
    
    # Get rank names
    rank_names <- colnames(rd_rank)
    # Convert names s that they have capital letters
    new_rank_names <- paste(toupper(substr(rank_names, 1, 1)), 
                            substr(rank_names, 2, nchar(rank_names)), sep = "")
    # Add new names to colnames of rd_rank
    colnames(rd_rank) <- new_rank_names
    
    # Combine columns
    rd <- cbind(rd_rank, rd_other)
    
    
    return(rd)
}

############################### .merge_coldata #################################
# This function merges colDatas,

# Input: Two TreeSEs and joining method
# Output: Merged colData
.merge_coldata <- function(tse_original, tse, join){
    # Take colDatas
    cd1 <- colData(tse_original)
    cd2 <- colData(tse)
    
    # Merge coldata
    cd <- .join_two_tables(cd1, cd2, join = "full")
    # Convert into DataFrame
    cd <- DataFrame(cd)
    
    return(cd)
}

############################## .join_two_tables ################################
# This general function is used to merge rowDatas, colDatas, and assays.

# Input: Two tables and joining method
# Output: One merged table

#' @importFrom dplyr coalesce
.join_two_tables <- function(df1, df2, join){
    # Get parameter based on join
    all.x <- switch(join,
                    full = TRUE,
                    inner = FALSE,
                    left = TRUE,
                    right = FALSE
    )
    all.y <- switch(join,
                    full = TRUE,
                    inner = FALSE,
                    left = FALSE,
                    right = TRUE
    )
    # Ensure that the data is in correct format
    df1 <- as.data.frame(df1)
    df2 <- as.data.frame(df2)
    
    # Get matching variables indices
    matching_variables_ids1 <- match( colnames(df2), colnames(df1) )
    # Get matching variable names
    matching_variables1 <- colnames(df1)[ matching_variables_ids1 ]
    # Remove NAs
    matching_variables1 <- matching_variables1[ !is.na(matching_variables1) ]
    
    # Get matching variables indices
    matching_variables_ids2 <- match( colnames(df1), colnames(df2) )
    # Get matching variable names
    matching_variables2 <- colnames(df2)[ matching_variables_ids2 ]
    # Remove NAs
    matching_variables2 <- matching_variables2[ !is.na(matching_variables2) ]
    
    # Make the matching variables unique
    matching_variables_mod1 <- paste0(matching_variables1, "_X")
    matching_variables_ids1 <- matching_variables_ids1[ !is.na(matching_variables_ids1) ]
    colnames(df1)[ matching_variables_ids1 ] <- matching_variables_mod1
    matching_variables_mod2 <- paste0(matching_variables2, "_Y")
    matching_variables_ids2 <- matching_variables_ids2[ !is.na(matching_variables_ids2) ]
    colnames(df2)[ matching_variables_ids2 ] <- matching_variables_mod2
    
    # Add rownames to one of the columns
    df1$rownames_merge_ID <- rownames(df1)
    df2$rownames_merge_ID <- rownames(df2)
    # Merge data frames into one data frame
    df <- merge(df1, df2, by = "rownames_merge_ID", all.x = all.x, all.y = all.y)
    # Add rownames and remove additional column
    rownames(df) <- df$rownames_merge_ID
    df$rownames_merge_ID <- NULL
    
    # Combine matching variables if found
    if( length(matching_variables1) > 0 ){
        for(i in 1:length(matching_variables1) ){
            # Get columns
            x <- matching_variables_mod1[i]
            y <- matching_variables_mod2[i]
            # Combine information from columns
            x_and_y_combined <- coalesce( df[ , x], df[ , y] )
            # Remove additional columns
            df[ , x ] <- NULL
            df[ , y ] <- NULL
            # Add column that has combined information
            df[ , matching_variables1[i] ] <- x_and_y_combined
        }
    }
    return(df)
}
