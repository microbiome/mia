#' Merge TreeSE objects into single TreeSE object.
#' 
#' @param x a \code{\link{SummarizedExperiment}} object or a list of 
#' \code{\link{SummarizedExperiment}} objects.
#' 
#' @param y a \code{\link{SummarizedExperiment}} object 
#' 
#' @param abund_values A single character value for selecting the
#' \code{\link[SummarizedExperiment:SummarizedExperiment-class]{assay}}
#' to be merged.
#' 
#' @param missing_values NA, 0, or a single character values specifying the notation
#' of missing values.
#'
#'
#' @param BPPARAM A
#'   \code{\link[BiocParallel:BiocParallelParam-class]{BiocParallelParam}}
#'   object specifying whether calculation of estimates should be parallelized.
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
#' @references
#'
#' @seealso
#'
#' @name mergeTreeSummarizedExperiment
#' @export
#'
#' @author Leo Lahti and Tuomas Borman. Contact: \url{microbiome.github.io}
#' 
#' @examples
#' data(GlobalPatterns)
#' tse <- GlobalPatterns
#' 
#' 
NULL

################################### Generic ####################################

#' @rdname mergeTreeSummarizedExperiment
#' @export
setGeneric("mergeTreeSummarizedExperiment", signature = c("x"),
        function(x, ... )
            standardGeneric("mergeTreeSummarizedExperiment"))

###################### Function for SimpleList of TreeSEs ######################

#' @rdname mergeTreeSummarizedExperiment
#' @export
setMethod("mergeTreeSummarizedExperiment", signature = c(x = "SimpleList"),
        function(x, abund_values = "counts", missing_values = 0, ..., BPPARAM = SerialParam() ){
            ################## Input check ##################
            # Can the abund_value the found form all the objects
            abund_values_bool <- lapply(x, .assay_cannot_be_found, abund_values = abund_values)
            if( any(abund_values_bool) ){
                stop("'abund_values' must specify an assay from assays. 'abund_values' ",
                     "cannot be found at least in one TreeSE.",
                     call. = FALSE)
            }
            # Is missing_values one of the allowed ones
            missing_values_bool <- (is.numeric(missing_values) && missing_values = 0) ||
                .is_a_string(missing_values) || is.na(missing_values)
            # If not then give error
            if(  missing_values_bool ){
                stop("'missing_values' must be 0, NA, or a single character value.",
                     call. = FALSE)
            }
            ################ Input check end ################
            # Take first element and remove it from the list
            tse <- x[[1]]
            x[[1]] <- NULL
            # Get rowTree to include if rownames match with all objects
            tse <- as(tse, "TreeSummarizedExperiment")
            rownames <- rownames(tse)
            rowtree <- rowTree(tse)
            # Remove all information but rowData, colData, and assay
            rowdata <- rowData(tse)
            colData(tse) <- colData(tse)
            assay <- assay(tse, abund_values)
            assays <- SimpleList(name = assay)
            names(assays) <- abund_values
            tse <- TreeSummarizedExperiment(assays = assays,
                                            rowData = rowdata,
                                            colData = colData
                                            )
            
            for(i in list ){
                tse <- .merge_TreeSE(i, tse_original = tse, abund_values = abund_values, 
                                     missing_values = missing_values)
            }
            # IF rownames match with original, add rowtree
            if( all(rownames(tse) == rownames) ){
                rowTree(tse) <- rowtree
            }
            
            return(tse)
        }
)

########################### Function for two TreeSEs ###########################

#' @rdname mergeTreeSummarizedExperiment
#' @export
setMethod("mergeTreeSummarizedExperiment", signature = c(x = "SummarizedExperiment"),
        function(x, y = NULL, ...){
            ################## Input check ##################
            # Check y
            class <- class(y)
            if( !(class == "SummarizedExperiment" || 
                  class = "TreeSummarizedExperiment") ){
                stop("'y' must be a 'SummarizedExperiment' object.",
                     call. = FALSE)
            } 
            ################ Input check end ################
            # Create a list based on TreeSEs
            list <- SimpleList(x, y)
            # Call the function for list
            mergeTreeSummarizedExperiment(list, ...)
        }
)

########################### Function for list TreeSEs ##########################

#' @rdname mergeTreeSummarizedExperiment
#' @export
setMethod("mergeTreeSummarizedExperiment", signature = c(x = "list"),
          function(x, ...){
              # Convert into a list
              x <- SimpleList(x)
              # Call the function for list
              mergeTreeSummarizedExperiment(x, ...)
          }
)

################################ HELP FUNCTIONS ################################

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

.merge_TreeSE <- function(tse_original, tse, abund_values, missing_values){
    
    # Merge rowData
    rowdata <- .merge_rowdata(tse_original, tse)
    # Merge colData
    coldata <- .merge_coldata(tse_original, tse)
    # Merge assay
    assay <- .merge_assay(tse_original, tse, abund_values, missing_values, rowdata, coldata)
    assays <- SimpleList(name = assay)
    names(assays) <- abund_values
    
    # Create TreeSE from the data
    tse <- TreeSummarizedExperiment(assays = assays,
                                    rowData = rowdata,
                                    colData = coldata)
    return(tse)
    
}

.merge_assay <- function(tse_original, tse, abund_values, missing_values, rd, cd){
    # Take assays
    assay1 <- assay(tse_original, abund_values)
    assay2 <- assay(tse, abund_values)
    
    # MErge two assays into one
    assay <- .full_join_two_tables(assay1, assay2)
    
    # Fill missing values
    assay[ is.na(assay) ] <- missing_values
    # Convert into matrix
    assay <- as.matrix(assay)
    
    # Order the assay based on rowData and colData
    assay <- assay[ match(rownames(rd), rownames(assay)), ]
    assay <- assay[ , match(rownames(cd), colnames(assay)) ]
    
    return(assay)
}

.merge_rowdata <- function(tse_original, tse){
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
    
    # Merge rowDatas
    rd <- .full_join_two_tables(rd1, rd2)
    # Convert into DataFrame
    rd <- DataFrame(rd)
    
    # Get column indices that match with taxonomy ranks
    ranks_ind <- match( TAXONOMY_RANKS, colnames(rd) )
    # Remove NAs
    ranks_ind <- ranks_ind[ !is.na(ranks_ind) ]
    # Get the data in correct order, take only column that have ranks
    rd_rank <- rd[ , ranks_ind]
    # Take other columns
    rd_other <- rd[ , !ranks_ind]
    
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

.merge_coldata <- function(tse_original, tse){
    # Take colDatas
    cd1 <- colData(tse_original)
    cd2 <- colData(tse)
    # Merge coldata
    cd <- .full_join_two_tables(cd1, cd2)
    # Convert into DataFrame
    cd <- DataFrame(cd)
    
    return(cd)
}

#' @importFrom dplyr coalesce
.full_join_two_tables <- function(df1, df2){
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
    colnames(df1)[ matching_variables_ids1 ] <- matching_variables_mod1
    matching_variables_mod2 <- paste0(matching_variables2, "_Y")
    colnames(df2)[ matching_variables_ids2 ] <- matching_variables_mod2
    
    # Add rownames to one of the columns
    df1$rownames_merge_ID <- rownames(df1)
    df2$rownames_merge_ID <- rownames(df2)
    # Merge data frames into one data frame
    df <- merge(df1, df2, by = "rownames_merge_ID", all = TRUE)
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
            df[ , x] <- NULL
            df[ , y] <- NULL
            # Add column that has combined information
            df[ , matching_variables1[i] ] <- x_and_y_combined
        }
    }
    return(df)
}
