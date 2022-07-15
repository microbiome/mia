#' Merge SE objects into single SE object.
#' 
#' @param x a \code{\link{SummarizedExperiment}} object or a list of 
#' \code{\link{SummarizedExperiment}} objects.
#' 
#' @param y a \code{\link{SummarizedExperiment}} object when \code{x} is a
#' \code{\link{SummarizedExperiment}} object. Disabled when \code{x} is a list.
#' 
#' @param assay_name A single character value for selecting the
#' \code{\link[SummarizedExperiment:SummarizedExperiment-class]{assay}}
#' to be merged. (By default: \code{assay_name = "counts"})
#' 
#' @param join A single character value for selecting the joining method.
#' Must be 'full', 'inner', 'left', or 'right'. 'left' and 'right' are disabled
#' when more than two objects are being merged.  (By default: \code{join = "full"})
#' 
#' @param missing_values NA, 0, or a single character values specifying the notation
#' of missing values. (By default: \code{missing_values = NA})
#' 
#' @param verbose A single boolean value to choose whether to show messages. 
#' (By default: \code{verbose = TRUE})
#'
#' @param ... optional arguments (not used).
#'
#' @return A single \code{SummarizedExperiment} object.
#'
#' @details
#' This function merges multiple \code{SummarizedExperiment} objects. It combines
#' \code{rowData}, \code{assays}, and \code{colData} so that the output includes
#' each unique row and column ones. If, for example, all rows are not shared with
#' individual objects, there are missing values in \code{assays}. The notation of missing
#' can be specified with the \code{missing_values} argument. If input consists of
#' \code{TreeSummarizedExperiment} objects, also \code{rowTree}, \code{colTree}, and
#' \code{referenceSeq} are preserved if possible.
#' 
#' Compared to \code{cbind} and \code{rbind} \code{mergeSEs} 
#' allows more freely merging since \code{cbind} and \code{rbind} expect 
#' that rows and columns are matching, respectively.
#' 
#' You can choose joining methods from \code{'full'}, \code{'inner'},
#'  \code{'left'}, and  \code{'right'}. In all the methods, all the samples are 
#'  included in the result object. However, with different methods, it is possible 
#'  to choose which rows are included.
#' 
#' \itemize{
#'   \item{\code{full} -- all unique features}
#'   \item{\code{inner} -- all shared features}
#'   \item{\code{left} -- all the features of the first object}
#'   \item{\code{right} -- all the features of the second object}
#' }
#' 
#' You can also doe e.g., a full join by using a function \code{full_join} which is 
#' an alias for \code{mergeSEs}. Also other joining methods have dplyr-like aliases.
#' 
#' The output depends on the input. If the input contains \code{SummarizedExperiment}
#' object, then the output will be \code{SummarizedExperiment}. When all the input
#' objects belong to \code{TreeSummarizedExperiment}, the output will be 
#' \code{TreeSummarizedExperiment}.
#'
#' @seealso
#' \itemize{
#'   \item{\code{TreeSummarizedExperiment::cbind}}
#'   \item{\code{TreeSummarizedExperiment::rbind}}
#'   \item{\code{\link[dplyr:full_join]{full_join}}}
#'   \item{\code{\link[dplyr:inner_join]{inner_join}}}
#'   \item{\code{\link[dplyr:left_join]{left_join}}}
#'   \item{\code{\link[dplyr:right_join]{right_join}}}
#' }
#'
#' @name mergeSEs
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
#' tse <- mergeSEs(tse1, tse2)
#' 
#' # Merge a list of TreeSEs
#' list <- SimpleList(tse1, tse2, tse3)
#' tse <- mergeSEs(list, assay_name = "counts", missing_values = 0)
#' tse
#' 
#' # With 'join', it is possible to specify the merging method. Subsets are used
#' # here just to show the functionality
#' tse_temp <- mergeSEs(tse[1:10, 1:10], tse[5:100, 11:20], join = "left")
#' tse_temp
#' 
#' # You can also do a left_join by using alias "left_join"
#' tse_temp <- left_join(tse[1:10, 1:10], tse[5:100, 11:20])
#' 
NULL

################################### Generic ####################################

#' @rdname mergeSEs
#' @export
setGeneric("mergeSEs", signature = c("x"),
        function(x, ... )
            standardGeneric("mergeSEs"))

###################### Function for SimpleList of TreeSEs ######################

#' @rdname mergeSEs
#' @export
setMethod("mergeSEs", signature = c(x = "SimpleList"),
        function(x, assay_name = "counts", join = "full", 
                 missing_values = NA, verbose = TRUE, ... ){
            ################## Input check ##################
            # Check the objects 
            class <- .check_objects_and_give_class(x)
            # Can the assay_name the found form all the objects
            assay_name_bool <- .assays_cannot_be_found(assay_name = assay_name, x)
            if( any(assay_name_bool) ){
                stop("'assay_name' must specify an assay from assays. 'assay_name' ",
                     "cannot be found at least in one SE object.",
                     call. = FALSE)
            }
            # Check join
            if( !(.is_a_string(join) &&
                join %in% c("full", "inner", "left", "right") ) ){
                stop("'join' must be 'full', 'inner', 'left', or 'right'.",
                     call. = FALSE)
            }
            # Check if join is not available
            if( length(x) > 2 &&
                join %in% c("left", "right") ){
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
            # Merge objects
            tse <- .merge_SE(x, class, join, assay_name, missing_values, verbose)
            return(tse)
        }
)

########################### Function for two TreeSEs ###########################

#' @rdname mergeSEs
#' @export
setMethod("mergeSEs", signature = c(x = "SummarizedExperiment"),
        function(x, y = NULL, ...){
            ################## Input check ##################
            # Check y
            if( !(is(y, "SummarizedExperiment")) ){
                stop("'y' must be a 'SummarizedExperiment' object.",
                     call. = FALSE)
            } 
            ################ Input check end ################
            # Create a list based on TreeSEs
            list <- SimpleList(x, y)
            # Call the function for list
            mergeSEs(list, ...)
        }
)

########################### Function for list TreeSEs ##########################

#' @rdname mergeSEs
#' @export
setMethod("mergeSEs", signature = c(x = "list"),
          function(x, ...){
              # Convert into a list
              x <- SimpleList(x)
              # Call the function for list
              mergeSEs(x, ...)
          }
)

################################# full_join ####################################

#' @rdname mergeSEs
#' @export
setGeneric("full_join", signature = c("x"),
    function(x, ...)
        standardGeneric("full_join"))

#' @rdname mergeSEs
#' @export
setMethod("full_join", signature = c(x = "ANY"),
    function(x, ...){
        mergeSEs(x, join = "full", ...)
    }
)

################################# inner_join ###################################

#' @rdname mergeSEs
#' @export
setGeneric("inner_join", signature = c("x"),
    function(x, ...)
        standardGeneric("inner_join"))

#' @rdname mergeSEs
#' @export
setMethod("inner_join", signature = c(x = "ANY"),
    function(x, ...){
        mergeSEs(x, join = "inner", ...)
    }
)

################################# left_join ####################################

#' @rdname mergeSEs
#' @export
setGeneric("left_join", signature = c("x"),
    function(x, ...)
        standardGeneric("left_join"))

#' @rdname mergeSEs
#' @export
setMethod("left_join", signature = c(x = "ANY"),
    function(x, ...){
        mergeSEs(x, join = "left", ...)
    }
)

################################# right_join ###################################

#' @rdname mergeSEs
#' @export
setGeneric("right_join", signature = c("x"),
    function(x, ...)
        standardGeneric("right_join"))

#' @rdname mergeSEs
#' @export
setMethod("right_join", signature = c(x = "ANY"),
    function(x, ...){
        mergeSEs(x, join = "right", ...)
    }
)

################################ HELP FUNCTIONS ################################

################################## .merge_SE ###################################
# This function merges SE objects into one SE

# Input: A list of SEs
# Output: SE
#' @importFrom utils txtProgressBar
.merge_SE <- function(x, class, join, assay_name, missing_values, verbose){
    # Take first element and remove it from the list
    tse <- x[[1]]
    x[[1]] <- NULL
    
    # Get the function based on class
    FUN <- switch(class,
                  TreeSummarizedExperiment = .get_TreeSummarizedExperiment_data,
                  SingleCellExperiment = .get_SingleCellExperiment_data,
                  SummarizedExperiment = .get_SummarizedExperiment_data,
                  )
    
    # Get the data in a list
    args <- do.call(FUN, args = list(tse = tse, assay_name = assay_name))
    
    # Get the function based on class
    FUN_constructor <- switch(class,
                              TreeSummarizedExperiment = TreeSummarizedExperiment,
                              SingleCellExperiment = SingleCellExperiment,
                              SummarizedExperiment = SummarizedExperiment
    )
    tse <- do.call(FUN_constructor, args = args)
    
    # Get the function based on class
    FUN <- switch(class,
                  TreeSummarizedExperiment = .merge_TreeSummarizedExperiments,
                  SingleCellExperiment = .merge_SingleCellExperiments,
                  SummarizedExperiment = .merge_SummarizedExperiments,
    )
    
    # Lopp through individual TreeSEs and add them to tse
    if( length(x) > 0 ){
        # Initialize progressbar if TRUE
        if( verbose ){
            pb <- txtProgressBar(min = 0, max = length(x), initial = 0,  style = 3)
        }
        for( i in 1:length(x) ){
            # Give information on progress if TRUE
            if( verbose ){
                setTxtProgressBar(pb, i)
            }
            # Get the ith object
            temp <- x[[i]]
            # Merge data
            args <- do.call(FUN, args = list(
                tse_original = tse,
                tse = temp,
                join = join,
                assay_name = assay_name,
                missing_values = missing_values
                ))
            # Create an object
            tse <- do.call(FUN_constructor, args = args)
        }
        # Close progress bar if TRUE
        if( verbose ){
            close(pb)
        }
    }
    return(tse)
}

###################### .get_TreeSummarizedExperiment_data ######################
# This function gets the desired data from one TreeSE and creates a list of 
# arguments containing the data

# Input; TreeSE
# Output: A list of arguments
.get_TreeSummarizedExperiment_data <- function(tse, assay_name){
    # Get rowTree and colTree
    row_tree <- rowTree(tse)
    col_tree <- colTree(tse)
    # Get a list of arguments of SCE object
    args <- .get_SingleCellExperiment_data(tse, assay_name)
    # Add TreeSE-specific slots
    args$rowTree <- row_tree
    args$colTree <- col_tree
    return(args)
}

######################## .get_SingleCellExperiment_data ########################
# This function gets the desired data from one SCE object and creates a list of 
# arguments containing the data

# Input; SCE
# Output: A list of arguments
.get_SingleCellExperiment_data <- function(tse, assay_name){
    # reducedDim is additional slot for SCE compared to SE. 
    # However, merging reducedDims leads to non-meaningful data
    # Get the arguments of SE object
    args <- .get_SummarizedExperiment_data(tse, assay_name)
    return(args)
}

######################## .get_SummarizedExperiment_data ########################
# This function gets the desired data from one SE object and creates a list of 
# arguments containing the data

# Input; SE
# Output: A list of arguments
.get_SummarizedExperiment_data <- function(tse, assay_name){
    # Remove all information but rowData, colData, metadata and assay
    row_data <- rowData(tse)
    col_data <- colData(tse)
    assay <- assay(tse, assay_name)
    assays <- SimpleList(name = assay)
    names(assays) <- assay_name
    metadata <- metadata(tse)
    # Create a list of arguments
    args <- list(assays = assays,
                rowData = row_data,
                colData = col_data,
                metadata = metadata
    )
    return(args)
    
}

######################## .check_objects_and_give_class #########################
# This function checks that the object are in correct format

# Input: a list of objects
# Output: A shared class of objects
.check_objects_and_give_class <- function(x){
    # Allowed classes
    allowed_classes <- c("TreeSummarizedExperiment", "SingleCellExperiment", "SummarizedExperiment")
    
    # Get the class based on hierarchy TreeSE --> SCE --> SE
    if( all( unlist( lapply(x, is, class2 = allowed_classes[[1]]) ) ) ){
        class <- allowed_classes[1]
    } else if( all( unlist( lapply(x, is, class2 = allowed_classes[[2]]) ) ) ){
        class <- allowed_classes[2]
    } else if( all( unlist( lapply(x, is, class2 = allowed_classes[[3]]) ) ) ){
        class <- allowed_classes[3]
    # If there is an object that does not belong to these classes give an error
    } else{
        stop("Input includes an object that is not 'SummarizedExperiment'.",
             call. = FALSE)
    }
    # If there are multiple classes, give a warning
    if( length(unique( unlist(lapply(x, function(y){ class(y)})) )) > 1 ){
        warning("The Input consist of multiple classes. ",
                "The output is '", class, "'.",
                call. = FALSE)
    }
    return(class)
}

########################### .assays_cannot_be_found #############################
# This function checks that the assay can be found from TreeSE objects of a list.

# Input: the name of the assay and a list of TreeSE objects
# Output: A list of boolean values
.assays_cannot_be_found <- function(assay_name, x){
    # Check if the assay_name can be found. If yes, then FALSE. If not, then TRUE
    list <- lapply(x, .assay_cannot_be_found, assay_name = assay_name)
    # Unlist the list
    result <- unlist(list)
    return(result)
}

############################ .assay_cannot_be_found #############################
# This function checks that the assay can be found from TreeSE. If it cannot be found
# --> TRUE, if it can be found --> FALSE

# Input: the name of the assay and TreSE object
# Output: TRUE or FALSE
.assay_cannot_be_found <- function(assay_name, tse){
    # Check if the assay_name can be found. If yes, then FALSE. If not, then TRUE
    tryCatch(
        {
            .check_assay_present(assay_name, tse)
            return(FALSE)
            
        },
        error = function(cond) {
            return(TRUE)
        }
    )
}

###################### .merge_TreeSummarizedExperiments ########################
# This function merges the data of two TreeSE objects into one set of arguments that
# can be feed to create a single object.

# Input: Two TreeSEs, the name of the assay, joining method, and the value to
# denote missing values that might occur when object do not share same features, e.g.
# Output: A list of arguments
.merge_TreeSummarizedExperiments <- function(tse_original, tse, join,  
                                             assay_name, missing_values){
    # Get row trees and rownode labels of 1st object
    row_tree1 <- rowTree(tse_original)
    row_node_labs1 <- rowLinks(tse_original)[ , "nodeLab" ]
    if( !is.null(row_node_labs1) ){
        names(row_node_labs1) <- rownames(tse_original)
    }
    # Get row trees and rownode labels of 2nd object
    row_tree2 <- rowTree(tse)
    row_node_labs2 <- rowLinks(tse)[ , "nodeLab" ]
    if( !is.null(row_node_labs2) ){
        names(row_node_labs2) <- rownames(tse)
    }
    
    # Get col trees and column node labels of 1st object
    col_tree1 <- colTree(tse_original)
    col_node_labs1 <- colLinks(tse_original)[ , "nodeLab" ]
    if( !is.null(col_node_labs1) ){
        names(col_node_labs1) <- colnames(tse_original)
    }
    # Get col trees and column node labels of 2nd object
    col_tree2 <- colTree(tse)
    col_node_labs2 <- colLinks(tse)[ , "nodeLab" ]
    if( !is.null(col_node_labs2) ){
        names(col_node_labs2) <- colnames(tse)
    }
    
    # Get reference sequences og both objects
    ref_seqs1 <- referenceSeq(tse_original)
    ref_seqs2 <- referenceSeq(tse)
    
    # Merge data to get a list of arguments
    args <- .merge_SingleCellExperiments(tse_original, tse, join,
                                        assay_name, missing_values)
    
    # rowTree
    # If labels of the 1st tree match with data, add tree1
    if( length(rownames(args$rowData)) > 0 && !is.null( row_tree1 ) &&
        all( rownames(args$rowData) %in% names(row_node_labs1) ) ){
        args$rowTree <- row_tree1
        # Get rowlinks into correct order
        row_node_labs1 <- row_node_labs1[ match(rownames(args$rowData), 
                                        names(row_node_labs1)) ]
        args$rowNodeLab <- row_node_labs1
    # If labels of the 2nd tree match with data, add tree2
    } else if( length(rownames(args$rowData)) > 0 && !is.null( row_tree2 ) &&
               all( rownames(args$rowData) %in% names(row_node_labs2) ) ){
        args$rowTree <- row_tree2
        # Get rowlinks into correct order
        row_node_labs2 <- row_node_labs2[ match(rownames(args$rowData), 
                                        names(row_node_labs2)) ]
        args$rowNodeLab <- row_node_labs2
    }
    # colTree
    # If labels of the 1st tree match with data, add tree1
    if( length(rownames(args$colData)) > 0 && !is.null( col_tree1 ) &&
        all( colnames(args$colData) %in% names(col_node_labs1) ) ){
        args$colTree <- col_tree1
        # Get column links into correct order
        col_node_labs1 <- col_node_labs1[ match(rownames(args$colData), 
                                        names(col_node_labs1)) ]
        args$colNodeLab <- col_node_labs1
    # If labels of the 2nd tree match with data, add tree2
    } else if( length(rownames(args$colData)) > 0  && !is.null( col_tree2 ) &&
               all( colnames(args$colData) %in% names(col_node_labs2) ) ){
        args$colTree <- row_tree2
        # Get column links into correct order
        col_node_labs2 <- col_node_labs2[ match(rownames(args$colData), 
                                        names(col_node_labs2)) ]
        args$colNodeLab <- col_node_labs2
    }
    # Reference sequences
    # If names of 1st sequences match with data, add refseq1
    if( length(rownames(args$rowData)) > 0 && !is.null( ref_seqs1 ) &&
        all( rownames(args$rowData) %in% names(ref_seqs1) ) ){
        # Put sequences into correct order
        ref_seqs1 <- ref_seqs1[ match(rownames(args$rowData), names(ref_seqs1)), ]
        args$referenceSeq <- ref_seqs1
    # If names of 2nd sequences match with data, add refseq2
    } else if( length(rownames(args$rowData)) > 0 && !is.null( ref_seqs2 ) &&
               all( rownames(args$rowData) %in% names(ref_seqs2) ) ){
        # Put sequences into correct order
        ref_seqs2 <- ref_seqs2[ match(rownames(args$rowData), names(ref_seqs2)), ]
        args$referenceSeq <- ref_seqs2
    }
    return(args)
}

######################## .merge_SingleCellExperiments ##########################
# This function merges the data of two SCE objects into one set of arguments that
# can be feed to create a single object.
# Input: Two SCEs
# Output: A single SCE
.merge_SingleCellExperiments <- function(tse_original, tse, join,  
                                         assay_name, missing_values){
    # reducedDim is additional slot for SCE compared to SE. 
    # However, merging reducedDims leads to non-meaningful data
    # Merge data to get a list of arguments
    args <- .merge_SummarizedExperiments(tse_original, tse, join,
                                        assay_name, missing_values)
    return(args)
}

######################## .merge_SummarizedExperiments ##########################
# This function merges the data of two SE objects into one set of arguments that
# can be feed to create a single object.

# Input: Two SEs
# Output: A list of arguments
.merge_SummarizedExperiments <- function(tse_original, tse, join,  
                                         assay_name, missing_values){
    # Merge rowData
    rowdata <- .merge_rowdata(tse_original, tse, join)
    # Merge colData
    coldata <- .merge_coldata(tse_original, tse, join)
    # Merge assay
    assay <- .merge_assay(tse_original, tse, assay_name, join, missing_values, rowdata, coldata)
    assays <- SimpleList(name = assay)
    names(assays) <- assay_name
    # Combine metadata
    metadata <- c( metadata(tse_original), metadata(tse) )
    
    # Create a list of data
    args <- list(assays = assays,
                rowData = rowdata,
                colData = coldata,
                metadata = metadata)
    return(args)
}

################################ .merge_assay ##################################
# This function merges assays.

# Input: Two TreeSEs, the name of the assay, joining method, value to denote
# missing values, merged rowData, and merged colData
# Output: Merged assay
.merge_assay <- function(tse_original, tse, assay_name, join,
                         missing_values, rd, cd){
    # Take assays
    assay1 <- assay(tse_original, assay_name)
    assay2 <- assay(tse, assay_name)
    
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
    # If ranks were found
    if( length(ranks_ind) != 0 ){
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
    }
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
