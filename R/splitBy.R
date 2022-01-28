#' Split \code{SummarizedExperiment} column-wise or row-wise based on grouping variable
#'
#' @param x A
#'   \code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#'   object
#'
#' @param grouping A single character value for selecting the grouping variable
#'   from \code{colData} or \code{rowData}.
#' 
#' @param MARGIN A single numeric value, 1 (row) or 2 (col),  for selecting from 
#'   where grouping variable should be searched.
#'   
#' @param ... 
#' \itemize{
#'   \item{\code{use_names} A boolean value to select whether to name elements of
#'   list by their group names.}
#' }
#'
#' @return
#' List of \code{SummarizedExperiment} objects in \code{SimpleList} format.
#'
#'
#' @details
#' \code{splitBy} split data based on grouping variable. Splitting can be done
#' column-wise or row-wise. Returned value is a list of \code{SummarizedExperiment}
#' objects; each element containing members of each group.
#'
#' @seealso
#' \code{\link[=splitByRanks]{splitByRanks}}
#' \code{\link[=unsplitByRanks]{unsplitByRanks}}
#' \code{\link[=merge-methods]{mergeRows}},
#' \code{\link[scuttle:sumCountsAcrossFeatures]{sumCountsAcrossFeatures}},
#' \code{\link[=agglomerate-methods]{agglomerateByRank}},
#' \code{\link[SingleCellExperiment:altExps]{altExps}},
#' \code{\link[SingleCellExperiment:splitAltExps]{splitAltExps}}
#'
#'
#' @name splitBy
#' @export
#' 
#' @author Leo Lahti and Tuomas Borman. Contact: \url{microbiome.github.io}
#'
#' @examples
#' data(GlobalPatterns)
#' tse <- GlobalPatterns
#' # Split data based on SampleType. 
#' se_list <- splitBy(tse, grouping = "SampleType")
#' 
#' se_list
#' 
#' # Create arbitrary groups
#' colData(tse)$group <- sample(1:10, ncol(tse), replace = TRUE)
#' rowData(tse)$group <- sample(1:10, nrow(tse), replace = TRUE)
#' 
#' # If variable named equally can be found from both colData and rowData, 
#' # MARGIN must be specified
#' se_list <- splitBy(tse, grouping = "group", MARGIN = 2)
#' 
#' # It is possible to split data also in row-wise
#' se_list <- splitBy(tse, grouping = "group", MARGIN = 1)
#' 
#' # Split data based on Phyla
#' se_list <- splitBy(tse, "Phylum")
#' 
#' se_list
#' 
#' # List of SE objects is returned. 
#' # Each element is named based on their group name. If you don't want to name
#' # elements, use use_name = FALSE
#' se_list <- splitBy(tse, grouping = "SampleType", use_name = FALSE)
#' 
#' # If you want to combine groups back together, you can use unsplitBy
#' unsplitBy(se_list)
#' 
NULL

#' @rdname splitBy
#' @export
setGeneric("splitBy",
           signature = "x",
           function(x, grouping = NULL, MARGIN = NULL, ...)
               standardGeneric("splitBy"))

#' @rdname splitBy
#' @export
setMethod("splitBy", signature = c(x = "ANY"),
    function(x, grouping = NULL, MARGIN = NULL, ...){
        ############################## INPUT CHECK #############################
        # Check grouping
        if( !.is_non_empty_string(grouping) ){
            stop("'grouping' must be a single non-empty character value ",
                 "and it must specify variable from colData or rowData.",
                 call. = FALSE)
        }
        # Check MARGIN
        if( !(MARGIN == "row" || MARGIN == "col" || is.null(MARGIN) ||
              (is.numeric(MARGIN) && (MARGIN == 1 || MARGIN == 2))) ){
            stop("'MARGIN' must be 1, 2, or NULL",
                 call. = FALSE)
        }
        ############################ INPUT CHECK END ###########################
        # Should data be splitted row-wise or column-wise?
        # Chek that variable can be found.
        MARGIN <- .split_by_rowise_or_colwise(x, grouping, MARGIN)
        # Split data
        se_list <- .split_by(x, grouping, MARGIN, ...)
        return(se_list)
    }
)

#' @param x A list of
#'   \code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#'   objects
#' 
#' @rdname splitBy
#' @export
setGeneric("unsplitBy",
           signature = c("x"),
           function(x, MARGIN = NULL, ...)
               standardGeneric("unsplitBy"))

#' @rdname splitBy
#' @export
setMethod("unsplitBy", signature = c(x = "list"),
    function(x, MARGIN = NULL, ...){
        ############################## INPUT CHECK #############################
        # Check MARGIN
        if( !(MARGIN == "row" || MARGIN == "col" || is.null(MARGIN) ||
              (is.numeric(MARGIN) && (MARGIN == 1 || MARGIN == 2))) ){
            stop("'MARGIN' must be 1, 2, or NULL",
                 call. = FALSE)
        }
        ############################ INPUT CHECK END ###########################
        # Can data be combined row-wise or column-wise
        # If row_wise, then column names should be equal
        row_wise <- all(sapply(x, FUN = function(obj){
            identical( colnames(obj), colnames(x[[1]])) 
        }))
        # If col-wise, then rownames should be equal
        col_wise <- all(sapply(x, FUN = function(obj){
            identical( rownames(obj), rownames(x[[1]])) 
        }))
        
        # Combine data
        # If data can be combined in both directions, and MARGIN is not specified
        if( row_wise && col_wise && is.null(MARGIN) ){
            stop("Both rownames and colnames match. Specify with 'MARGIN' ",
                 "whether to combine data column-wise or row-wise.",
                 call. = FALSE)
            # If data is combined column-wise
        } else if( col_wise && (MARGIN == "col" || MARGIN == 2 || is.null(MARGIN)) ){
            se <- do.call(cbind, x)
            # If data is combined row-wise
        } else if( row_wise && (MARGIN == "row" || MARGIN == 1 || is.null(MARGIN)) ){
            se <- do.call(rbind, x)
            # If MARGIN is specified but rownames do not match
        } else if ( !col_wise && (MARGIN == "col" || MARGIN == 2) ){
            stop("rownames should match when data is combined column-wise",
                 call. = FALSE)
            # If MARGIN is specified but colnames do not match
        } else if ( !row_wise && (MARGIN == "row" || MARGIN == 1) ){
            stop("colnames should match when data is combined row-wise",
                 call. = FALSE)
            # If rownames nor colnames do not match, data cannot be combined
        } else{
            stop("Data cannot be combined since rownames nor colnames do not match.",
                 call. = FALSE)
        }
        return(se)
    }
)
   
#' @rdname splitBy
#' @export
setMethod("unsplitBy", signature = c(x = "SimpleList"),
    function(x, MARGIN = NULL, ...){
        # Convert into a list
        x <- as.list(x)
        unsplitBy(x, MARGIN, ...)
    }
)
################################ HELP FUNCTIONS ################################
# Split data in column-wise or row-wise based on grouping variable. 
.split_by <- function(x, grouping, MARGIN, use_names = TRUE, ...){
    # Check use_names
    if( !.is_a_bool(use_names) ){
        stop("'use_names' must be a boolean value.", 
             call. = FALSE)
    }
    
    # Split data in column-wise
    if( MARGIN == "col" || MARGIN == 2 ){
        # Get sample indices for each group
        group_indices <- split(seq_len(ncol(x)), colData(x)[[grouping]])
        # Split data
        se_list <- lapply(seq_along(group_indices), function(i){ 
            temp <- x[ , group_indices[[i]] ]
            return(temp)
        })
        # Split data in row-wise
    } else if( MARGIN == "row" || MARGIN == 1 ){
        # Get sample indices for each group
        group_indices <- split(seq_len(nrow(x)), rowData(x)[[grouping]])
        # Split data
        se_list <- lapply(seq_along(group_indices), function(i){ 
            temp <- x[ group_indices[[i]], ]
            return(temp)
        })
    }
    # If TRUE, give group names for elements of list
    if( use_names ){
        names(se_list) <- names(group_indices)
    }
    # Convert into SimpleList
    se_list <- SimpleList(se_list)
    return(se_list)
}

# This function returns the MARGIN: where the grouping variable is found. 
# It also checks, that variable can be found. 
# If user has specified MARGIN, it checks that it is correct.
.split_by_rowise_or_colwise <- function(x, grouping, MARGIN){
    # Check if variable can be found
    if( is.null(MARGIN) ){
        # Is variable from colData or rowData?
        is_coldata_variable <- any( grouping %in% colnames(colData(x)) )
        is_rowdata_variable <- any( grouping %in% colnames(rowData(x)) )
        
        # If variable is can be found from both, but MARGIN is not specified
        if( is_coldata_variable && is_rowdata_variable ){
            stop("'grouping' defines variable from both colData and rowData.",
                 " Use 'MARGIN' to specify the MARGIN.",
                 call. = FALSE)
        }
        # If variable cannot be found
        if( !is_coldata_variable && !is_rowdata_variable ){
            stop("Variable defined by 'grouping' cannot be found.",
                 call. = FALSE)
        }
        # Specify MARGIN
        MARGIN <- ifelse(is_coldata_variable, "col", "row")
    } else if( MARGIN == "col" || MARGIN == 2 ){
        # Can variable be found from colData?
        is_coldata_variable <- any( grouping %in% colnames(colData(x)) )
        is_rowdata_variable <- FALSE
        
        # If variable cannot be found
        if( !is_coldata_variable ){
            stop("Variable defined by 'grouping' cannot be found from colData.",
                 call. = FALSE)
        }
    } else if( MARGIN == "row" || MARGIN == 1 ){
        # Can variable be found from rowData?
        is_rowdata_variable <- any( grouping %in% colnames(rowData(x)) )
        is_coldata_variable <- FALSE
        
        # If variable cannot be found
        if( !is_rowdata_variable ){
            stop("Variable defined by 'grouping' cannot be found from rowData.",
                 call. = FALSE)
        }
    }
    return(MARGIN)
}
