#' Split \code{SummarizedExperiment} column-wise or row-wise based on grouping variable
#'
#' @param x A
#'   \code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#'   object
#'
#' @param grouping A single character value for selecting the grouping variable
#'   from \code{colData} or \code{rowData}.
#' 
#' @param MARGIN A single character or numeric value for selecting from where grouping
#'   variable should be searched. Must be "row"/1 or "col"/2. 
#'   
#' @param ... 
#' \itemize{
#'   \item{\code{use_names} A boolean value to select whether to name elements of
#'   list by their group names.}
#' }
#'
#' @return
#' List of \code{SummarizedExperiment} objects.
#'
#'
#' @details
#' \code{splitBy} split data based on grouping variable. Splitting can be done
#' column-wise or row-wise. Returned value is a list of \code{SummarizedExperiment}
#' objects; each element containing members of each group.
#'
#' @seealso
#' \code{\link[mia::splitByRanks]{splitByRanks}},
#' \code{\link[mia::unSplitByRanks]{unSplitByRanks}},
#' \code{\link[=merge-methods]{mergeRows}},
#' \code{\link[scuttle:sumCountsAcrossFeatures]{sumCountsAcrossFeatures}},
#' \code{\link[=agglomerate-methods]{agglomerateByRank}},
#' \code{\link[SingleCellExperiment:altExps]{altExps}},
#' \code{\link[SingleCellExperiment:splitAltExps]{splitAltExps}}
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
#' se_list[1:2]
#' 
#' # Create arbitrary groups
#' colData(tse)$group <- sample(1:10, ncol(tse), replace = TRUE)
#' rowData(tse)$group <- sample(1:10, nrow(tse), replace = TRUE)
#' 
#' # If variable named equally can be found from both colData and rowData, 
#' # MARGIN must be specified
#' se_list <- splitBy(tse, grouping = "SampleType", MARGIN = "col")
#' 
#' # It is possible to split data also in row-wise
#' se_list <- splitBy(tse, grouping = "SampleType", MARGIN = "row")
#' 
#' # However, if you want to split data based on ranks, use splitByRanks
#' se_list <- splitByRanks(tse)
#' 
#' # List of SE objects is returned. 
#' # Each element is named based on their group name. If you don't want to name
#' # elements, use use_name = FALSE
#' se_list <- splitBy(tse, grouping = "SampleType", use_name = FALSE)
#' 
#' # If you want to combine groups back together, you can use cbind or rbind
#' do.call(cbind, se_list)
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
            stop("'MARGIN' must be 'row', 1, 'col', 2, or NULL",
                 call. = FALSE)
        }
        ############################ INPUT CHECK END ###########################
        MARGIN <- .split_by_rowise_or_colwise(x, grouping, MARGIN)
        se_list <- .split_by(x, grouping, MARGIN, ...)
        return(se_list)
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
        group_indices <- split(1:ncol(x), colData(x)[[grouping]])
        # Split data
        se_list <- lapply(seq_along(group_indices), function(i){ 
            temp <- x[ , group_indices[[i]] ]
            return(temp)
        })
        # Split data in row-wise
    } else if( MARGIN == "row" || MARGIN == 1 ){
        # Get sample indices for each group
        group_indices <- split(1:nrow(x), rowData(x)[[grouping]])
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
    return(se_list)
}

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
    # If variable is taxonomy rank
    if( is_rowdata_variable && any(grouping %in% taxonomyRanks(x)) ){
        stop("'grouping' defines a taxonomy rank. Please use splitByRank instead.",
             call. = FALSE)
    }
    return(MARGIN)
}
