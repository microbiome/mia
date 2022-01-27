#' Split/Unsplit a \code{SingleCellExperiment} by grouping variable.s
#'
#'
#' @param x a
#'   \code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#'   object
#'
#' @param ... arguments passed to \code{agglomerateByRank} function for
#'   \code{SummarizedExperiment} objects and other functions.
#'   See \code{\link[=agglomerate-methods]{agglomerateByRank}} for more details.
#'
#' @return
#' For \code{splitBy}: \code{x}, with objects of \code{x} agglomerated for
#' selected ranks as \code{altExps}.
#'
#'
#' @details
#' \code{splitBy} will use by default ...
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
#'
#' @examples
#' data(GlobalPatterns)
NULL


#' @rdname splitBy
#' @export
setGeneric("splitBy",
           signature = "x",
           function(x, grouping = NULL, direction = NULL, ...)
               standardGeneric("splitBy"))



#' @rdname splitBy
#' @export
setMethod("splitBy", signature = c(x = "any"),
    function(x, grouping = NULL, direction = NULL, ...){
        ############################## INPUT CHECK #############################
        # Check grouping
        if( !.is_non_empty_string(grouping) ){
            stop("'grouping' must be a single non-empty character value ",
                 "and it must specify variable from colData or rowData.",
                 call. = FALSE)
        }
        # Check direction
        if( !(direction == "row" || direction == "col" || is.null(direction)) ){
            stop("'direction' must be 'row', 'col', or NULL",
                 call. = FALSE)
        }
        
        # Check if variable can be found
        if( direction == "col" ){
            # Can variable be found from colData?
            is_coldata_variable <- any( grouping %in% colnames(colData(x)) )
            is_rowdata_variable <- FALSE

            # If variable cannot be found
            if( !is_coldata_variable ){
                stop("Variable defined by 'grouping' cannot be found from colData.",
                     call. = FALSE)
            }
        } else if( direction == "row" ){
            # Can variable be found from rowData?
            is_rowdata_variable <- any( grouping %in% colnames(rowData(x)) )
            is_coldata_variable <- FALSE
            
            # If variable cannot be found
            if( !is_rowdata_variable ){
                stop("Variable defined by 'grouping' cannot be found from rowData.",
                call. = FALSE)
            }
        } else{
            # Is variable from colData or rowData?
            is_coldata_variable <- any( grouping %in% colnames(colData(x)) )
            is_rowdata_variable <- any( grouping %in% colnames(rowData(x)) )
            
            # If variable is can be found from both, but direction is not specified
            if( is_coldata_variable && is_rowdata_variable ){
                stop("'grouping' defines variable from both colData and rowData.",
                     " Use 'direction' to specify the direction.",
                     call. = FALSE)
            }
            # If variable cannot be found
            if( !is_coldata_variable && !is_rowdata_variable ){
                stop("Variable defined by 'grouping' cannot be found.",
                     call. = FALSE)
            }
        }
        # If variable is taxonomy rank
        if( is_rowdata_variable && any(grouping %in% taxonomyRanks(x)) ){
            stop("'grouping' defines a taxonomy rank. Please use splitByRank instead.",
                 call. = FALSE)
        }
        ############################ INPUT CHECK END ###########################
        
        
        
    }
)

