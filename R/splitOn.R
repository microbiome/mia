#' Split \code{TreeSummarizedExperiment} column-wise or row-wise based on grouping variable
#'
#' @param x A
#'   \code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#'   object or a list of 
#'   \code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#'   objects.
#'
#' @param f A single character value for selecting the grouping variable
#'   from \code{rowData} or \code{colData} or a \code{factor} or \code{vector} 
#'   with the same length as one of the dimensions. If \code{f} matches with both
#'   dimensions, \code{MARGIN} must be specified. 
#'   Split by cols is not encouraged, since this is not compatible with 
#'   storing the results in \code{altExps}.
#'
#' @param keep.dimred \code{TRUE} or \code{FALSE}: Should the
#'   \code{reducedDims(x)} be transferred to the result? Please note, that this
#'   breaks the link between the data used to calculate the reduced dims.
#'   (By default: \code{keep.dimred = FALSE})
#' 
#' @param keep_reducedDims Deprecated. Use \code{keep.dimred} instead.
#'   
#' @param update.tree \code{TRUE} or \code{FALSE}: Should the rowTree be updated
#'   based on splitted data? Option is enabled when \code{x} is a 
#'   \code{TreeSummarizedExperiment} object or a list of such objects. 
#'   (By default: \code{update.tree = FALSE})
#' 
#' @param update_rowTree Deprecated. Use \code{update.tree } instead.
#'   
#' @param altexp a \code{character} vector specifying the alternative experiments
#'   to be unsplit. (By default: \code{altexp = names(altExps(x))})
#' 
#' @param altExpNames Deprecated. Use \code{altexp} instead.
#'   
#' @param ... Arguments passed to \code{agglomerateByVariable} function for
#'   \code{SummarizedExperiment} objects and other functions.
#'   See \code{\link[=agglomerate-methods]{agglomerateByVariable}} for more 
#'   details.
#'   \itemize{
#'     \item{\code{use.names} A single boolean value to select whether to name elements of
#'     list by their group names.}
#'   }
#'
#'
#' @details
#' \code{splitOn} split data based on grouping variable. Splitting can be done
#' column-wise or row-wise. The returned value is a list of
#' \code{SummarizedExperiment} objects; each element containing members of each
#' group.
#'
#' @return
#' For \code{splitOn}: \code{SummarizedExperiment} objects in a \code{SimpleList}.
#'
#' For \code{unsplitOn}: \code{x}, with \code{rowData} and \code{assay}
#' data replaced by the unsplit data. \code{colData} of x is kept as well
#' and any existing \code{rowTree} is dropped as well, since existing
#' \code{rowLinks} are not valid anymore.
#' 
#' @name splitOn
#' @seealso
#' \code{\link[=agglomerate-methods]{agglomerateByRanks}}
#' \code{\link[=agglomerate-methods]{agglomerateByVariable}},
#' \code{\link[scuttle:sumCountsAcrossFeatures]{sumCountsAcrossFeatures}},
#' \code{\link[=agglomerate-methods]{agglomerateByRank}},
#' \code{\link[SingleCellExperiment:altExps]{altExps}},
#' \code{\link[SingleCellExperiment:splitAltExps]{splitAltExps}}
#'
#' @export
#' @author Leo Lahti and Tuomas Borman. Contact: \url{microbiome.github.io}
#' 
#' @examples
#' data(GlobalPatterns)
#' tse <- GlobalPatterns
#' # Split data based on SampleType. 
#' se_list <- splitOn(tse, f = "SampleType")
#' 
#' # List of SE objects is returned. 
#' se_list
#' 
#' # Create arbitrary groups
#' rowData(tse)$group <- sample(1:3, nrow(tse), replace = TRUE)
#' colData(tse)$group <- sample(1:3, ncol(tse), replace = TRUE)
#' 
#' # Split based on rows
#' # Each element is named based on their group name. If you don't want to name
#' # elements, use use_name = FALSE. Since "group" can be found from rowdata and colData
#' # you must use MARGIN.
#' se_list <- splitOn(tse, f = "group", use.names = FALSE, MARGIN = 1)
#' 
#' # When column names are shared between elements, you can store the list to altExps
#' altExps(tse) <- se_list
#' 
#' altExps(tse)
#' 
#' # If you want to split on columns and update rowTree, you can do
#' se_list <- splitOn(tse, f = colData(tse)$group, update.tree = TRUE)
#' 
#' # If you want to combine groups back together, you can use unsplitBy
#' unsplitOn(se_list)
#' 
NULL

#' @rdname splitOn
#' @export
setGeneric("splitOn",
            signature = "x",
            function(x, ...)
                standardGeneric("splitOn"))

# This function collects f (grouping variable), MARGIN, and 
# use.names and returns them as a list.
.norm_args_for_split_by <- function(x, f, MARGIN = NULL, use.names = use_names,
                                    use_names = TRUE, ...){
    # input check
    # Check f
    if(is.null(f)){
        stop("'f' must either be a single non-empty character value or",
            " vector coercible to factor alongside the one of the dimensions of 'x'",
            call. = FALSE)
    }
    # Check MARGIN
    if( !(is.null(MARGIN) || (is.numeric(MARGIN) && (MARGIN == 1 || MARGIN == 2 ))) ){
        stop("'MARGIN' must be NULL, 1, or 2.", call. = FALSE )
    }
    # If f is a vector containing levels
    if( !.is_non_empty_string(f) ){
        # Convert into factors
        f <- factor(f, unique(f))
        # Check if the length of f matches with one of the dimensions
        if(!length(f) %in% dim(x)){
            stop("'f' must either be a single non-empty character value or",
                " vector coercible to factor alongside the on of the ",
                "dimensions of 'x'.",
                call. = FALSE)
        # If it matches with both dimensions, give error if MARGIN is not specified
        } else if( is.null(MARGIN) && all(length(f) == dim(x)) ){
            stop("The length of 'f' matches with nrow and ncol. ",
                "Please specify 'MARGIN'.", call. = FALSE)
        # If MARGIN is specified but it does not match with length of f
        } else if( !is.null(MARGIN) && (length(f) !=  dim(x)[[MARGIN]]) ){
            stop("'f' does not match with ", 
                ifelse(MARGIN==1, "nrow", "ncol"), ". Please check 'MARGIN'.",
                call. = FALSE)
        # IF f matches with nrow
        } else if(length(f) == dim(x)[[1]] && is.null(MARGIN)  ){
            MARGIN <- 1L
        # If f matches with ncol
        } else if( is.null(MARGIN) ){
            MARGIN <- 2L
        }
    # Else if f is a character specifying column from rowData or colData  
    } else {
        # If MARGIN is specified
        if( !is.null(MARGIN) ){
            # Search from rowData or colData based on MARGIN
            dim_name <- switch(MARGIN,
                                "1" = "rowData",
                                "2" = "colData")
            # Specify right function
            dim_FUN <- switch(MARGIN,
                                "1" = retrieveFeatureInfo,
                                "2" = retrieveCellInfo)
            # Try to get information
            tmp <- try({dim_FUN(x, f, search = dim_name)},
                        silent = TRUE)
            # Give error if it cannot be found
            if(is(tmp,"try-error")){
                stop("'f' is not found. ",
                    "Please check that 'f' specifies a column from ", dim_name, ".", 
                    call. = FALSE)
            }
            # Get values
            f <- tmp$value
        # Else if MARGIN is not specified
        } else{
            # Try to get information from rowData
            tmp_row <- try({retrieveFeatureInfo(x, f, search = "rowData")},
                            silent = TRUE)
            # Try to get information from colData
            tmp_col <- try({retrieveCellInfo(x, f, search = "colData")}, 
                            silent = TRUE)
            
            # If it was not found 
            if( is(tmp_row, "try-error") && is(tmp_col, "try-error") ){
                stop("'f' is not found. ",
                    "Please check that 'f' specifies a column from ",
                    "rowData or colData.", 
                    call. = FALSE)
                # If f was found from both
            } else if( !is(tmp_row, "try-error") && !is(tmp_col, "try-error") ){
                stop("'f' can be found from both rowData and colData. ",
                    "Please specify 'MARGIN'.",
                    call. = FALSE)
                # If it was found from rowData
            } else if( !is(tmp_row, "try-error") ){
                MARGIN <- 1L
                # Get values
                f <- tmp_row$value
                # Otherwise, it was found from colData
            } else{
                MARGIN <- 2L
                # Get values
                f <- tmp_col$value
            }
        }
        # Convert values into factors
        f <- factor(f, unique(f))
        
        # If there are NAs, add NA as level
        if( any(is.na(f)) ){
            f <- addNA(f)
        }
    }
    # Check use.names
    if( !.is_a_bool(use.names) ){
        stop("'use.names' must be TRUE or FALSE.",
            call. = FALSE)
    }
    # Create a list from arguments
    list(f = f,
        MARGIN = MARGIN,
        use.names = use.names)
}

# PErform the split
.split_on <- function(x, args, ...){
    # Get grouping variable and its values
    f <- args[["f"]]
    # Choose nrow or ncol based on MARGIN
    dim_FUN <- switch(args[["MARGIN"]],
                        "1" = nrow,
                        "2" = ncol)
    # Get indices from 1 to nrow/ncol
    idx <- seq_len(dim_FUN(x))
    # Split indices into groups based on grouping variable
    idxs <- split(idx, f)
    # Subset function takes SE and list of groups which have indices
    # It divides the data into groups
    subset_FUN <- function(x, i = TRUE, j = TRUE){
        x[i, j]
    }
    # Based on MARGIN, divide data in row-wise or column-wise
    if(args[["MARGIN"]] == 1){
        ans <- SimpleList(lapply(idxs, subset_FUN, x = x))
    } else {
        ans <- SimpleList(lapply(idxs, subset_FUN, x = x, i = TRUE))
    }
    # If user do not want to use names, unname
    if(!args[["use.names"]]){
        ans <- unname(ans)
    # Otherwise convert NAs to "NA", if there is a level that do not have name
    } else{
        names(ans)[ is.na(names(ans)) ] <- "NA"
    }
    ans
}

#' @rdname splitOn
#' @export
setMethod("splitOn", signature = c(x = "SummarizedExperiment"),
    function(x, f = NULL,  ...){
        # Get arguments
        args <- .norm_args_for_split_by(x, f = f, ...)
        # Split data
        .split_on(x, args, ...)
    }
)

#' @rdname splitOn
#' @export
setMethod("splitOn", signature = c(x = "SingleCellExperiment"),
    function(x, f = NULL, ...){
        # Get arguments
        args <- .norm_args_for_split_by(x, f = f, ...)
        # Should alternative experiment be removed? --> yes
        args[["altexp.rm"]] <- TRUE
        # Split data
        .split_on(x, args, ...)
    }
)

#' @rdname splitOn
#' @export
setMethod("splitOn", signature = c(x = "TreeSummarizedExperiment"),
    function(x, f = NULL, update.tree = update_rowTree, update_rowTree = FALSE,
            ...){
        # Input check
        # Check update.tree
        if( !.is_a_bool(update.tree) ){
            stop("'update.tree' must be TRUE or FALSE.",
                call. = FALSE)
        }
        # Input check end
        # Split data
        x <- callNextMethod()
        # Manipulate rowTree or not?
        if( update.tree ){
            # If the returned value is a list, go through all of them
            if( is(x, 'SimpleList') ){
                x <- SimpleList(lapply(x, .agglomerate_trees))

            } else {
                # Otherwise, the returned value is TreeSE
                x <- .agglomerate_trees(x)
            }
        }
        x
    }
)

################################################################################
# unsplitOn

#' @rdname splitOn
#' @export
setGeneric("unsplitOn",
            signature = c("x"),
            function(x, ...)
                standardGeneric("unsplitOn"))

# Perform the unsplit
.list_unsplit_on <- function(ses, update.tree = FALSE, MARGIN = NULL, ...){
    # Input check
    is_check <- vapply(ses,is,logical(1L),"SummarizedExperiment")
    if(!all(is_check)){
        stop("Input must be a list of SummarizedExperiment or derived objects ",
            "only.",
            call. = FALSE)
    }
    # Check update.tree
    if( !.is_a_bool(update.tree) ){
        stop("'update.tree' must be TRUE or FALSE.",
            call. = FALSE)
    }
    if( !(is.null(MARGIN) || (is.numeric(MARGIN) && (MARGIN == 1 || MARGIN == 2 ))) ){
        stop("'MARGIN' must be NULL, 1, or 2.", call. = FALSE )
    }
    # Input check end
    # If list contains only one element, return it
    if( length(ses) == 1 ){
        return(ses[[1]])
    }
    # Get dimensions of each SE in the list
    dims <- vapply(ses, dim, integer(2L))
    # Based on which dimension SE objects share, select MARGIN.
    # If they share rows, then MARGIN is col, and vice versa
    if( is.null(MARGIN) ){
        if( length(unique(dims[1L,])) == 1 && length(unique(dims[2L,])) == 1 ){
            stop("The dimensions match with row and column-wise. ",
                "Please specify 'MARGIN'.", call. = FALSE)
        } else if(length(unique(dims[1L,])) == 1L){
            MARGIN <- 2L
        } else if(length(unique(dims[2L,])) == 1L) {
            MARGIN <- 1L
        } else {
            stop("The dimensions are not equal across all elements. ", 
                "Please check that either number of rows or columns match.", 
                call. = FALSE)
        }
    } else{
        # Get correct dimension, it is opposite of MARGIN
        dim <- ifelse(MARGIN == 1, 2, 1)
        if( length(unique(dims[dim,])) != 1L ){
            stop("The dimensions are not equal across all elements.", call. = FALSE)
        }
    }
    
    # Get the class of objects SCE, SE or TreeSE
    class_x <- class(ses[[1L]])
    # Combine assays
    args <- list(assays = .unsplit_assays(ses, MARGIN = MARGIN))
    # Combine rowData if data share columns
    if(MARGIN == 1L){
        rd <- .combine_rowData(ses)
        # Add rownames since they are missing after using combining
        rownames(rd) <- unlist(unname(lapply(ses, rownames)))
        rr <- .combine_rowRanges(ses)
        args$rowRanges <- rr
        args$colData <- colData(ses[[1L]])
    # Combine colData if data share rows
    } else {
        args$colData <- .combine_colData(ses)
        args$rowRanges <- rowRanges(ses[[1L]])
        rd <- rowData(ses[[1L]])
    }
    # Create a object specified by class_x from the data
    ans <- do.call(class_x, args)
    # Add rowData
    rowData(ans) <- rd
    # Update rownames
    rownames(ans) <- rownames(rd)
    
    # IF the object is TreeSE. add rowTree
    if( class_x == "TreeSummarizedExperiment" ){
        # Update or add old tree from the first element of list
        if( update.tree ){
            ans <- addHierarchyTree(ans)
        } else{
            rowTree(ans) <- rowTree(ses[[1L]])
        }
    }
    ans
}

#' @importFrom SummarizedExperiment colData
#' @importFrom BiocGenerics rbind
.combine_colData <- function(ses) {
    # Get colDatas of objects
    cds <- lapply(ses, colData)
    # Bind them together row-wise
    cd <- do.call(rbind,unname(cds))
    # Add sample names
    rownames(cd) <- unlist(unname(lapply(ses, colnames)))
    cd
}


#' @rdname splitOn
#' @importFrom SingleCellExperiment altExpNames altExp altExps
#' @export
setMethod("unsplitOn", signature = c(x = "list"),
    function(x, update.tree = update_rowTree, update_rowTree = FALSE, ...){
        # Unsplit list and create SCE, SE, or TreeSE from it
        .list_unsplit_on(x, update.tree, ...)
    }
)
#' @rdname splitOn
#' @importFrom SingleCellExperiment altExpNames altExp altExps
#' @export
setMethod("unsplitOn", signature = c(x = "SimpleList"),
    function(x, update.tree = update_rowTree, update_rowTree = FALSE, ...){
        unsplitOn(as.list(x), update.tree, ...)
    }
)

#' @rdname splitOn
#' @importFrom SingleCellExperiment altExpNames altExp altExps reducedDims<-
#' @export
setMethod("unsplitOn", signature = c(x = "SingleCellExperiment"),
    function(x, altexp = altExpNames, altExpNames = names(altExps(x)),
            keep.dimred = keep_reducedDims,
            keep_reducedDims = FALSE, ...){
        # input check
        if(!.is_a_bool(keep.dimred)){
            stop("'keep.dimred' must be TRUE or FALSE.", call. = FALSE)
        }
        # Get alternative experiment names since data is located there
        ae_names <- names(altExps(x))
        # Get only those experiments that user has specified
        ae_names <- ae_names[ae_names %in% altexp]
        if(length(ae_names) == 0L){
            stop("No altExp matching 'altexp' in name.", call. = FALSE)
        }
        # Get alternative experiments as a list
        ses <- altExps(x)[ae_names]
        # And unsplit the data
        se <- .list_unsplit_on(ses, ...)
        # Add reducedDims if specified
        if( keep.dimred ){
            reducedDims(se) <- reducedDims(x)
        }
        return(se)
    }
)
