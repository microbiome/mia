#' Split \code{SummarizedExperiment} column-wise or row-wise based on grouping variable
#'
#' @param x A
#'   \code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#'   object or a list of 
#'   \code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#'   objects.
#'
#' @param f A single character value for selecting the grouping variable
#'   from \code{rowData} or \code{colData} or a \code{factor} or \code{vector} 
#'   with the same length as one of the dimensions. Rows take precedence.
#'   Split by cols is not encouraged, since this is not compatible with 
#'   storing the results in \code{altExps}.
#'
#' @param keep_reducedDims \code{TRUE} or \code{FALSE}: Should the
#'   \code{reducedDims(x)} be transferred to the result? Please note, that this
#'   breaks the link between the data used to calculate the reduced dims.
#'   (By default: \code{keep_reducedDims = FALSE})
#'   
#' @param update_rowTree \code{TRUE} or \code{FALSE}: Should the rowTree be updated
#'   based on splitted data? Option is enabled when \code{x} is a 
#'   \code{TreeSummarizedExperiment} object or a list of such objects. 
#'   (By default: \code{update_rowTree = FALSE})
#'   
#' @param altExpNames a \code{character} vector specifying the alternative experiments
#'   to be unsplit. (By default: \code{altExpNames = altExpNames(x)})
#'   
#' @param ... Arguments passed to \code{mergeRows}/\code{mergeCols} function for
#'   \code{SummarizedExperiment} objects and other functions.
#'   See \code{\link[=agglomerate-methods]{mergeRows}} for more details.
#'   \itemize{
#'     \item{\code{use_names} A single boolean value to select whether to name elements of
#'     list by their group names.}
#'   }
#'
#' @return
#' For \code{splitBy}: \code{SummarizedExperiment} objects in a \code{SimpleList}.
#'
#' For \code{unsplitBy}: \code{x}, with \code{rowData} and \code{assay}
#' data replaced by the unsplit data. \code{colData} of x is kept as well
#' and any existing \code{rowTree} is dropped as well, since existing
#' \code{rowLinks} are not valid anymore.
#'
#' @details
#' \code{splitBy} split data based on grouping variable. Splitting can be done
#' column-wise or row-wise. The returned value is a list of
#' \code{SummarizedExperiment} objects; each element containing members of each
#' group.
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
#' @name splitOn
#' @export
#'
#' @examples
#' data(GlobalPatterns)
#' tse <- GlobalPatterns
#' # Split data based on SampleType. 
#' se_list <- splitOnCols(tse, f = "SampleType")
#' 
#' # List of SE objects is returned. 
#' se_list
#' 
#' # Create arbitrary groups
#' rowData(tse)$group <- sample(1:10, nrow(tse), replace = TRUE)
#' 
#' # Split based on rows
#' # Each element is named based on their group name. If you don't want to name
#' # elements, use use_name = FALSE
#' se_list <- splitOnRows(tse, f = "groyup")
#' se_list
#' 
#' # If you want to combine groups back together, you can use unsplitBy
#' unsplitOn(se_list)
#' 
NULL

#' @rdname splitOn
#' @export
setGeneric("splitOnRows",
           signature = "x",
           function(x, ...)
               standardGeneric("splitOnRows"))

#' @rdname splitOn
#' @export
setGeneric("splitOnCols",
           signature = "x",
           function(x, ...)
               standardGeneric("splitOnCols"))

# This function collects f (grouping variable), MARGIN, and 
# use_names and returns them as a list.
.norm_args_for_split_by <- function(x, f, MARGIN, use_names = TRUE, ...){
    # input check
    if(is.null(f)){
        stop("'f' must either be a single non-empty character value or",
             " vector coercible to factor alongside the specified MARGIN of 'x'.",
             "dimensions of 'x'",
             call. = FALSE)
    }
    # Check f or extract the factor from rowData or colData
    if( !.is_non_empty_string(f) ){
        if( length(f) > 1L ){
            f <- factor(f, unique(f))
        }
        # Check if the length of f matches with one of the dimensions
        if(!length(f) %in% dim(x)[[MARGIN]]){
            stop("'f' must either be a single non-empty character value or",
                 " vector coercible to factor alongside the specified MARGIN of 'x'.",
                 call. = FALSE)
        }
    } else {
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
        # Give error if it cannot be found from neither
        if(is(tmp,"try-error")){
            stop("'f' is not found. . ",
                 "Please check that 'f' specifies a column from ", dim_name, ".", 
                 call. = FALSE)
        } 
        # Get values and convert them into factors
        f <- tmp$value
        f <- factor(f, unique(f))
    }
    # Check use_names
    if( !.is_a_bool(use_names) ){
        stop("'use_names' must be TRUE or FALSE.",
             call. = FALSE)
    }
    # Create a list from arguments
    list(f = f,
         MARGIN = MARGIN,
         use_names = use_names)
}

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
    if(!args[["use_names"]]){
        ans <- unname(ans)
    }
    ans
}

#' @rdname splitOnRows
#' @export
setMethod("splitOnRows", signature = c(x = "SummarizedExperiment"),
    function(x, f = NULL,  ...){
        # Get arguments
        args <- .norm_args_for_split_by(x, f = f, MARGIN = 1)
        # Split data
        .split_on(x, args, ...)
    }
)

#' @rdname splitOnRows
#' @export
setMethod("splitOnRows", signature = c(x = "SingleCellExperiment"),
    function(x, f = NULL, ...){
        # Get arguments
        args <- .norm_args_for_split_by(x, f = f, MARGIN = 1)
        # Should alternative experiment be removed? --> yes
        args[["strip_altexp"]] <- TRUE
        # Split data
        .split_on(x, args, ...)
    }
)

#' @rdname splitonRows
#' @export
setMethod("splitOnRows", signature = c(x = "TreeSummarizedExperiment"),
    function(x, f = NULL, update_rowTree = FALSE,
             ...){
        # Input check
        # Check update_rowTree
        if( !.is_a_bool(update_rowTree) ){
            stop("'update_rowTree' must be TRUE or FALSE.",
                 call. = FALSE)
        }
        # Input check end
        # Split data
        x <- callNextMethod()
        # Manipulate rowTree or not?
        if( update_rowTree ){
            # If the returned value is a list, go through all of them
            if( class(x) == "SimpleList" ){
                x <- SimpleList(lapply(x, addTaxonomyTree))
            } else {
                # Otherwise, the returned value is TreeSE
                x <- addTaxonomyTree(x)
            }
        }
        x
    }
)

#' @rdname splitOn
#' @export
setMethod("splitOnCols", signature = c(x = "SummarizedExperiment"),
    function(x, f = NULL,  ...){
        # Get arguments
        args <- .norm_args_for_split_by(x, f = f, MARGIN = 2)
        # Split data
        .split_on(x, args, ...)
        }
)

#' @rdname splitOn
#' @export
setMethod("splitOnCols", signature = c(x = "SingleCellExperiment"),
    function(x, f = NULL, ...){
        # Get arguments
        args <- .norm_args_for_split_by(x, f = f, MARGIN = 2)
        # Should alternative experiment be removed? --> yes
        args[["strip_altexp"]] <- TRUE
        # Split data
        .split_on(x, args, ...)
    }
)

################################################################################
# unsplitBy

#' @rdname splitOn
#' @export
setGeneric("unsplitBy",
           signature = c("x"),
           function(x, ...)
               standardGeneric("unsplitBy"))

.list_unsplit_by <- function(ses, update_rowTree, ...){
    # Input check
    is_check <- vapply(ses,is,logical(1L),"SummarizedExperiment")
    if(!all(is_check)){
        stop("Input must be a list of SummarizedExperiment or derived objects ",
             "only.",
             call. = FALSE)
    }
    # Check update_rowTree
    if( !.is_a_bool(update_rowTree) ){
        stop("'update_rowTree' must be TRUE or FALSE.",
             call. = FALSE)
    }
    # Input check end
    # Get dimensions of each SE in the list
    dims <- vapply(ses, dim, integer(2L))
    # Based on which dimension SE objects share, select MARGIN.
    # If they share rows, then MARGIN is col, and vice versa
    if(length(unique(dims[1L,])) == 1L){
        MARGIN <- 2L
    } else if(length(unique(dims[2L,])) == 1L) {
        MARGIN <- 1L
    } else {
        stop("No dimensions are equal across all elmenents.", call. = FALSE)
    }
    # Get the class of objects SCE, SE or TreeSE
    class_x <- class(ses[[1L]])
    # Combine assays
    args <- list(assays = .unsplit_assays(ses, MARGIN = MARGIN))
    # Combine rowData if data share columns, and vice versa
    if(MARGIN == 1L){
        rd <- .combine_rowData(ses)
        rr <- .combine_rowRanges(ses)
        args$rowRanges <- rr
    } else {
        args$colData <- .combine_colData(ses)
        args$rowRanges <- rowRanges(ses[[1L]])
        rd <- rowData(ses[[1L]])
    }
    # Create a object specified by class_x from the data
    ans <- do.call(class_x, args)
    # Add rowData
    rowData(ans) <- rd
    if( class_x == "TreeSummarizedExperiment" && update_rowTree ){
        ans <- addTaxonomyTree(ans)
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


#' @rdname splitBy
#' @importFrom SingleCellExperiment altExpNames altExp altExps
#' @export
setMethod("unsplitBy", signature = c(x = "list"),
    function(x, update_rowTree = FALSE, ...){
        # Unsplit list and create SCE, SE, or TreeSE from it
        .list_unsplit_by(x, update_rowTree, ...)
    }
)
#' @rdname splitBy
#' @importFrom SingleCellExperiment altExpNames altExp altExps
#' @export
setMethod("unsplitBy", signature = c(x = "SimpleList"),
    function(x, update_rowTree = FALSE, ...){
        unsplitBy(as.list(x), update_rowTree, ...)
    }
)

#' @rdname splitBy
#' @importFrom SingleCellExperiment altExpNames altExp altExps
#' @export
setMethod("unsplitBy", signature = c(x = "SingleCellExperiment"),
    function(x, altExpNames = altExpNames(x), keep_reducedDims = FALSE, ...){
        # input check
        if(!.is_a_bool(keep_reducedDims)){
            stop("'keep_reducedDims' must be TRUE or FALSE.", call. = FALSE)
        }
        # Get alternative experiment names since data is located there
        ae_names <- altExpNames(x)
        # Get only those experiments that user has specified
        ae_names <- ae_names[ae_names %in% altExpNames]
        if(length(ae_names) == 0L){
            stop("No altExp matching 'altExpNames' in name.", call. = FALSE)
        }
        # Get alternative experiments as a list
        ses <- altExps(x)[ae_names]
        # And unsplit the data
        .unsplit_by(x, ses, keep_reducedDims, ...)
    }
)
   
#' @rdname splitBy
#' @export
setMethod("unsplitBy", signature = c(x = "TreeSummarizedExperiment"),
    function(x, altExpNames = altExpNames(x), keep_reducedDims = FALSE, 
             update_rowTree = FALSE, ...){
        # input check
        if(!.is_a_bool(update_rowTree)){
            stop("'update_rowTree' must be TRUE or FALSE.", call. = FALSE)
        }
        ans <- callNextMethod()
        #
        if( update_rowTree ){
            ans <- addTaxonomyTree(ans)
        }
        ans
    }
)
