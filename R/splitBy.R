#' Split \code{SummarizedExperiment} column-wise or row-wise based on grouping variable
#'
#' @param x A
#'   \code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#'   object
#'
#' @param f A single character value for selecting the grouping variable
#'   from \code{rowData} or \code{colData} or a \code{factor} or \code{vector} 
#'   with the same length as one of the dimensions. Rows take precedence.
#'   Split by cols is not encouraged, since this is not compatible with 
#'   storing the results in \code{altExps}.
#'   
#' @param skip_agglomerate \code{TRUE} or \code{FALSE}: Should the agglomerate
#'   be skipped? (By default: \code{skip_agglomerate = FALSE})
#'
#' @param keep_reducedDims \code{TRUE} or \code{FALSE}: Should the
#'   \code{reducedDims(x)} be transferred to the result? Please note, that this
#'   breaks the link between the data used to calculate the reduced dims.
#'   (default: \code{keep_reducedDims = FALSE})
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
#' For \code{splitBy}: if \code{skip_agglomerate == FALSE} a
#' \code{SummarizedExperiment}. Otherwise \code{SummarizedExperiment} objects in
#' a \code{SimpleList}.
#'
#' For \code{unsplitBy}: \code{x}, with \code{rowData} and \code{assay}
#' data replaced by the unsplit data. \code{colData} of x is kept as well
#' and any existing \code{rowTree} is dropped as well, since existing
#' \code{rowLinks} are not valid anymore.
#'
#' @details
#' \code{splitBy} split data based on grouping variable. Splitting can be done
#' column-wise or row-wise. You can specify if you want to agglomerate the data
#' or not.
#' 
#' If data is agglomerated the returned value is a \code{SummarizedExperiment}
#' object that has values ummed-up based on grouping variable.
#' 
#' If data is not agglomerated, the returned value is a list of
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
#' @name splitBy
#' @export
#'
#' @examples
#' data(GlobalPatterns)
#' tse <- GlobalPatterns
#' # Split data based on SampleType. 
#' se_list <- splitBy(tse, f = "SampleType", skip_agglomerate = TRUE)
#' 
#' se_list
#' 
#' # You can also agglomerate (sum-up) data based on grouping variable
#' se <- splitBy(tse, grouping = "SampleType")
#' 
#' # Create arbitrary groups
#' colData(tse)$group <- sample(1:10, ncol(tse), replace = TRUE)
#' rowData(tse)$group <- sample(1:10, nrow(tse), replace = TRUE)
#' 
#' # If variable named equally can be found from both colData and rowData, 
#' # row takes precedence
#' se_list <- splitBy(tse, f = "group", skip_agglomerate = TRUE)
#' se_list
#' 
#' # List of SE objects is returned. 
#' # Each element is named based on their group name. If you don't want to name
#' # elements, use use_name = FALSE
#' se_list <- splitBy(tse, f = "SampleType", use_names = FALSE,
#'                    skip_agglomerate = TRUE)
#' 
#' # If you want to combine groups back together, you can use unsplitBy
#' unsplitBy(se_list)
#' 
NULL

#' @rdname splitBy
#' @export
setGeneric("splitBy",
           signature = "x",
           function(x, ...)
               standardGeneric("splitBy"))

# This function collects f (grouping variable), MARGIN, skip_aggglomerate, and 
# use_names and returns them as a list.
.norm_args_for_split_by <- function(x, f, skip_agglomerate = FALSE,
                                    use_names = TRUE, ...){
    # input check
    if(is.null(f)){
        stop("'f' must either be a single non-empty character value or",
             " vector coercible to factor alongside one of the ",
             "dimensions of 'x'",
             call. = FALSE)
    }
    # Check f or extract the factor from rowData or colData
    if( !.is_non_empty_string(f) ){
        if( length(f) > 1L ){
            f <- factor(f, unique(f))
        }
        # Check if the length of f matches with one of the dimensions
        if(!length(f) %in% dim(x)){
            stop("'f' must either be a single non-empty character value or",
                 " vector coercible to factor alongside one of the ",
                 "dimensions of 'x'",
                 call. = FALSE)
        }
        # Get the dimension that matches
        MARGIN <- which(dim(x) %in% length(f))
        # If it matches with both directions
        if(length(MARGIN) > 1 ){
            # Get 2 if it can be found from colData, otherwise get 1
            MARGIN <- ifelse(any( sapply(colData(x), 
                                         function(var){all.equal(as.character(var), 
                                                                 as.character(f))}) == TRUE ), 
                             2, 1)
        }
    } else {
        # Try to get informaton from rowData
        tmp <- try({retrieveFeatureInfo(x, f, search = "rowData")},
                 silent = TRUE)
        if(is(tmp,"try-error")){
            # If it cannot be found from rowData, try to find it from colData
            tmp <- try({retrieveCellInfo(x, f, search = "colData")},
                     silent = TRUE)
            # Give error if it cannot be found from neither
            if(is(f,"try-error")){
                stop("", call. = FALSE)
            } 
            # Margin is columns
            MARGIN <- 2L
        } else {
            # If info can be found from rowData, margin is rows
            MARGIN <- 1L
        }
        # Get values and convert them into factors
        f <- tmp$value
        f <- factor(f, unique(f))
    }
    # Check skip_agglomerate
    if( !.is_a_bool(skip_agglomerate) ){
        stop("'skip_agglomerate must be TRUE or FALSE.",
             call. = FALSE)
    }
    # Check use_names
    if( !.is_a_bool(use_names) ){
        stop("'use_names must be TRUE or FALSE.",
             call. = FALSE)
    }
    # Create a list from arguments
    list(f = f,
         MARGIN = MARGIN,
         skip_agglomerate = skip_agglomerate,
         use_names = use_names)
}

.split_by <- function(x, args, ...){
    # Get grouping variable and its values
    f <- args[["f"]]
    # If user want to agglomerate data, agglomerate it based on MARGIN and 
    # and grouping variable
    if(!args[["skip_agglomerate"]]){
        merge_FUN <- switch(args[["MARGIN"]],
                            "1" = mergeRows,
                            "2" = mergeCols )
        return(merge_FUN(x, f, ...))
    }
    # Choose nrow or ncol based on MARGIN
    dim_FUN <- switch(args[["MARGIN"]],
                      "1" = nrow,
                      "2" = ncol)
    # Get indices from 1 to nrow/ncol
    idx <- seq_len(dim_FUN(x))
    # Split indices into groups based on grouping variable
    idxs <- split(idx, f)
    # SUbset function takes SE and list of groups which have indices
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

#' @rdname splitBy
#' @export
setMethod("splitBy", signature = c(x = "SummarizedExperiment"),
    function(x, f = NULL, skip_agglomerate = FALSE, ...){
        # Get arguments
        args <- .norm_args_for_split_by(x, f = f, 
                                        skip_agglomerate = skip_agglomerate)
        # Split data
        .split_by(x, args, ...)
    }
)

#' @rdname splitBy
#' @export
setMethod("splitBy", signature = c(x = "SingleCellExperiment"),
    function(x, f = NULL, skip_agglomerate = FALSE, ...){
        # Get arguments
        args <- .norm_args_for_split_by(x, f = f,
                                      skip_agglomerate = skip_agglomerate)
        # Should alternative experiment be removed? --> yes
        args[["strip_altexp"]] <- TRUE
        # Split data
        .split_by(x, args, ...)
    }
)

#' @rdname splitBy
#' @export
setMethod("splitBy", signature = c(x = "TreeSummarizedExperiment"),
    function(x, f = NULL, skip_agglomerate = FALSE, ...){
        callNextMethod()
    }
)

################################################################################
# unsplitBy

#' @rdname splitBy
#' @export
setGeneric("unsplitBy",
           signature = c("x"),
           function(x, ...)
               standardGeneric("unsplitBy"))

.list_unsplit_by <- function(ses, ...){
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
    rowData(ans) <- rd
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
    function(x, ...){
        # Unsplit list and create SCE, SE, or TreeSE from it
        .list_unsplit_by(x, ...)
    }
)
#' @rdname splitBy
#' @importFrom SingleCellExperiment altExpNames altExp altExps
#' @export
setMethod("unsplitBy", signature = c(x = "SimpleList"),
    function(x, ...){
        unsplitBy(as.list(x), ...)
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
        # Get alternative experiment names, since data is located there
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
   
#' @rdname unsplitBy
#' @export
setMethod("unsplitBy", signature = c(x = "TreeSummarizedExperiment"),
    function(x, altExpNames = altExpNames(x), keep_reducedDims = FALSE, ...){
        callNextMethod()
    }
)
