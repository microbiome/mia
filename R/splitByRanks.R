#' Agglomerate a \code{SummarizedExperiment} based on several taxonomic ranks
#'
#' \code{agglomerateByRanks} takes a \code{SummarizedExperiment}, splits it along the
#' taxonomic ranks, aggregates the data per rank, converts the input to a 
#' \code{SingleCellExperiment} objects and stores the aggregated data as 
#' alternative experiments. \code{unsplitByRanks} takes these alternative 
#' experiments and flattens them again into a single 
#' \code{SummarizedExperiment}.
#'
#' @param x a
#'   \code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#'   object
#'
#' @param ranks a character vector defining taxonomic ranks. Must all be values
#'   of \code{taxonomyRanks()} function.
#'
#' @param na.rm \code{TRUE} or \code{FALSE}: Should taxa with an empty rank be
#'   removed? Use it with caution, since results with NA on the selected rank
#'   will be dropped. This setting can be tweaked by defining
#'   \code{empty.fields} to your needs. (default: \code{na.rm = TRUE})
#'
#' @param keep_reducedDims \code{TRUE} or \code{FALSE}: Should the
#'   \code{reducedDims(x)} be transferred to the result? Please note, that this
#'   breaks the link between the data used to calculate the reduced dims.
#'   (default: \code{keep_reducedDims = FALSE})
#'   
#' @param as.list \code{TRUE} or \code{FALSE}: Should the list of 
#'   \code{SummarizedExperiment} objects be returned by the function 
#'   \code{agglomerateByRanks} as a SimpleList or stored in altExps?
#'   (default: \code{as.list = FALSE})
#'
#' @param ... arguments passed to \code{agglomerateByRank} function for
#'   \code{SummarizedExperiment} objects and other functions.
#'   See \code{\link[=agglomerate-methods]{agglomerateByRank}} for more details.
#'
#' @return
#' For \code{agglomerateByRanks}: 
#' If \code{as.list = TRUE} : \code{SummarizedExperiment} objects in a 
#' \code{SimpleList} 
#' If \code{as.list = FALSE} : The \code{SummarizedExperiment} passed as a 
#' parameter and now containing the \code{SummarizedExperiment} objects in its
#' altExps
#'
#' For \code{unsplitByRanks}: \code{x}, with \code{rowData} and \code{assay}
#' data replaced by the unsplit data. \code{colData} of x is kept as well
#' and any existing \code{rowTree} is dropped as well, since existing
#' \code{rowLinks} are not valid anymore.
#'
#' @details
#' \code{agglomerateByRanks} will use by default all available taxonomic ranks, but
#' this can be controlled by setting \code{ranks} manually. \code{NA} values
#' are removed by default, since they would not make sense, if the result
#' should be used for \code{unsplitByRanks} at some point. The input data 
#' remains unchanged in the returned \code{SingleCellExperiment} objects.
#'
#' \code{unsplitByRanks} will remove any \code{NA} value on each taxonomic rank
#' so that no ambiguous data is created. In additional, a column
#' \code{taxonomicLevel} is created or overwritten in the \code{rowData} to
#' specify from which alternative experiment this originates from. This can also
#' be used for \code{\link[SingleCellExperiment:splitAltExps]{splitAltExps}} to
#' split the result along the same factor again. The input data from the base
#' objects is not returned, only the data from the \code{altExp()}. Be aware that
#' changes to \code{rowData} of the base object are not returned, whereas only 
#' the \code{colData} of the base object is kept. 
#'
#' @rdname agglomerate-methods
#'
#' @examples
#' data(GlobalPatterns)
#' # print the available taxonomic ranks
#' taxonomyRanks(GlobalPatterns)
#'
#' # agglomerateByRanks
#' # 
#' tse <- agglomerateByRanks(GlobalPatterns)
#' altExps(tse)
#' altExp(tse,"Kingdom")
#' altExp(tse,"Species")
#'
#' # unsplitByRanks
#' tse <- unsplitByRanks(tse)
#' tse
#'
#' @aliases splitByRanks
#' @export
setGeneric("agglomerateByRanks",
           signature = "x",
           function(x, ...)
               standardGeneric("agglomerateByRanks"))

.norm_args_for_split_by_ranks <- function(na.rm, ...){
    args <- list(...)
    if(missing(na.rm)){
        na.rm <- TRUE
    }
    args[["na.rm"]] <- na.rm
    args
}

.split_by_ranks <- function(x, ranks, args){
    # input check
    if(!.is_non_empty_character(ranks)){
        stop("'ranks' must be character vector.",
             call. = FALSE)
    }
    if(nrow(x) == 0L){
        stop("'x' has nrow(x) == 0L.",call. = FALSE)
    }
    .check_taxonomic_ranks(ranks,x)
    #
    FUN <- function(rank){
        do.call(agglomerateByRank,
                c(list(x = x, rank = rank), args))
    }
    ans <- lapply(ranks,FUN)
    names(ans) <- ranks
    SimpleList(ans)
}

#' @rdname agglomerate-methods
#' @aliases splitByRanks
#' @export
setMethod("agglomerateByRanks", signature = c(x = "SummarizedExperiment"),
    function(x, ranks = taxonomyRanks(x), na.rm = TRUE, as.list = FALSE, ...){
        if ( is(x, "SummarizedExperiment") ){
            args <- .norm_args_for_split_by_ranks(na.rm = na.rm, ...)
            x <- as(x, "TreeSummarizedExperiment")
            warning("SummarizedExperiment does not have altExps slot. ",
                    "Therefore, it has been converted to ",
                    "TreeSummarizedExperiment.")
        }
        res <- .split_by_ranks()
        # Add to altExp if user has specified to do so
        if( !as.list ){
            # Create a general function that adds results to altExps (That could be used in other functions also in the future)
            res <- .add_to_altExps(x, res)
        }
        return(res)
    }
)

#' @rdname agglomerate-methods
#' @aliases splitByRanks
#' @export
setMethod("agglomerateByRanks", signature = c(x = "SingleCellExperiment"),
            function(x, ranks = taxonomyRanks(x), na.rm = TRUE, as.list = FALSE,
                    ...){
                args <- .norm_args_for_split_by_ranks(na.rm = na.rm, ...)
                args[["strip_altexp"]] <- TRUE
                callNextMethod()
            }
)

#' @rdname agglomerate-methods
#' @aliases splitByRanks
#' @export
setMethod("agglomerateByRanks", signature = c(x = "TreeSummarizedExperiment"),
          function(x, ranks = taxonomyRanks(x), na.rm = TRUE, as.list = FALSE,
                   ...){
              callNextMethod()
          }
)

################################################################################
# splitByRanks

#' @rdname agglomerate-methods
#' @export 
setGeneric("splitByRanks",
           signature = "x",
           function(x, ...)
               standardGeneric("splitByRanks"))

#' @rdname agglomerate-methods
#' @export 
setMethod("splitByRanks", signature = c(x = "SummarizedExperiment"),
          function(x,...){
              .Deprecated(msg = paste0("'splitByRanks' is deprecated. ",
                                       "Use 'agglomerateByRanks' instead."))
              agglomerateByRanks(x, as.list = FALSE,...)
          }
)

################################################################################
# unsplitByRanks

#' @rdname agglomerate-methods
#' @export
setGeneric("unsplitByRanks",
           signature = "x",
           function(x, ...)
               standardGeneric("unsplitByRanks"))


#' @importFrom SingleCellExperiment reducedDims
#' @importFrom SummarizedExperiment colData
.unsplit_by <- function(x, ses, keep_reducedDims, ...){
    class_x <- class(x)
    #
    args <- list(assays = .unsplit_assays(ses),
                 colData = colData(x))
    if(keep_reducedDims){
        args$reducedDims <- reducedDims(x)
    }
    rd <- .combine_rowData(ses)
    rr <- .combine_rowRanges(ses)
    args$rowRanges <- rr
    ans <- do.call(class_x, args)
    rowData(ans) <- rd
    ans
}

#' @importFrom SingleCellExperiment altExpNames altExp altExps
.unsplit_by_ranks <- function(x, ranks, keep_reducedDims, ...){
    ae_names <- altExpNames(x)
    ae_names <- ae_names[ae_names %in% ranks]
    if(length(ae_names) == 0L){
        stop("No altExp matching 'ranks' in name.", call. = FALSE)
    }
    ses <- altExps(x)[ae_names]
    # remove any empty information on the given ranks
    for(i in seq_along(ses)){
        ses[[i]] <-
            .remove_with_empty_taxonomic_info(ses[[i]], names(ses)[i], NA)
    }
    ans <- .unsplit_by(x, ses, keep_reducedDims, ...)
    rownames(ans) <- getTaxonomyLabels(ans, make_unique = FALSE)
    ans
}

#' @rdname agglomerate-methods
#' @export
setMethod("unsplitByRanks", signature = c(x = "SingleCellExperiment"),
    function(x, ranks = taxonomyRanks(x), keep_reducedDims = FALSE, ...){
        # input check
        if(!.is_a_bool(keep_reducedDims)){
          stop("'keep_reducedDims' must be TRUE or FALSE.", call. = FALSE)
        }
        #
        .unsplit_by_ranks(x, ranks = ranks, keep_reducedDims = keep_reducedDims,
                          ...)
    }
)

#' @rdname agglomerate-methods
#' @export
setMethod("unsplitByRanks", signature = c(x = "TreeSummarizedExperiment"),
    function(x, ranks = taxonomyRanks(x), keep_reducedDims = FALSE, ...){
        callNextMethod()
    }
)


#' @importFrom SummarizedExperiment assayNames assay
#' @importFrom DelayedArray DelayedArray
.get_assay <- function(se, assay.type){
    if (assay.type %in% assayNames(se)) {
        assay <- DelayedArray(assay(se, assay.type))
    } else {
        dim <- dim(se)
        assay <- DelayedArray(matrix(NA_real_,dim[[1L]],dim[[2L]]))
    }
    assay
}

#' @importFrom BiocGenerics rbind cbind
.combine_assays <- function(ses, assay.type, MARGIN = 1L){
    bind_FUN <- switch(MARGIN, "1" = rbind, "2" = cbind)
    current <- lapply(ses, .get_assay, assay.type)
    combined <- do.call(bind_FUN, current)
    rownames(combined) <- NULL
    colnames(combined) <- NULL
    combined
}

#' @importFrom SummarizedExperiment assayNames assay
.unsplit_assays <- function(ses, MARGIN = 1L) {
    assay.types <- unique(unlist(lapply(ses, assayNames)))
    combined <- lapply(assay.types,
                       .combine_assays,
                       ses = ses,
                       MARGIN = MARGIN)
    names(combined) <- assay.types
    combined
}

#' @importFrom SummarizedExperiment rowRanges
.combine_rowRanges <- function(ses, prefix) {
    rr <- lapply(ses, rowRanges)
    unname(do.call(c,unname(rr)))
}


#' @importFrom SummarizedExperiment rowData
#' @importFrom BiocGenerics rbind
.combine_rowData <- function(ses, prefix) {
    rr <- lapply(ses, rowData)
    rd <- do.call(rbind,unname(rr))
    rownames(rd) <- NULL
    cn <- colnames(rd)
    # add column
    if("taxonomicLevel" %in% cn){
        warning("'taxonomicLevel' column in rowData overwritten.",
                call. = FALSE)
    }
    tl <- mapply(rep,
                 names(ses),
                 vapply(ses,nrow,integer(1)))
    tl <- unlist(unname(tl))
    rd$taxonomicLevel <- factor(tl, unique(tl))
    rd
}
