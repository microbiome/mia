#' Split/Unsplit a \code{SingleCellExperiment} by taxonomic ranks
#'
#' \code{splitByRanks} takes a \code{SummarizedExperiment}, splits it along the
#' taxonomic ranks, aggregates the data per rank, converts the input to a 
#' \code{SingleCellExperiment} objects and stores the aggregated data as 
#' alternative experiments.
#' 
#' \code{unsplitByRanks} takes these alternative experiments and flattens them 
#' again into a single \code{SummarizedExperiment}.
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
#' @param ... arguments passed to \code{agglomerateByRank} function for
#'   \code{SummarizedExperiment} objects and other functions.
#'   See \code{\link[=agglomerate-methods]{agglomerateByRank}} for more details.
#'
#' @return
#' For \code{splitByRanks}: \code{x}, with objects of \code{x} agglomerated for
#' selected ranks as \code{altExps}.
#'
#' For \code{unsplitByRanks}: \code{x}, with \code{rowData} and \code{assay}
#' data replaced by the unsplit data. \code{colData} of x is kept as well
#' and any existing \code{rowTree} is dropped as well, since existing
#' \code{rowLinks} are not valid anymore.
#'
#' @details
#' \code{splitByRanks} will use by default all available taxonomic ranks, but
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
#' @seealso
#' \code{\link[=merge-methods]{mergeRows}},
#' \code{\link[scuttle:sumCountsAcrossFeatures]{sumCountsAcrossFeatures}},
#' \code{\link[=agglomerate-methods]{agglomerateByRank}},
#' \code{\link[SingleCellExperiment:altExps]{altExps}},
#' \code{\link[SingleCellExperiment:splitAltExps]{splitAltExps}}
#'
#' @name splitByRanks
#'
#' @examples
#' data(GlobalPatterns)
#' # print the available taxonomic ranks
#' taxonomyRanks(GlobalPatterns)
#'
#' # splitByRanks
#' altExps(GlobalPatterns) <- splitByRanks(GlobalPatterns)
#' altExps(GlobalPatterns)
#' altExp(GlobalPatterns,"Kingdom")
#' altExp(GlobalPatterns,"Species")
#'
#' # unsplitByRanks
#' x <- unsplitByRanks(GlobalPatterns)
#' x
NULL

################################################################################
# splitByRanks

#' @rdname splitByRanks
#' @export
setGeneric("splitByRanks",
           signature = "x",
           function(x, ...)
               standardGeneric("splitByRanks"))

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

#' @rdname splitByRanks
#' @export
setMethod("splitByRanks", signature = c(x = "SummarizedExperiment"),
    function(x, ranks = taxonomyRanks(x), na.rm = TRUE, ...){
        args <- .norm_args_for_split_by_ranks(na.rm = na.rm, ...)
        .split_by_ranks(x, ranks, args)
    }
)

#' @rdname splitByRanks
#' @export
setMethod("splitByRanks", signature = c(x = "SingleCellExperiment"),
          function(x, ranks = taxonomyRanks(x), na.rm = TRUE, ...){
              args <- .norm_args_for_split_by_ranks(na.rm = na.rm, ...)
              args[["strip_altexp"]] <- TRUE
              .split_by_ranks(x, ranks, args)
          }
)

#' @rdname splitByRanks
#' @export
setMethod("splitByRanks", signature = c(x = "TreeSummarizedExperiment"),
          function(x, ranks = taxonomyRanks(x), na.rm = TRUE, ...){
              callNextMethod()
          }
)


################################################################################
# splitByAlternativeRanks

#' @rdname splitByRanks
#' @export
setGeneric("splitByAlternativeRanks",
           signature = "x",
           function(x, ...)
               standardGeneric("splitByAlternativeRanks"))

.merge_rows_wrapper <- function(x, col, ...){
    f <- rowData(x)[,col,drop=TRUE]
    f <- factor(f, unqiue(f))
    mergeRows(x, f = f,) 
}

.split_by_alt_ranks <- function(x, alt_ranks, ...){
    # input check
    if(!.is_non_empty_character(alt_ranks)){
        stop("'alt_ranks' must be character vector.",
             call. = FALSE)
    }
    if(nrow(x) == 0L){
        stop("'x' has nrow(x) == 0L.",call. = FALSE)
    }
    #
    ans <- lapply(alt_ranks, .merge_rows_wrapper, x = x)
    names(ans) <- alt_ranks
    SimpleList(ans)
}

#' @rdname splitByRanks
#' @export
setMethod("splitByAlternativeRanks", signature = c(x = "SummarizedExperiment"),
          function(x, alt_ranks, ...){
              .split_by_alt_ranks(x, ranks, ...)
          }
)

#' @rdname splitByRanks
#' @export
setMethod("splitByAlternativeRanks", signature = c(x = "SingleCellExperiment"),
          function(x, alt_ranks, ...){
              altExps(x) <- NULL
              x <- .split_by_alt_ranks(x, ranks, ...)
              x
          }
)

#' @rdname splitByRanks
#' @export
setMethod("splitByAlternativeRanks", signature = c(x = "TreeSummarizedExperiment"),
          function(x, alt_ranks, ...){
              callNextMethod()
          }
)

################################################################################
# unsplitByRanks

#' @rdname splitByRanks
#' @export
setGeneric("unsplitByRanks",
           signature = "x",
           function(x, ...)
               standardGeneric("unsplitByRanks"))


#' @rdname splitByRanks
#' @importFrom SingleCellExperiment altExpNames altExp altExps reducedDims
#' @export
setMethod("unsplitByRanks", signature = c(x = "SingleCellExperiment"),
    function(x, ranks = taxonomyRanks(x), keep_reducedDims = FALSE, ...){
        # input check
        if(!.is_a_bool(keep_reducedDims)){
          stop("'keep_reducedDims' must be TRUE or FALSE.", call. = FALSE)
        }
        #
        class_x <- class(x)
        ae_names <- altExpNames(x)
        ae_names <- ae_names[ae_names %in% ranks]
        if(length(ae_names) == 0L){
          stop("No altExp matching 'ranks' in name.", call. = FALSE)
        }
        ses <- altExps(x)[ae_names]
        # remove any empty information on the given ranke
        for(i in seq_along(ses)){
          ses[[i]] <-
              .remove_with_empty_taxonomic_info(ses[[i]], names(ses)[i], NA)
        }
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
        rownames(ans) <- getTaxonomyLabels(ans, make_unique = FALSE)
        ans
    }
)

#' @rdname splitByRanks
#' @importFrom SingleCellExperiment altExpNames altExp altExps reducedDims
#' @export
setMethod("unsplitByRanks", signature = c(x = "TreeSummarizedExperiment"),
    function(x, ranks = taxonomyRanks(x), keep_reducedDims = FALSE, ...){
        callNextMethod()
    }
)


#' @importFrom SummarizedExperiment assayNames assay
#' @importFrom DelayedArray DelayedArray
.get_assay <- function(se, assay_name){
    if (assay_name %in% assayNames(se)) {
        assay <- DelayedArray(assay(se, assay_name))
    } else {
        dim <- dim(se)
        assay <- DelayedArray(matrix(NA_real_,dim[[1L]],dim[[2L]]))
    }
    assay
}

#' @importFrom BiocGenerics rbind
.combine_assays <- function(ses, assay_name){
    current <- lapply(ses, .get_assay, assay_name)
    combined <- do.call(rbind, current)
    rownames(combined) <- NULL
    colnames(combined) <- NULL
    combined
}

#' @importFrom SummarizedExperiment assayNames assay
.unsplit_assays <- function(ses) {
    assay_names <- unique(unlist(lapply(ses, assayNames)))
    combined <- lapply(assay_names, .combine_assays,
                       ses = ses)
    names(combined) <- assay_names
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
