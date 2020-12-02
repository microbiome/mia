#' Split a \code{SingleCellExperiment} be taxonomic ranks
#'
#' \code{splitByRanks}
#'
#' @param x \code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#'
#' @param ranks a character vector defining taxonomic ranks. Must all be values
#'   of \code{taxonomicRanks()} function.
#'
#' @param ... arguments passed to \code{agglomerateByRank} function for
#'   \code{SummarizedExperiment} objects and other functions.
#'   See \code{\link[=agglomerate-methods]{agglomerateByRank}} for more details.
#'
#' @return \code{x}, with objects of \code{x} agglomerated for selected ranks
#'   as \code{altExps}.
#'
#' @seealso
#' \code{\link[=merge-methods]{mergeRows}},
#' \code{\link[scuttle:sumCountsAcrossFeatures]{sumCountsAcrossFeatures}},
#' \code{\link[=agglomerate-methods]{agglomerateByRank}},
#' \code{\link[SingleCellExperiment:altExps]{altExps}}
#'
#' @name splitByRanks
#'
#' @examples
#' data(GlobalPatterns)
#' # print the available taxonomic ranks
#' taxonomyRanks(GlobalPatterns)
#'
#' altExps(GlobalPatterns) <- splitByRanks(GlobalPatterns)
#' altExps(GlobalPatterns)
#' altExp(GlobalPatterns,"Kingdom")
#' altExp(GlobalPatterns,"Species")
NULL

setGeneric("splitByRanks",
           signature = "x",
           function(x, ...)
               standardGeneric("splitByRanks"))

setGeneric("unsplitByRanks",
           signature = "x",
           function(x, ...)
               standardGeneric("unsplitByRanks"))

#' @rdname splitByRanks
#' @export
setMethod("splitByRanks", signature = c(x = "TreeSummarizedExperiment"),
    function(x, ranks = taxonomyRanks(x), na.rm = TRUE, ...){
        # input check
        if(!.is_non_empty_character(ranks)){
            stop("'ranks' must be character vector.",
                 call. = FALSE)
        }
        .check_taxonomic_ranks(ranks,x)
        #
        args <- list(...)
        if(missing(na.rm)){
            na.rm <- TRUE
        }
        args[["na.rm"]] <- na.rm
        args[["strip_altexp"]] <- TRUE
        FUN <- function(rank){
            do.call(agglomerateByRank,
                    c(list(x = x, rank = rank), args))
        }
        ans <- lapply(ranks,FUN)
        names(ans) <- ranks
        SimpleList(ans)
    }
)

#' @rdname splitByRanks
#' @export
setMethod("unsplitByRanks", signature = c(x = "TreeSummarizedExperiment"),
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

        args <- list(assays = .unsplit_assays(ses),
                     colData = colData(x),
                     rowTree = rowTree(x))
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
    rd
}
