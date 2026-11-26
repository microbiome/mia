#' These functions will be deprecated. Please use other functions instead.
#' 
#' @param x A
#'   \code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#'   object.
#' 
#' @param ... Additional parameters. See dedicated function.
#' 
#' @seealso \code{\link{mia-defunct}} for the functions that are defunct.
#' 
#' @name deprecate
NULL

#' @rdname deprecate
#' @export
calculateDMN <- function(x, ...){
    if( is(x, "SummarizedExperiment") ){
        .Defunct(msg = paste0("'calculateDMN' is defunct for ",
            "SummarizedExperiment. Use 'addCluster' with ",
            "'bluster::DmmParam()' instead."))
    }
    .Deprecated(msg = paste0("'calculateDMN' is deprecated. Use ",
        "'bluster::clusterRows' with 'bluster::DmmParam()' instead."))
    .calculate_DMN(x, ...)
}

#' @rdname deprecate
#' @export
setGeneric("agglomerateByModule",
            signature = "x",
            function(x, ...)
                standardGeneric("agglomerateByModule"))

#' @rdname deprecate
#' @export
setMethod("agglomerateByModule", signature = c(x = "SummarizedExperiment"),
            function(x, ...){
                .Deprecated(msg = paste0("'agglomerateByModule' is deprecated.",
                                        " Use 'agglomerateBySet' instead."))
                agglomerateBySet(x, ...)
            }
)

#' @rdname deprecate
#' @export
setGeneric("subsetByPrevalent", signature = "x", function(x, ...)
    standardGeneric("subsetByPrevalent"))

#' @rdname deprecate
#' @export
setMethod("subsetByPrevalent", signature = c(x = "SummarizedExperiment"),
    function(x, ...){
        .Deprecated(msg = paste0(
            "'subsetByPrevalent' is deprecated. Use 'subsetBy' instead, ",
            "e.g., subsetBy(x, rows = rowMeans(counts > 0) > 0.5). To ",
            "agglomerate by rank first, call agglomerateByRank() before ",
            "subsetBy()."))
        .subset_by_prevalence(x, ...)
    }
)

#' @rdname deprecate
#' @export
setGeneric("agglomerateByPrevalence", signature = "x", function(x, ...)
    standardGeneric("agglomerateByPrevalence"))

#' @rdname deprecate
#' @export
setMethod("agglomerateByPrevalence", signature = c(x = "SummarizedExperiment"),
    function(x, ...){
        .Deprecated(msg = paste0(
            "'agglomerateByPrevalence' is deprecated. Use 'subsetBy' with ",
            "collapse = TRUE instead, e.g., x |> agglomerateByRank(rank) |> ",
            "subsetBy(rows = rowMeans(counts > 0) > 0.5, collapse = TRUE)."))
        .subset_by_prevalence(x, collapse = TRUE, ...)
    }
)

#' @rdname deprecate
#' @export
setGeneric("subsetByRare", signature = "x", function(x, ...)
    standardGeneric("subsetByRare"))

#' @rdname deprecate
#' @export
setMethod("subsetByRare", signature = c(x = "SummarizedExperiment"),
    function(x, ...){
        .Deprecated(msg = paste0(
            "'subsetByRare' is deprecated. Use 'subsetBy' instead, ",
            "e.g., subsetBy(x, rows = rowMeans(counts > 0) <= 0.5). To ",
            "agglomerate by rank first, call agglomerateByRank() before ",
            "subsetBy()."))
        .subset_by_prevalence(x, rare = TRUE, ...)
    }
)

# Shared implementation of the deprecated subsetByPrevalent(), subsetByRare(),
# agglomerateByPrevalence() and their aliases. The data is agglomerated by
# 'rank', prevalence is calculated with getPrevalence(), and the data is
# subsetted with subsetBy(). With collapse = TRUE, the features that are not
# prevalent are merged into one feature named 'other.name'.
.subset_by_prevalence <- function(
        x, rank = NULL, prevalence = 50/100, include.lowest = FALSE,
        rare = FALSE, collapse = FALSE, other.name = other_label,
        other_label = "Other", ...){
    if( !.is_numeric_string(prevalence) ){
        stop("'prevalence' must be a single numeric value or coercible to ",
            "one.", call. = FALSE)
    }
    prevalence <- as.numeric(prevalence)
    if( !.is_a_bool(include.lowest) ){
        stop("'include.lowest' must be TRUE or FALSE.", call. = FALSE)
    }
    if( collapse && !.is_a_string(other.name) ){
        stop("'other.name' must be a single character value.", call. = FALSE)
    }
    x <- .merge_features(x, rank = rank, ...)
    # Prevalence must be in the same order as the rows
    args <- list(...)
    args <- args[ !names(args) %in% c("sort", "rank") ]
    pr <- do.call(getPrevalence, c(list(x = x, sort = FALSE), args))
    keep <- if( include.lowest ) pr >= prevalence else pr > prevalence
    keep <- keep %in% TRUE
    if( rare ){
        keep <- !keep
    }
    # The condition is passed as a value with do.call(), so that rowData
    # variables cannot mask it.
    if( !collapse ){
        return(do.call(subsetBy, list(x = x, rows = keep)))
    }
    # Arguments for merging the features into "Other"
    merge_args <- args[ names(args) %in% c("update.tree", "update.refseq") ]
    res <- do.call(subsetBy, c(
        list(x = x, rows = keep, collapse = TRUE, collapse.name = other.name),
        merge_args))
    # The merged feature does not have taxonomy, except for the rank that was
    # used to agglomerate the data.
    if( !all(keep) && ncol(rowData(res)) > 0L ){
        other <- nrow(res)
        rd <- rowData(res)
        rd[other, ] <- NA
        if( !is.null(rank) && rank %in% colnames(rd) ){
            rd[other, rank] <- other.name
        }
        rowData(res) <- rd
    }
    return(res)
}
