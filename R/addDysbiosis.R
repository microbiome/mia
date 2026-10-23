#' Calculate dysbiosis scores
#'
#' Dysbiosis scores measure how much the community of a sample differs from
#' the communities of control samples.
#'
#' @param x a \code{\link[SummarizedExperiment]{SummarizedExperiment}} object.
#'
#' @param col.var \code{Character scalar}. Specifies a variable from
#' \code{colData(x)} that defines the groups.
#'
#' @param reference \code{NULL} or \code{Character scalar}. Specifies the group
#' in \code{col.var} that is used as controls. If \code{NULL}, the first level
#' of \code{factor(x[[col.var]])} is used. (Default: \code{NULL})
#'
#' @param assay.type \code{Character scalar}. Specifies the name of assay
#' used in calculation. (Default: \code{"counts"})
#'
#' @param index \code{Character vector}. Specifies the dysbiosis scores to
#' calculate. (Default: \code{"median_distance"})
#'
#' @param name \code{Character vector}. Specifies the column names for the
#' scores in \code{colData(x)}. (Default: \code{index})
#'
#' @param group \code{NULL} or \code{Character scalar}. Specifies a variable
#' from \code{colData(x)}, e.g., subject. Control samples that share the value
#' with the sample are excluded from its comparison. (Default: \code{NULL})
#'
#' @param method \code{Character scalar}. Specifies the dissimilarity measure
#' passed to \code{\link[=getDissimilarity]{getDissimilarity}}.
#' (Default: \code{"bray"})
#'
#' @param ... additional arguments passed to
#' \code{\link[=getDissimilarity]{getDissimilarity}}.
#'
#' @details
#' The following scores are available:
#' \itemize{
#'   \item \code{median_distance}: Median dissimilarity between the sample and
#'   the control samples. The sample itself is excluded from the controls. With
#'   Bray-Curtis dissimilarity, this is the dysbiosis score of Lloyd-Price et
#'   al. (2019).
#' }
#'
#' @return
#' \code{getDysbiosis} returns a \code{DataFrame} with one row per sample.
#' \code{addDysbiosis} returns \code{x} with the scores added to
#' \code{colData(x)}.
#'
#' @references
#' Lloyd-Price J, et al. (2019) Multi-omics of the gut microbial ecosystem in
#' inflammatory bowel diseases. \emph{Nature} 569, 655-662.
#'
#' @seealso
#' \code{\link[=getDissimilarity]{getDissimilarity}},
#' \code{\link[=addDivergence]{addDivergence}}
#'
#' @name addDysbiosis
#'
#' @examples
#' data(Tengeler2020)
#' tse <- Tengeler2020
#'
#' tse <- addDysbiosis(tse, col.var = "patient_status", reference = "Control")
#' colData(tse)[, c("patient_status", "median_distance")]
#'
NULL

#' @rdname addDysbiosis
#' @export
setMethod("addDysbiosis", signature = c(x = "SummarizedExperiment"),
    function(x, ...){
        res <- getDysbiosis(x, ...)
        x <- .add_values_to_colData(
            x, unname(as.list(res)), colnames(res),
            altexp = list(...)[["altexp"]])
        return(x)
    }
)

#' @rdname addDysbiosis
#' @export
setMethod("getDysbiosis", signature = c(x = "SummarizedExperiment"),
    function(
        x, col.var, reference = NULL, assay.type = "counts",
        index = "median_distance", name = index, group = NULL,
        method = "bray", ...){
        ############################## Input check #############################
        # Support altExp hiddenly
        x <- .check_and_get_altExp(x, ...)
        .check_assay_present(assay.type, x)
        if( !(.is_a_string(col.var) && col.var %in% colnames(colData(x))) ){
            stop("'col.var' must specify a column from colData(x).",
                call. = FALSE)
        }
        groups <- factor(x[[col.var]])
        if( nlevels(groups) < 2L ){
            stop("'col.var' must define at least two groups.", call. = FALSE)
        }
        if( is.null(reference) ){
            reference <- levels(groups)[[1L]]
        }
        if( !(.is_a_string(reference) && reference %in% levels(groups)) ){
            stop("'reference' must be NULL or specify a group in 'col.var'.",
                call. = FALSE)
        }
        supported <- c("median_distance")
        if( !(.is_non_empty_character(index) && all(index %in% supported)) ){
            stop("'index' must be one of the following: '",
                paste0(supported, collapse = "', '"), "'", call. = FALSE)
        }
        if( !(.is_non_empty_character(name) &&
                length(name) == length(index) && !anyDuplicated(name)) ){
            stop("'name' must be unique and have the same length as 'index'.",
                call. = FALSE)
        }
        if( !(is.null(group) ||
                (.is_a_string(group) && group %in% colnames(colData(x)))) ){
            stop("'group' must be NULL or specify a column from colData(x).",
                call. = FALSE)
        }
        ############################ Input check end ###########################
        is_control <- !is.na(groups) & groups == reference
        res <- lapply(index, function(i){
            switch(i, median_distance = .calc_median_distance_to_controls(
                x, assay.type, is_control, group, method, ...))
        })
        names(res) <- name
        res <- DataFrame(res, row.names = colnames(x), check.names = FALSE)
        return(res)
    }
)

# For each sample, this function calculates the median dissimilarity to the
# control samples, excluding the sample itself and samples from the same group.
#' @importFrom stats median
.calc_median_distance_to_controls <- function(
        x, assay.type, is_control, group, method, altexp = NULL, ...){
    diss <- getDissimilarity(
        x, assay.type = assay.type, method = method, ...) |> as.matrix()
    diss <- diss[colnames(x), colnames(x), drop = FALSE]
    diag(diss) <- NA
    if( !is.null(group) ){
        same <- outer(x[[group]], x[[group]], "==")
        same[is.na(same)] <- FALSE
        diss[same] <- NA
    }
    diss <- diss[, is_control, drop = FALSE]
    res <- apply(diss, 1L, median, na.rm = TRUE)
    return(unname(res))
}
