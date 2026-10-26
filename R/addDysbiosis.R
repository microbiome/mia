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
#' @param assay.type \code{Character scalar}. Specifies the name of assay
#' used in calculation. (Default: \code{"counts"})
#'
#' @param index \code{Character vector}. Specifies the dysbiosis scores to
#' calculate. (Default: \code{"median_distance"})
#'
#' @param name \code{Character vector}. Specifies the column names for the
#' scores in \code{colData(x)}. (Default: \code{index})
#'
#' @param method \code{Character scalar}. Specifies the dissimilarity measure
#' passed to \code{\link[=getDissimilarity]{getDissimilarity}}.
#' (Default: \code{"bray"})
#'
#' @param ... Additional arguments. The following optional arguments are
#'   supported:
#'   \itemize{
#'     \item \code{reference}: \code{NULL} or \code{Character scalar}. Specifies
#'     the group in \code{col.var} that is used as controls. If \code{NULL}, the
#'     first level of \code{factor(x[[col.var]])} is used. (Default:
#'     \code{NULL})
#'
#'     \item \code{group}: \code{NULL} or \code{Character scalar}. Specifies a
#'     variable from \code{colData(x)} (e.g., subject ID). Control samples that
#'     share the same value with the current sample are excluded from its
#'     comparison. (Default: \code{NULL})
#'
#'     \item \code{altexp}: \code{NULL} or \code{Character scalar}. If provided,
#'     use this \code{altExp} assay container of \code{x} for the computation.
#'     (Default: \code{NULL})
#'
#'     \item Further arguments are passed to
#'     \code{\link[=getDissimilarity]{getDissimilarity}}.
#'   }
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
        x, col.var,
        assay.type = "counts",
        index = "median_distance", name = index,
        method = "bray",
        ...
){
        ############################## Input check #############################
        # Support altExp hiddenly
        x <- .check_and_get_altExp(x, ...)
        .check_assay_present(assay.type, x)
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
        ############################ Input check end ###########################
        # Factor defining the group membership of each sample; reference group
        # is placed as the first level.
        reference_groups <- .get_reference_for_dysbiosis(x, col.var, ...)

        # Optional grouping variable (e.g., patient ID) used to exclude
        # distances between samples from the same subject.
        sample_groups <- .get_grouping_for_dysbiosis(x, col.var, ...)

        mat <- assay(x, assay.type)

        res <- lapply(index, function(i){
            FUN <- switch(
                i,
                median_distance = .calc_median_distance_to_controls
            )
            FUN(
                mat = mat,
                groups = reference_groups,
                sample_groups = sample_groups,
                method = method,
                ...
            )
        })
        names(res) <- name
        res <- DataFrame(res, row.names = colnames(x), check.names = FALSE)

        return(res)
    }
)

# Determine reference group ordering for dysbiosis calculation
#
# Creates a factor from colData(x)[[col.var]] and ensures that the chosen
# reference level appears first.
.get_reference_for_dysbiosis <- function(x, col.var, reference = NULL, ...){
    if( !(.is_a_string(col.var) && col.var %in% colnames(colData(x))) ){
        stop("'col.var' must specify a column from colData(x).", call. = FALSE)
    }

    groups <- factor(x[[col.var]])
    if( nlevels(groups) < 2L ){
        stop("'col.var' must define at least two groups.", call. = FALSE)
    }

    lvls <- levels(groups)
    if( is.null(reference) ){
        reference <- lvls[[1L]]
    }
    if( !(.is_a_string(reference) && reference %in% lvls) ){
        stop("'reference' must be NULL or specify a group in 'col.var'.",
             call. = FALSE)
    }

    # Ensure reference is first level
    groups <- factor(groups, levels = c(reference, lvls) |> unique())

    return(groups)
}

# Optional grouping variable for excluding within-subject distances
.get_grouping_for_dysbiosis <- function(x, col.var, group = NULL, ...){
    if( is.null(group) ){
        return(NULL)
    }
    if( !(.is_a_string(group) &&
            group %in% colnames(colData(x)) && group != col.var) ){
        stop("'group' must specify a different column from colData(x) than ",
            "'col.var'.",
            call. = FALSE)
    }
    group <- x[[group]]
    return(group)
}

# Median distance to reference (control) samples
#
# For each sample, compute the median dissimilarity to samples in the reference
# group. Self-distances are excluded. If sample_groups is provided, distances
# between samples with the same value of sample_groups are excluded (e.g.,
# repeated measures from the same patient).
# @importFrom DelayedMatrixStats rowMedians
.calc_median_distance_to_controls <- function(
        mat, groups, sample_groups, method, altexp = NULL, ...){

    # Calculate dissimilarities
    diss <- getDissimilarity(mat |> t(), method = method, ...)
    diss <- diss |> as.matrix()
    diag(diss) <- NA_real_

    # Optionally exclude distances between samples from the same subject/group
    if( !is.null(sample_groups) ){
        same <- outer(sample_groups, sample_groups, "==")
        same[is.na(same)] <- FALSE
        diss[same] <- NA_real_
    }

    # Subset to reference samples (reference group is level 1)
    is_control <- as.integer(groups) == 1L
    diss <- diss[, is_control, drop = FALSE]

    res <- rowMedians(diss, na.rm = TRUE)
    res <- res |> unname()

    return(res)
}
