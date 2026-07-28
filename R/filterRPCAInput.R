#' @name
#' filterRPCAInput
#'
#' @title
#' Filter raw count assays before RPCA or Joint-RPCA
#'
#' @description
#' \code{filterRPCAInput()} applies Gemelli-style preprocessing filters to raw
#' count assays before robust centered log-ratio transformation and RPCA.
#'
#' The filters are intended for raw count data, not transformed assays such as
#' rclr. Samples are filtered by total count. Features are filtered by total
#' count and prevalence across samples.
#'
#' @details
#' This function mirrors Gemelli-style RPCA table preprocessing:
#'
#' * sample total count must be `> min.sample.count`
#' * feature total count must be `> min.feature.count`
#' * feature prevalence percentage must be `> min.feature.frequency`
#'
#' The strict `>` comparison is intentional.
#'
#' For exact comparison with the current Gemelli script, use:
#'
#' ```
#' min.sample.count = 0
#' min.feature.count = 0
#' min.feature.frequency = 0
#' ```
#'
#' @return
#' Filtered \code{SummarizedExperiment} object.
#'
#' @inheritParams addAlpha
#'
#' @param min.sample.count \code{Numeric scalar} or \code{NULL}. Specifies the
#' minimum sample total count. Samples with total counts less than or equal to
#' this value are removed. If \code{NULL}, this filter is skipped.
#' (Default: \code{0})
#'
#' @param min.feature.count \code{Numeric scalar} or \code{NULL}. Specifies the
#' minimum feature total count. Features with total counts less than or equal
#' to this value are removed. If \code{NULL}, this filter is skipped.
#' (Default: \code{0})
#'
#' @param min.feature.frequency \code{Numeric scalar} or \code{NULL}. Specifies
#' the minimum feature prevalence across samples as a proportion between 0 and
#' 1. Features with prevalence less than or equal to this value are removed.
#' If \code{NULL}, this filter is skipped. (Default: \code{0})
#'
#' @param ... additional arguments.
#'
#' @examples
#'
#' data(GlobalPatterns)
#' tse <- GlobalPatterns
#'
#' tse_sub <- filterRPCAInput(
#'     tse,
#'     min.sample.count = 5,
#'     min.feature.count = 10,
#'     min.feature.frequency = 0.3
#' )
#'
#' @seealso
#' \code{\link[=addRPCA]{addRPCA}} and
#' \code{\link[=addJointRPCA]{addJointRPCA}}
#'
#' @references
#'
#' Martino, C. and Shenhav, L. et al. (2020)
#' Context-aware dimensionality reduction deconvolutes gut microbial community
#' dynamics.
#' _Nat. Biotechnol._ doi:10.1038/s41587-020-0660-7
#'
NULL

#' @rdname filterRPCAInput
#' @export
setMethod("filterRPCAInput",
    signature = c(x = "SummarizedExperiment"),
    function(
        x,
        assay.type = "counts",
        min.sample.count = 0,
        min.feature.count = 0,
        min.feature.frequency = 0,
        ...
    ) {
        .check_assay_present(assay.type, x)
        if (!(is.null(min.sample.count) ||
            (.is_an_integer(min.sample.count) && min.sample.count >= 0))) {
            stop("'min.sample.count' must be a single positive integer value.",
                call. = FALSE
            )
        }
        if (!(is.null(min.feature.count) ||
            (.is_an_integer(min.feature.count) && min.feature.count >= 0))) {
            stop("'min.feature.count' must be a single positive integer value.",
                call. = FALSE
            )
        }
        if (!(is.null(min.feature.frequency) ||
            (.is_a_numeric(min.feature.frequency) &&
                min.feature.frequency >= 0 &&
                min.feature.frequency <= 1))) {
            stop("'min.feature.frequency' must be a single numeric value ",
                "between 0 and 1.",
                call. = FALSE
            )
        }
        #
        x_sub <- .filter_based_on_abundance(
            x,
            assay.type = assay.type,
            min.sample.count = min.sample.count,
            min.feature.count = min.feature.count,
            min.feature.frequency = min.feature.frequency,
            ...
        )
        return(x_sub)
    }
)

############################### HELPER FUNCTIONS ###############################

.filter_based_on_abundance <- function(
    tse,
    assay.type = "counts",
    min.sample.count = 0,
    min.feature.count = 0,
    min.feature.frequency = 0,
    na.rm = FALSE,
    ...
) {
    if (!.is_a_bool(na.rm)) {
        stop("'na.rm' must be TRUE or FALSE.", call. = FALSE)
    }
    # Calculate stats
    mat <- assay(tse, assay.type)
    col_sums <- mat |> colSums(na.rm = na.rm)
    row_sums <- mat |> rowSums(na.rm = na.rm)
    row_frequency <- rowSums(mat > 0, na.rm = na.rm) / ncol(mat)

    # Get samples and features that exceed the thresholds
    col_index <- rep(TRUE, ncol(mat))
    row_index <- rep(TRUE, nrow(mat))

    if (!is.null(min.sample.count)) {
        col_index <- col_index & col_sums > min.sample.count
    }
    if (!is.null(min.feature.count)) {
        row_index <- row_index & row_sums > min.feature.count
    }
    if (!is.null(min.feature.frequency)) {
        row_index <- row_index & row_frequency > min.feature.frequency
    }

    # Do filtering
    tse <- tse[row_index, col_index]

    return(tse)
}
