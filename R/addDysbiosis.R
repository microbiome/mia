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
#' supported:
#' \itemize{
#'   \item \code{reference}: \code{NULL} or \code{Character scalar}. Specifies
#'   the group in \code{col.var} that is used as controls. If \code{NULL}, the
#'   first level of \code{factor(x[[col.var]])} is used. (Default:
#'   \code{NULL})
#'   \item \code{group}: \code{NULL} or \code{Character scalar}. Specifies a
#'   variable from \code{colData(x)} (e.g., subject ID). Control samples that
#'   share the same value with the current sample are excluded from its
#'   comparison. (Default: \code{NULL})
#'   \item \code{k}: \code{Integer scalar}. Number of nearest control
#'   samples used by \code{cloud}. (Default: 5% of the control samples, at
#'   least 1)
#'   \item \code{altexp}: \code{NULL} or \code{Character scalar}. If provided,
#'   use this \code{altExp} assay container of \code{x} for the computation.
#'   (Default: \code{NULL})
#'   \item Further arguments are passed to
#'   \code{\link[=getDissimilarity]{getDissimilarity}}.
#' }
#'
#' @details
#' The following scores are available:
#' \itemize{
#'   \item \code{median_distance}: Median dissimilarity between the sample and
#'   the control samples. The sample itself is excluded from the controls. With
#'   Bray-Curtis dissimilarity, this is the dysbiosis score of Lloyd-Price et
#'   al. (2019).
#'
#'   \item \code{centroid_distance}: Distance from the sample to the centroid
#'   of the control samples minus its distance to the centroid of the other
#'   samples (AlShawaqfeh et al. 2017). Positive values mean that the sample is
#'   closer to the other samples than to the controls. Centroids are computed
#'   from the dissimilarities, and each sample is included in the centroid of
#'   its own group.
#'
#'   \item \code{shannon_jsd}: Mean difference in Shannon diversity between
#'   the sample and the control samples, multiplied by its mean Jensen-Shannon
#'   divergence to them (Santiago et al. 2019). The sample itself is excluded
#'   from the controls. \code{method} is not used.
#'
#'   \item \code{reference_plane}: Euclidean distance, in the first three
#'   principal coordinates of the dissimilarities, from the sample to a plane
#'   fitted to the control samples and restricted to their range (Halfvarson
#'   et al. 2017).
#'
#'   \item \code{cloud}: CLOUD statistic (Montassier et al. 2018): the mean
#'   dissimilarity from the sample to its \code{k} nearest control samples,
#'   divided by the same quantity averaged over the control samples (each
#'   compared with the other controls). Values above 1 mean that the sample is
#'   further from the controls than controls are from each other.
#' }
#'
#' The \code{centroid_distance}, \code{shannon_jsd} and
#' \code{reference_plane} scores are adapted from the dysbiosisR package
#' (Shetty and de Steenhuijsen Piters).
#'
#' @return
#' \code{getDysbiosis} returns a \code{DataFrame} with one row per sample.
#' \code{addDysbiosis} returns \code{x} with the scores added to
#' \code{colData(x)}.
#'
#' @references
#' \insertRef{AlShawaqfeh2017}{mia}
#'
#' \insertRef{Halfvarson2017}{mia}
#'
#' \insertRef{Lloyd-Price2019}{mia}
#'
#' \insertRef{Montassier2018}{mia}
#'
#' \insertRef{Santiago2019}{mia}
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
#' # Other scores
#' res <- getDysbiosis(
#'     tse, col.var = "patient_status", reference = "Control",
#'     index = c("centroid_distance", "reference_plane", "cloud"))
#' head(res)
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
        supported <- c(
            "median_distance", "centroid_distance", "shannon_jsd",
            "reference_plane", "cloud")
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
                median_distance = .calc_median_distance_to_controls,
                centroid_distance = .calc_centroid_distance,
                shannon_jsd = .calc_shannon_jsd,
                reference_plane = .calc_reference_plane_distance,
                cloud = .calc_cloud
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

# Distance to the control centroid minus distance to the centroid of the other
# samples. The distance from sample i to the centroid of group G is
# sqrt(mean_j d_ij^2 - sum_jk d_jk^2 / (2 |G|^2)), for j, k in G.
.calc_centroid_distance <- function(
        mat, groups, method, sample_groups = NULL, altexp = NULL, ...){
    diss <- getDissimilarity(mat |> t(), method = method, ...) |> as.matrix()
    d2 <- diss^2
    is_control <- as.integer(groups) == 1L
    to_centroid <- function(in_group){
        within <- sum(d2[in_group, in_group]) / (2 * sum(in_group)^2)
        sqrt(pmax(rowMeans(d2[, in_group, drop = FALSE]) - within, 0))
    }
    res <- to_centroid(is_control) - to_centroid(!is_control)
    return(unname(res))
}

# Mean difference in Shannon diversity to the controls times the mean
# Jensen-Shannon divergence to them; the sample itself is excluded.
.calc_shannon_jsd <- function(
        mat, groups, sample_groups = NULL, method = NULL, altexp = NULL, ...){
    rel <- t(t(mat) / colSums(mat))
    shannon <- apply(rel, 2, function(p){
        p <- p[p > 0]
        -sum(p * log(p))
    })
    jsd <- getDissimilarity(mat |> t(), method = "jsd") |> as.matrix()
    is_control <- as.integer(groups) == 1L
    res <- vapply(seq_len(ncol(mat)), function(i){
        ref <- is_control
        ref[i] <- FALSE
        mean(shannon[i] - shannon[ref]) * mean(jsd[i, ref])
    }, numeric(1))
    return(res)
}

# Distance from each sample to a plane fitted to the controls in the first
# three principal coordinates. The plane is restricted to the range of the
# controls: the projection of the sample on the plane is clipped to that range
# and moved back onto the plane before the distance is measured.
#' @importFrom stats cmdscale lm coef
.calc_reference_plane_distance <- function(
        mat, groups, method, sample_groups = NULL, altexp = NULL, ...){
    diss <- getDissimilarity(mat |> t(), method = method, ...)
    pc <- cmdscale(diss, k = 3L)
    is_control <- as.integer(groups) == 1L
    ref <- pc[is_control, , drop = FALSE]
    fit <- coef(lm(ref[, 3L] ~ ref[, 1L] + ref[, 2L]))
    # Plane a*x + b*y - z + d = 0
    normal <- c(fit[[2L]], fit[[3L]], -1)
    d <- fit[[1L]]
    ranges <- apply(ref, 2L, range)
    res <- apply(pc, 1L, function(point){
        proj <- point - (sum(normal * point) + d) / sum(normal^2) * normal
        proj <- pmin(pmax(proj, ranges[1L, ]), ranges[2L, ])
        proj[3L] <- d + normal[1L] * proj[1L] + normal[2L] * proj[2L]
        sqrt(sum((point - proj)^2))
    })
    return(unname(res))
}

# CLOUD statistic: mean distance to the k nearest controls, divided by the
# average of the same quantity over the controls. The sample itself is never
# its own neighbour.
.calc_cloud <- function(
        mat, groups, method, k = NULL, sample_groups = NULL, altexp = NULL,
        ...){
    diss <- getDissimilarity(mat |> t(), method = method, ...) |> as.matrix()
    diag(diss) <- NA_real_
    is_control <- as.integer(groups) == 1L
    if( is.null(k) ){
        k <- max(1L, round(0.05 * sum(is_control)))
    }
    if( !(.is_an_integer(k) && k >= 1L && k < sum(is_control)) ){
        stop("'k' must be a positive integer smaller than the number of ",
            "control samples.", call. = FALSE)
    }
    to_controls <- diss[, is_control, drop = FALSE]
    nearest <- apply(to_controls, 1L, function(d) mean(sort(d)[seq_len(k)]))
    res <- nearest / mean(nearest[is_control])
    return(unname(res))
}
