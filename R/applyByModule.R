#' Apply a function by feature or sample group
#'
#' @name applyByModule
#'
#' @param x A \code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}} object.
#'
#' @param by \code{Character scalar}. Whether \code{group} information pertains
#'   the \code{"rows"} or \code{"cols"} of \code{x}.
#'
#' @param group \code{Character vector}. Names of the groups or modules by which
#'   \code{FUN} is applied.
#'
#' @param FUN \code{Function scalar}. A function which takes \code{x} as the
#'   first argument and \code{...} as additional arguments. To update \code{x}
#'   in place, use \code{get} and not \code{add} functions, as the latter will
#'   populate the altExp slot.
#'
#' @param ... Additional arguments passed to \code{FUN}.
#' \itemize{
#'   \item \code{meta.name}: \code{Character scalar}. The name of the metadata
#'   slot where results are stored. (Default: \code{"mod.res"})
#'
#'   \item \code{min.group.size}: \code{Numeric scalar}. The minimum number of
#'   features or samples in a group for \code{FUN} to be applied.
#'   (Default: \code{NULL})
#' }
#'
#' @returns \code{x} updated with results stored in the appropriate slot
#' depending on the given \code{FUN}.
#'
#' @examples
#' library(miaViz)
#'
#' # Import dataset
#' data("Tengeler2020", package = "mia")
#' tse <- Tengeler2020
#'
#' # Compute alpha diversity indices for each Order
#' tse <- applyByModule(
#'     tse,
#'     by = "rows",
#'     group = "Order",
#'     FUN = getAlpha,
#'     index = c("shannon", "faith"),
#'     min.group.size = 2
#' )
#'
#' # Compute beta diversity measures for each Order
#' tse <- applyByModule(
#'     tse,
#'     by = "rows",
#'     group = "Order",
#'     FUN = getMDS,
#'     method = "unifrac",
#'     min.group.size = 2
#' )
#'
#' # Transform assay for each Order, storing results as altExps
#' tse <- applyByModule(
#'     tse,
#'     by = "rows",
#'     group = "Order",
#'     FUN = transformAssay,
#'     method = "clr",
#'     pseudocount = TRUE
#' )
#'
#' rowData(tse)$Var1 <- rowData(tse)$Family == "Ruminococcaceae"
#' rowData(tse)$Var2 <- rowData(tse)$Order == "Clostridiales"
#'
#' # Plot abundance for each module, specifying metadata name
#' tse <- applyByModule(
#'     tse,
#'     by = "rows",
#'     group = c("Var1", "Var2"),
#'     FUN = plotAbundance,
#'     meta.name = "abund_plots"
#' )
#'
NULL

#' @export
#' @rdname applyByModule
applyByModule <- function(x, by, group, FUN, ...) {
    if (length(by) == 0L || !by %in% c("rows", "cols")) {
        stop("'by' must be either 'rows' or 'cols'.", call. = FALSE)
    }
    data_fun <- if (by == "rows") rowData else colData
    if (!.is_non_empty_character(group) || !all(group %in% colnames(data_fun(x)))) {
        stop("All elements of 'group' must be variables of 'x'.", call. = FALSE)
    }

    calc <- .applyByModule_calculate(
        x = x, by = by, group = group, FUN = FUN, ...
    )

    x <- .applyByModule_store(
        x = x,
        res = calc$res,
        by = by,
        group = calc$group,
        FUN = FUN,
        ...
    )

    return(x)
}

#' @importFrom BiocParallel bplapply
.applyByModule_calculate <- function(x, by, group, FUN, ..., min.group.size = NULL) {
    if (!(is.null(min.group.size) || (is.numeric(min.group.size) && length(min.group.size) == 1L))) {
        stop("'min.group.size' must be a single numeric value.", call. = FALSE)
    }
    # Select side information based on margin
    data_fun <- if (by == "rows") rowData else colData

    # Build indices per group
    if (length(group) == 1L) {
        gvar <- data_fun(x)[[group]]
        gvar <- as.factor(gvar)

        group_levels <- levels(gvar)
        gnum <- as.numeric(gvar)

        idx <- lapply(seq_along(group_levels), function(g) which(gnum == g))
        names(idx) <- group_levels
        group <- group_levels
    } else {
        df <- data_fun(x)[group]
        idx <- apply(df, 2L, function(col) which(col != 0))
    }

    # Apply min group size filter
    if (!is.null(min.group.size)) {
        to_keep <- which(lengths(idx, use.names = FALSE) >= min.group.size)
        idx <- idx[to_keep]
        group <- group[to_keep]
    }

    # Compute
    res <- bplapply(idx, function(i) FUN(x[i, ], ...))

    return(list(res = res, group = group))
}

.applyByModule_store <- function(x, res, by, group, FUN, meta.name = "mod_res", ...) {
    if (!.is_a_string(meta.name)) {
        stop("'meta.name' must be a non-empty single character value.",
            call. = FALSE
        )
    }
    fun_name <- FUN |>
        substitute() |>
        deparse()

    kwargs <- list(...)

    if (by == "rows" && inherits(res[[1L]], "SummarizedExperiment")) {
        x <- .applyByModule_store_altExps(x, res)
        return(x)
    }

    if (fun_name == "getAlpha") {
        x <- .applyByModule_store_alpha(x, res, group)
        return(x)
    }

    if (fun_name %in% c("getDivergence", "getDominant")) {
        x <- .applyByModule_store_div_dom(x, res, group, fun_name, kwargs)
        return(x)
    }

    if (fun_name %in% c(
        "getMDS", "getNMDS", "getCCA", "getRDA", "getLDA",
        "getNMF", "getDPCoA", "calculatePCA"
    )) {
        x <- .applyByModule_store_reducedDims(x, res)
        return(x)
    }

    if (fun_name == "getCluster") {
        x <- .applyByModule_store_cluster(x, res, group, kwargs)
        return(x)
    }

    x <- .applyByModule_store_metadata(x, res, meta.name)
    return(x)
}

.applyByModule_store_altExps <- function(x, res) {
    altExps(x) <- c(altExps(x), res)
    return(x)
}

#' @importFrom dplyr bind_cols
.applyByModule_store_alpha <- function(x, res, group) {
    num_index <- ncol(res[[1L]])
    col_labs <- rep(group, each = num_index)

    res_df <- res |>
        lapply(as.data.frame) |>
        bind_cols(.name_repair = "minimal")

    colnames(res_df) <- paste0(colnames(res_df), ".", col_labs)
    colData(x) <- .update_side_info(colData(x), res_df)

    return(x)
}

#' @importFrom dplyr bind_cols
.applyByModule_store_div_dom <- function(x, res, group, fun_name, kwargs) {
    if ("name" %in% names(kwargs)) {
        lab <- kwargs$name
    } else {
        lab <- switch(fun_name,
            getDivergence = "divergence",
            getDominant   = "dominant"
        )
    }

    res_df <- res |>
        lapply(as.data.frame) |>
        dplyr::bind_cols(.name_repair = "minimal")

    colnames(res_df) <- paste0(lab, "_", group)
    colData(x) <- .update_side_info(colData(x), res_df)

    return(x)
}

.applyByModule_store_reducedDims <- function(x, res) {
    reducedDims(x) <- c(reducedDims(x), res)
    return(x)
}

#' @importFrom dplyr bind_cols
.applyByModule_store_cluster <- function(x, res, group, kwargs) {
    res_df <- res |>
        lapply(as.data.frame) |>
        dplyr::bind_cols(.name_repair = "minimal")

    clust_by <- if ("MARGIN" %in% names(kwargs)) kwargs$MARGIN else "rows"
    clust_col <- if ("clust.col" %in% names(kwargs)) kwargs$clust.col else "cluster"

    colnames(res_df) <- paste0(clust_col, ".", group)

    if (clust_by == "rows") {
        rowData(x) <- .update_side_info(rowData(x), res_df)
    } else {
        colData(x) <- .update_side_info(colData(x), res_df)
    }

    return(x)
}

.applyByModule_store_metadata <- function(x, res, meta.name) {
    metadata(x)[[meta.name]] <- res
    return(x)
}

.update_side_info <- function(old, new) {
    to_remove <- which(colnames(old) %in% colnames(new))
    if (length(to_remove) != 0L) {
        warning("Some columns were replaced.", call. = FALSE)
        old[to_remove] <- NULL
    }
    updated <- cbind(old, new)
    return(updated)
}
