#' Apply a function by feature or sample group
#'
#' @name applyByModule
#' 
#' @description
#' applyByModule provides a convenient method to apply the same operation on
#' multiple, possibly overlapping groups of features or samples, representing
#' for example functional microbial modules or sample categories, respectively.
#' 
#' @param x A \code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}} object.
#'
#' @param .by \code{Character scalar}. Whether \code{.group} information pertains
#'   the \code{"rows"} or \code{"cols"} of \code{x}.
#'
#' @param .group \code{Character vector}. Names of the groups or modules by which
#'   \code{FUN} is applied.
#'
#' @param FUN \code{Function scalar}. A function which takes \code{x} as the
#'   first argument and \code{...} as additional arguments. To update \code{x}
#'   in place, use \code{get} and not \code{add} functions, as the latter will
#'   populate the altExp slot.
#'
#' @param min.group.size \code{Numeric scalar}. The minimum number of
#'   features or samples in a group for \code{FUN} to be applied.
#'   (Default: \code{NULL})
#'
#' @param meta.name \code{Character scalar}. The name of the metadata
#'   slot where results are stored. (Default: \code{"mod.res"})
#'
#' @param ... Additional arguments passed to \code{FUN}.
#'
#' @returns \code{x} updated with results stored in the appropriate slot
#' depending on the given \code{FUN}.
#' 
#' @details
#' It is recommended to use \code{get*} functions to store results in the
#' appropriate slot of \code{x}, whereas using the \code{add*} functions applied
#' to features will save the results as separate experiments in the altExp slot.
#' In summary, function results are stored as follows:
#' 
#' \itemize{
#'   \item colData: getAlpha, getDominance, getPrevalence, getCluster
#'   \item reducedDims: getMDS, getNMDS, getCCA, getRDA, getLDA, getNMF, getDPCoA, calculatePCA
#'   \item altExps: transformAssay and add* functions
#'   \item metadata: all else (miaViz ggplot2 objects, statistical summaries, etc.)
#' }
#' 
#' @examples
#' library(bluster)
#' library(miaViz)
#' library(patchwork)
#'
#' # Import dataset
#' data("Tengeler2020", package = "mia")
#' tse <- Tengeler2020
#'
#' # Compute alpha diversity indices for each Order
#' tse <- applyByModule(
#'     tse,
#'     .by = "rows",
#'     .group = "Order",
#'     FUN = getAlpha,
#'     index = c("shannon", "faith"),
#'     min.group.size = 2
#' )
#' 
#' # View colData names for alpha diversity results
#' names(colData(tse))
#' 
#' # Compute beta diversity measures for each Order
#' tse <- applyByModule(
#'     tse,
#'     .by = "rows",
#'     .group = "Order",
#'     FUN = getMDS,
#'     method = "unifrac",
#'     min.group.size = 2
#' )
#' 
#' # View reducedDim names for beta diversity results
#' reducedDimNames(tse)
#' 
#' # Transform assay for each Order, storing results as altExps
#' tse <- applyByModule(
#'     tse,
#'     .by = "rows",
#'     .group = "Order",
#'     FUN = transformAssay,
#'     method = "clr",
#'     pseudocount = TRUE
#' )
#' 
#' # View altExp names for newly transformed assays
#' altExpNames(tse)
#' 
#' # Create two example microbial modules
#' rowData(tse)$Mod1 <- rowData(tse)$Family == "Ruminococcaceae"
#' rowData(tse)$Mod2 <- rowData(tse)$Order == "Clostridiales"
#' 
#' # Find sample-wise most dominant genus for each module
#' tse <- applyByModule(
#'     tse,
#'     .by = "rows",
#'     .group = c("Mod1", "Mod2"),
#'     FUN = getDominant,
#'     group = "Genus"
#' )
#' 
#' # View colData names for dominance results
#' names(colData(tse))
#' 
#' # Cluster samples by features for each module
#' tse <- applyByModule(
#'     tse,
#'     .by = "rows",
#'     .group = c("Mod1", "Mod2"),
#'     FUN = getCluster,
#'     assay.type = "counts",
#'     by = "cols",
#'     BLUSPARAM = KmeansParam(centers = 3)
#' )
#' 
#' # View colData names for clustering results
#' names(colData(tse))
#' 
#' # Plot tree for each module, specifying metadata name
#' tse <- applyByModule(
#'     tse,
#'     .by = "rows",
#'     .group = c("Mod1", "Mod2"),
#'     FUN = plotRowTree,
#'     edge.colour.by = "Genus",
#'     meta.name = "tree_plots"
#' )
#' 
#' # Visualize module-wise abundance plots
#' metadata(tse)$tree_plots |>
#'     wrap_plots()
NULL

#' @export
#' @rdname applyByModule
setMethod("applyByModule", "SummarizedExperiment", function(x, .by, .group, FUN,
    min.group.size = NULL, meta.name = "mod_res", ...){
    if (length(.by) == 0L || !.by %in% c("rows", "cols")) {
        stop("'.by' must be either 'rows' or 'cols'.", call. = FALSE)
    }
    data_fun <- if (.by == "rows") rowData else colData
    if (!.is_non_empty_character(.group) || !all(.group %in% colnames(data_fun(x)))) {
        stop("All elements of '.group' must be variables of 'x'.", call. = FALSE)
    }
    # Deduce function name
    fun_name <- FUN |>
        substitute() |>
        deparse()
    # Calculate group-wise results for function
    calc <- .applyByModule_calculate(
        x = x,
        .by = .by,
        .group = .group,
        FUN = FUN,
        FUN.name = fun_name,
        min.group.size = min.group.size,
        ...
    )
    # Store results in the appropriate slot
    x <- .applyByModule_store(
        x = x,
        res = calc$res,
        .by = .by,
        .group = calc$group,
        FUN = FUN,
        FUN.name = fun_name,
        meta.name = meta.name,
        ...
    )

    return(x)
})

#' @importFrom BiocParallel bplapply
.applyByModule_calculate <- function(x, .by, .group, FUN, FUN.name, min.group.size, ...) {
    if (!(is.null(min.group.size) || (is.numeric(min.group.size) && length(min.group.size) == 1L))) {
        stop("'min.group.size' must be a single numeric value.", call. = FALSE)
    }
    # Select side information based on margin
    data_fun <- if (.by == "rows") rowData else colData

    # Build indices per group
    if (length(.group) == 1L) {
        # Retrieve grouping variable
        gvar <- data_fun(x)[[.group]]
        # Factorise grouping variable
        gvar <- as.factor(gvar)
        # Find levels of grouping variable
        group_levels <- levels(gvar)
        gnum <- as.numeric(gvar)
        # Find group membership for each feature
        idx <- lapply(seq_along(group_levels), function(g) which(gnum == g))
        names(idx) <- group_levels
        .group <- group_levels
    } else {
        # Retrieve grouping variables
        df <- data_fun(x)[.group]
        # Find group memberships for each feature
        idx <- apply(df, 2L, function(col) which(col != 0))
    }

    # Apply min group size filter
    if (!is.null(min.group.size)) {
        to_keep <- which(lengths(idx, use.names = FALSE) >= min.group.size)
        idx <- idx[to_keep]
        .group <- .group[to_keep]
    }
    
    if (length(idx) == 0L) {
        stop("No '.group' has size above 'min.group.size'.", call. = FALSE)
    }

    # Compute
    res <- bplapply(idx, function(i) FUN(x[i, ], ...))

    return(list(res = res, group = .group))
}

.applyByModule_store <- function(x, res, .by, .group, FUN, FUN.name, meta.name, ...) {
    if (!.is_a_string(meta.name)) {
        stop("'meta.name' must be a non-empty single character value.",
            call. = FALSE
        )
    }
    # Retrieve named arguments for function
    kwargs <- list(...)
    
    if (.by == "rows" && inherits(res[[1L]], "SummarizedExperiment")) {
        # Store add* results in altExp
        altExps(x) <- .applyByModule_update_slot(altExps(x), res)
    } else if (FUN.name == "getAlpha") {
        # Store alpha results in colData
        x <- .applyByModule_store_alpha(x, res, .group)
    } else if (FUN.name %in% c("getDivergence", "getDominant")) {
        # Store divergence and dominance results in colData
        x <- .applyByModule_store_div_dom(x, res, .group, FUN.name, kwargs)
    } else if (FUN.name %in% c(
        "getMDS", "getNMDS", "getCCA", "getRDA", "getLDA",
        "getNMF", "getDPCoA", "calculatePCA"
    )) {
        # Store beta diversity results in reducedDims
        reducedDims(x) <- .applyByModule_update_slot(reducedDims(x), res)
    } else if (FUN.name == "getCluster") {
        # Store clustering results in colData
        x <- .applyByModule_store_cluster(x, res, .group, kwargs)
    }else{
        # Store other results in meta data
        metadata(x)[[meta.name]] <- res
    }
    return(x)
}

#' @importFrom dplyr bind_cols
.applyByModule_store_alpha <- function(x, res, .group) {
    # Find number of indices
    num_index <- ncol(res[[1L]])
    # Replicate index names
    col_labs <- rep(.group, each = num_index)
    # Bind resulting data frames
    res_df <- res |>
        lapply(as.data.frame) |>
        bind_cols(.name_repair = "minimal")
    # Prepend group names with metric labels
    colnames(res_df) <- paste0(colnames(res_df), ".", col_labs)
    # Update side information with new results
    colData(x) <- .applyByModule_update_slot(colData(x), res_df)

    return(x)
}

#' @importFrom dplyr bind_cols
.applyByModule_store_div_dom <- function(x, res, .group, FUN.name, kwargs) {
    # Retrieve result column label
    if ("name" %in% names(kwargs)) {
        lab <- kwargs$name
    } else {
        lab <- switch(FUN.name,
            getDivergence = "divergence",
            getDominant   = "dominant"
        )
    }
    # Bind resulting data frames
    res_df <- res |>
        lapply(as.data.frame) |>
        bind_cols(.name_repair = "minimal")
    # Prepend group names with metric label
    colnames(res_df) <- paste0(lab, ".", .group)
    # Update side information with new results
    colData(x) <- .applyByModule_update_slot(colData(x), res_df)

    return(x)
}

#' @importFrom dplyr bind_cols
.applyByModule_store_cluster <- function(x, res, .group, kwargs) {
    # Bind resulting data frames
    res_df <- res |>
        lapply(as.data.frame) |>
        bind_cols(.name_repair = "minimal")
    # Retrieve margin and clust.col args
    clust_by <- if ("by" %in% names(kwargs)) kwargs$by else "rows"
    clust_col <- if ("clust.col" %in% names(kwargs)) kwargs$clust.col else "cluster"
    # Prepend group names with clust.col
    colnames(res_df) <- paste0(clust_col, ".", .group)
    # Update side information with new results
    if (clust_by == "rows") {
        rowData(x) <- .applyByModule_update_slot(rowData(x), res_df)
    } else {
        colData(x) <- .applyByModule_update_slot(colData(x), res_df)
    }

    return(x)
}

.applyByModule_update_slot <- function(old, new) {
    # Find variables to replace
    to_remove <- which(names(old) %in% names(new))
    # Replace old variables with new ones
    if (length(to_remove) != 0L) {
        warning("Some elements were replaced.", call. = FALSE)
        old[to_remove] <- NULL
    }
    # Select bind function based on info type
    bind <- switch(class(old), DFrame = cbind, SimpleList = c)
    # Bind unique variables together
    updated <- bind(old, new)
    return(updated)
}
