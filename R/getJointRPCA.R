## Joint RPCA front-end helpers
##
## Exported:
##   - jointRPCAuniversal()
##   - getJointRPCA()
##

#' Universal Joint RPCA Wrapper
#'
#' @param x Input object: \code{MultiAssayExperiment}, \code{SummarizedExperiment}
#'   (including \code{TreeSummarizedExperiment}), list of matrices, or single matrix.
#' @param experiments Character vector of experiment names to extract when \code{x}
#'   is a \code{MultiAssayExperiment} (i.e. \code{names(experiments(x))}).
#'   If \code{NULL}, all experiments are used.
#'
#'   For \code{MultiAssayExperiment} inputs, \strong{one assay per experiment} is
#'   used: by default the first assay returned by
#'   \code{assayNames()} (or index \code{1L} if unnamed).
#'   The actually used assay names are recorded in \code{$assay_names_used} in
#'   the result. If you need a different assay (e.g. \code{"relab"} instead of
#'   \code{"counts"}), subset or reorder assays in \code{x} before calling
#'   \code{jointRPCAuniversal()}.
#' @param transform Character string specifying preprocessing applied to each
#'   input table before ordination. Use \code{"rclr"} to apply the robust CLR
#'   transform (via \code{decostand(method = "rclr")}) or \code{"none"} to
#'   disable transformation (data are used as-is after masking non-finite values).
#' @param optspace.tol Numeric tolerance passed to \code{optspace()}.
#' @param center Logical; whether to center the reconstructed low-rank matrix
#'   (double-centering) prior to SVD/PCA steps.
#' @param scale Logical; whether to scale the reconstructed matrix prior to
#'   SVD/PCA steps. Defaults to \code{FALSE}.
#' @param ... Additional arguments passed to \code{.joint_rpca()}.
#'
#' @return Output from \code{.joint_rpca()} with extra fields when \code{x} is a
#'   \code{MultiAssayExperiment}:
#'   \itemize{
#'     \item \code{$experiment_names}: character vector of experiments used.
#'     \item \code{$assay_names_used}: named character vector giving, for each
#'           experiment, the assay name that was used (typically the first in
#'           \code{assayNames()}).
#'   }
#' @importFrom SummarizedExperiment assayNames
#' @importFrom SummarizedExperiment assay
#' @importFrom MultiAssayExperiment experiments
#' @importFrom vegan decostand
#' @importFrom vegan optspace
#' @export

jointRPCAuniversal <- function(x, experiments = NULL,
                               transform = c("rclr", "none"),
                               optspace.tol = 1e-5,
                               center = TRUE,
                               scale = FALSE,
                               ...) {
    
    transform <- match.arg(transform)
    
    assay_names_used <- NULL
    
    if (inherits(x, "MultiAssayExperiment")) {
        mae <- .extract_mae_tables(x, experiments)
        tables <- mae$tables
        experiments <- mae$experiments
        assay_names_used <- mae$assay_names_used
        
    } else if (inherits(x, "SummarizedExperiment")) {
        tables <- list(assay(x))
        anm <- assayNames(x)
        nm <- if (length(anm) && !is.na(anm[1])) anm[1] else "assay1"
        names(tables) <- nm
        
    } else if (is.list(x) && all(vapply(x, is.matrix, logical(1)))) {
        tables <- x
        if (is.null(names(tables))) {
            names(tables) <- paste0("view", seq_along(tables))
        }
        
    } else if (is.matrix(x)) {
        tables <- list(x)
        names(tables) <- "assay1"
        
    } else {
        stop(
            "Unsupported input type for jointRPCAuniversal(): ",
            paste(class(x), collapse = ", "),
            call. = FALSE
        )
    }
    
    res <- .joint_rpca(
        tables = tables,
        transform = transform,
        optspace.tol = optspace.tol,
        center = center,
        scale = scale,
        ...
    )
    
    if (inherits(x, "MultiAssayExperiment")) {
        res$experiment_names <- experiments
        res$assay_names_used <- assay_names_used
    }
    
    return(res)
}

#' Run Joint-RPCA and store embedding in reducedDim
#'
#' Convenience wrapper that runs Joint Robust PCA on one or more compositional
#' tables and stores the resulting sample embedding in \code{reducedDim(x, name)},
#' similar to \code{runMDS()} and \code{runPCA()}.
#'
#' @param x A \code{SummarizedExperiment}, \code{TreeSummarizedExperiment},
#'   \code{MultiAssayExperiment}, or a related object supported by
#'   \code{jointRPCAuniversal()}.
#' @param experiments Optional character vector of experiment names to use when
#'   \code{x} is a \code{MultiAssayExperiment} (i.e. \code{names(experiments(x))}).
#'   Ignored for \code{SummarizedExperiment} inputs.
#' @param altexp Optional name of an alternative experiment. If supplied,
#'   Joint-RPCA is run on \code{altExp(x, altexp)} instead of \code{x}.
#' @param name Character scalar giving the name of the \code{reducedDim} slot
#'   in which to store the joint sample embedding. Defaults to \code{"JointRPCA"}.
#' @param transform Character string specifying preprocessing applied to each
#'   input table before ordination. Use \code{"rclr"} to apply the robust CLR
#'   transform (via \code{decostand(method = "rclr")}) or \code{"none"} to
#'   disable transformation (data are used as-is after masking non-finite values).
#' @param optspace.tol Numeric tolerance passed to \code{optspace()}.
#' @param center Logical; whether to center the reconstructed low-rank matrix
#'   (double-centering) prior to SVD/PCA steps.
#' @param scale Logical; whether to scale the reconstructed matrix prior to
#'   SVD/PCA steps. Defaults to \code{FALSE}.
#' @param ... Additional arguments passed to \code{jointRPCAuniversal()} and then
#'   to the internal \code{.joint_rpca()} engine (e.g. \code{n.components},
#'   \code{min.sample.count}, \code{min.feature.count}, \code{min.feature.frequency},
#'   \code{max.iterations}, \code{sample.metadata}).
#'
#' @return The input object \code{x} with a new entry in
#'   \code{reducedDim(x, name)} containing the Joint-RPCA sample embedding.
#'   The full Joint-RPCA result (including distances, cross-validation
#'   statistics and transformed tables) is stored in
#'   \code{metadata(x)$JointRPCA[[name]]}.
#'
#' @export

getJointRPCA <- function(x,
                         experiments = NULL,
                         altexp = NULL,
                         name   = "JointRPCA",
                         transform = c("rclr", "none"),
                         optspace.tol = 1e-5,
                         center = TRUE,
                         scale = FALSE,
                         ...) {
    
    transform <- match.arg(transform)
    
    # Select the object to operate on
    y <- x
    if (!is.null(altexp)) {
        y <- altExp(x, altexp)
    }
    
    # Use universal front-end to build tables + run .joint_rpca()
    res <- jointRPCAuniversal(
        y,
        experiments = experiments,
        transform   = transform,
        optspace.tol = optspace.tol,
        center      = center,
        scale       = scale,
        ...
    )
    
    # Extract sample embedding
    emb <- res[["ord_res"]][["samples"]]
    
    if (is.null(emb)) {
        stop(
            "Internal error: JointRPCA did not return a sample embedding. ",
            "Please report this and include sessionInfo().",
            call. = FALSE
        )
    }
    
    emb <- as.matrix(emb)
    
    # Ensure embedding rownames match colnames of the target object (Bioconductor requirement)
    target_cols <- colnames(x)
    if (!is.null(altexp)) {
        target_cols <- colnames(x)
    }
    
    if (is.null(target_cols)) {
        stop("Cannot store reducedDim: 'x' has no colnames().", call. = FALSE)
    }
    
    if (is.null(rownames(emb))) {
        stop("Cannot store reducedDim: embedding has no rownames().", call. = FALSE)
    }
    
    # Require the same set of samples/cells
    if (!setequal(rownames(emb), target_cols)) {
        missing_in_emb <- setdiff(target_cols, rownames(emb))
        extra_in_emb   <- setdiff(rownames(emb), target_cols)
        stop(
            "Cannot store reducedDim: embedding rownames do not match colnames(x).\n",
            "Missing in embedding: ", paste(missing_in_emb, collapse = ", "), "\n",
            "Extra in embedding: ", paste(extra_in_emb, collapse = ", "),
            call. = FALSE
        )
    }
    
    # Reorder embedding to exactly match colnames(x)
    emb <- emb[target_cols, , drop = FALSE]
    rownames(emb) <- target_cols
    
    if (nrow(emb) == 0L || ncol(emb) == 0L || is.null(rownames(emb))) {
        stop(
            "Internal error: JointRPCA returned an invalid sample embedding. ",
            "Please report this and include sessionInfo().",
            call. = FALSE
        )
    }
    
    # Store embedding in reducedDim only if supported (SCE / TreeSE / mia-specific)
    cls <- class(x)
    if (any(cls %in% c("SingleCellExperiment", "TreeSummarizedExperiment"))) {
        reducedDim(x, name) <- emb
    }
    
    # Store full result in metadata
    if (is.null(metadata(x)$JointRPCA)) {
        metadata(x)$JointRPCA <- list()
    }
    metadata(x)$JointRPCA[[name]] <- res
    
    return(x)
}