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
#'   \code{SummarizedExperiment::assayNames()} (or index \code{1L} if unnamed).
#'   The actually used assay names are recorded in \code{$assay.names.used} in
#'   the result. If you need a different assay (e.g. \code{"relab"} instead of
#'   \code{"counts"}), subset or reorder assays in \code{x} before calling
#'   \code{jointRPCAuniversal()}.
#' @param ... Additional arguments passed to \code{.joint_rpca()}.
#'
#' @return Output from \code{.joint_rpca()} with extra fields when \code{x} is a
#'   \code{MultiAssayExperiment}:
#'   \itemize{
#'     \item \code{$experiment.names}: character vector of experiments used.
#'     \item \code{$assay.names.used}: named character vector giving, for each
#'           experiment, the assay name that was used (typically the first in
#'           \code{assayNames()}).
#'   }
#' @export

jointRPCAuniversal <- function(x, experiments = NULL, ...) {
    assay_names_used <- NULL
    
    if (inherits(x, "MultiAssayExperiment")) {
        
        exps <- MultiAssayExperiment::experiments(x)
        
        if (is.null(experiments)) {
            experiments <- names(exps)
        }
        if (length(experiments) == 0L) {
            stop("No experiments found in 'x'.")
        }
        
        assay_names_used <- setNames(
            character(length(experiments)),
            experiments
        )
        
        tables <- vector("list", length(experiments))
        names(tables) <- experiments
        
        for (i in seq_along(experiments)) {
            e <- experiments[[i]]
            exp_se <- exps[[e]]
            if (is.null(exp_se)) {
                stop(sprintf("Experiment '%s' not found in 'x'.", e))
            }
            
            anm <- SummarizedExperiment::assayNames(exp_se)
            if (length(anm)) {
                default_assay <- anm[1]
            } else {
                default_assay <- 1L
            }
            
            assay_names_used[[e]] <- if (is.character(default_assay)) {
                default_assay
            } else {
                as.character(default_assay)
            }
            
            tables[[e]] <- SummarizedExperiment::assay(exp_se, default_assay)
        }
        
        names(tables) <- experiments
        
    } else if (inherits(x, "SummarizedExperiment")) {
        #covers TreeSummarizedExperiment as it inherits from SE
        tables <- list(SummarizedExperiment::assay(x))
        nm <- tryCatch(SummarizedExperiment::assayNames(x)[1], error = function(e) NULL)
        names(tables) <- if (is.null(nm) || is.na(nm)) "assay1" else nm
        
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
            paste(class(x), collapse = ", ")
        )
    }
    
    res <- .joint_rpca(tables = tables, ...)
    
    if (inherits(x, "MultiAssayExperiment")) {
        res$experiment.names <- experiments
        res$assay.names.used <- assay_names_used
    }
    
    res
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
#' @param ... Additional arguments passed to \code{jointRPCAuniversal()},
#'   and from there to the internal \code{.joint_rpca()} engine (e.g. \code{n.components},
#'   \code{min.sample.count}, \code{max.iterations}, etc.).
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
                         ...) {
    #select the object to operate on
    y <- x
    if (!is.null(altexp)) {
        y <- altExp(x, altexp)
    }
    
    #use universal front-end to build tables + run .joint_rpca()
    res <- jointRPCAuniversal(y, experiments = experiments, ...)
    
    #extract sample embedding
    emb <- res$ord.res$samples
    if (is.null(emb) || !is.matrix(emb)) {
        stop("Joint-RPCA did not return a valid sample embedding in ord.res$samples.")
    }
    
    #store embedding in reducedDim only if supported (SCE / TreeSE / mia-specific)
    cls <- class(x)
    if (any(cls %in% c("SingleCellExperiment", "TreeSummarizedExperiment"))) {
        reducedDim(x, name) <- emb
    }
    
    #store full result in metadata
    if (is.null(metadata(x)$JointRPCA)) {
        metadata(x)$JointRPCA <- list()
    }
    metadata(x)$JointRPCA[[name]] <- res
    
    x
}