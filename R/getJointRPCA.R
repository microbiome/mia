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

#' Joint Robust PCA on Multiple Compositional Tables
#'
#' Internal engine for Joint Robust Principal Component Analysis (RPCA) using
#' OptSpace on multiple compositional tables.
#'
#' This function assumes a list of already extracted tables and is typically
#' called via \code{jointRPCAuniversal()} or \code{getJointRPCA()}.
#'
#' @param tables A list of compositional data tables (matrices or data frames).
#' @param n.test.samples Integer specifying the number of samples to hold out for testing
#'   (only used if \code{sample.metadata} is \code{NULL}). Default is 10.
#' @param sample.metadata Optional data frame containing sample-level metadata.
#' @param train.test.column The name of the column in \code{sample.metadata}
#'   that defines training vs test samples.
#' @param n.components Integer specifying the number of principal components to compute.
#' @param transform Character string specifying preprocessing applied to each
#'   input table before ordination: \code{"rclr"} or \code{"none"}.
#' @param optspace.tol Numeric tolerance passed to \code{vegan::optspace()}.
#' @param center,scale Logical; whether to center/scale the reconstructed matrix
#'   prior to SVD/PCA steps.
#' @param min.sample.count Minimum total count required for a sample to be retained.
#' @param min.feature.count Minimum total count required for a feature to be retained.
#' @param min.feature.frequency Minimum percentage (0–100) of samples in which a
#'   feature must be non-zero to be retained.
#' @param max.iterations Maximum number of optimization iterations.
#'
#' @return A list with \code{ord_res}, \code{dist}, \code{cv_stats}, and
#'   \code{rclr_tables}.
#'
#' @keywords internal
#' @noRd

.joint_rpca <- function(tables,
                        n.test.samples = 10,
                        sample.metadata = NULL,
                        train.test.column = NULL,
                        n.components = 3,
                        transform = c("rclr", "none"),
                        min.sample.count = 0,
                        min.feature.count = 0,
                        min.feature.frequency = 0,
                        max.iterations = 5,
                        optspace.tol = 1e-5,
                        center = TRUE,
                        scale = FALSE) {
    
    transform <- match.arg(transform)
    
    if (is.null(names(tables))) {
        names(tables) <- paste0("view", seq_along(tables))
    }
    
    if (n.components < 2) {
        stop("n.components must be at least 2.", call. = FALSE)
    }
    if (max.iterations < 1) {
        stop("max.iterations must be at least 1.", call. = FALSE)
    }
    
    # Filtering (always done, independent of rclr)
    tables <- lapply(tables, function(tbl) {
        out <- .rpca_table_processing(
            tbl,
            min.sample.count      = min.sample.count,
            min.feature.count     = min.feature.count,
            min.feature.frequency = min.feature.frequency
        )
        if (nrow(out) == 0L || ncol(out) == 0L) {
            stop(
                "Filtering removed all data in at least one table (0 features or 0 samples). ",
                "Try relaxing filtering thresholds (min.sample.count / min.feature.count / ",
                "min.feature.frequency) or check your input.",
                call. = FALSE
            )
        }
        out
    })
    
    # Find shared samples across views
    sample.sets <- lapply(tables, colnames)
    shared.all.samples <- Reduce(intersect, sample.sets)
    if (length(shared.all.samples) == 0L) {
        stop(
            "No samples overlap between all tables. ",
            "Check that colnames are consistent across views, or",
            "if you are using pre-transformed tables set transform = 'none'.",
            call. = FALSE
        )
    }
    if (length(shared.all.samples) < (n.components + 1L)) {
        stop(
            "Too few shared samples across all tables after filtering (",
            length(shared.all.samples), "). Need at least n.components + 1 shared samples. ",
            "Try lowering n.components or relaxing filtering thresholds.",
            call. = FALSE
        )
    }
    unshared.samples <- setdiff(unique(unlist(sample.sets)), shared.all.samples)
    if (length(unshared.samples) > 0) {
        warning(sprintf("Removing %d sample(s) that do not overlap in tables.", length(unshared.samples)))
    }
    
    # Restrict each table to the shared sample set
    tables <- lapply(tables, function(tbl) {
        tbl[, shared.all.samples, drop = FALSE]
    })
    shared.all.samples <- Reduce(intersect, lapply(tables, colnames))
    
    # Transform tables: rCLR or masking
    rclr_tables <- lapply(tables, function(tbl) {
        mat <- as.matrix(tbl)
        rown <- rownames(mat)
        coln <- colnames(mat)
        
        if (transform == "rclr") {
            mat[!is.finite(mat)] <- 0
            mat[mat < 0]         <- 0
            out <- vegan::decostand(mat, method = "rclr", MARGIN = 2)
            dimnames(out) <- list(rown, coln)
            out
        } else {
            out <- .mask_value_only(mat)$data
            dimnames(out) <- list(rown, coln)
            out
        }
    })
    names(rclr_tables) <- names(tables)
    
    # Determine train/test split
    if (!is.null(sample.metadata) && !is.null(train.test.column)) {
        md <- as.data.frame(sample.metadata)
        md <- md[shared.all.samples, , drop = FALSE]
        train.samples <- rownames(md)[md[[train.test.column]] == "train"]
        test.samples  <- rownames(md)[md[[train.test.column]] == "test"]
    } else {
        ord.tmp <- .optspace_helper(
            rclr.table      = t(rclr_tables[[1]]),
            feature.ids     = rownames(rclr_tables[[1]]),
            sample.ids      = colnames(rclr_tables[[1]]),
            n.components    = n.components,
            max.iterations  = max.iterations,
            tol             = optspace.tol,
            center          = center,
            scale           = scale
        )$ord_res
        sorted.ids <- rownames(ord.tmp$samples[order(ord.tmp$samples[, 1]), ])
        idx <- round(seq(1, length(sorted.ids), length.out = n.test.samples))
        test.samples <- sorted.ids[idx]
        train.samples <- setdiff(shared.all.samples, test.samples)
    }
    
    # Run joint OptSpace
    result <- .joint_optspace_helper(
        tables         = rclr_tables,
        n.components   = n.components,
        max.iterations = max.iterations,
        test.samples   = test.samples,
        train.samples  = train.samples,
        sample.order   = shared.all.samples,
        tol            = optspace.tol,
        center         = center,
        scale          = scale
    )
    
    return(list(
        ord_res      = result$ord_res,
        dist         = result$dist,
        cv_stats     = result$cv_stats,
        rclr_tables  = rclr_tables
    ))
}

#' Joint RPCA Ordination Across Multiple Compositional Tables
#'
#' Internal function that performs Robust PCA via joint OptSpace decomposition across multiple compositional tables.
#' It splits each table into train/test sets, applies joint factorization, reconstructs sample and feature embeddings,
#' optionally projects test samples, computes a sample distance matrix, and returns cross-validation error statistics.
#'
#' @param tables A list of compositional matrices or data frames with features as rows and samples as columns.
#' @param n.components Number of principal components to compute.
#' @param max.iterations Maximum number of optimization iterations for OptSpace.
#' @param test.samples Character vector of sample IDs to be projected into the ordination space.
#' @param train.samples Character vector of sample IDs used to fit the ordination.
#'
#' @return A list with:
#' \describe{
#'   \item{ord_res}{An \code{OrdinationResults} object containing embeddings, loadings, and variance explained.}
#'   \item{dist}{A \code{DistanceMatrix} object for sample embeddings.}
#'   \item{cv_stats}{A data frame summarizing reconstruction error across iterations and tables.}
#' }
#'
#' @keywords internal
#' @noRd

.joint_optspace_helper <- function(tables,
                                   n.components,
                                   max.iterations,
                                   test.samples,
                                   train.samples,
                                   sample.order = NULL,
                                   tol = 1e-5,
                                   center = TRUE,
                                   scale = FALSE) {
    
    # Coerce to matrices and enforce colnames presence
    tables <- lapply(tables, function(tbl) {
        mat <- as.matrix(tbl)
        if (is.null(colnames(mat))) {
            stop("[.joint_optspace_helper] Input table is missing column names (sample IDs).")
        }
        mat
    })
    
    # Global set of samples present in all views
    all_samples <- Reduce(intersect, lapply(tables, colnames))
    
    # Align train/test to actually available samples
    test.samples  <- intersect(test.samples,  all_samples)
    train.samples <- intersect(train.samples, all_samples)
    
    if (!length(test.samples) || !length(train.samples)) {
        stop("[.joint_optspace_helper] Empty train/test split after aligning sample IDs.")
    }
    
    # Split and transpose training/test data per table
    tables.split <- lapply(tables, function(tbl) {
        list(
            t(tbl[, test.samples,  drop = FALSE]),
            t(tbl[, train.samples, drop = FALSE])
        )
    })
    
    # Format input for solver
    tables.for.solver <- lapply(tables.split, function(pair) {
        lapply(pair, as.matrix)
    })
    
    # Run joint OptSpace solver
    opt.result <- .joint_optspace_solve(
        train.test.pairs = tables.for.solver,
        n.components     = n.components,
        max.iter         = max.iterations,
        tol              = tol
    )
    
    U <- opt.result$U
    S <- opt.result$S
    V_list <- opt.result$V_list
    dists <- opt.result$dists
    
    # Assign row/column names to loadings
    pc.names <- paste0("PC", seq_len(n.components))
    
    # Combine feature loadings with table-derived row names
    vjoint <- do.call(rbind, Map(function(tbl, V) {
        rownames(V) <- rownames(tbl)
        colnames(V) <- pc.names
        V
    }, tables, V_list))
    
    U <- U[seq_along(train.samples), , drop = FALSE]
    rownames(U) <- train.samples
    colnames(U) <- pc.names
    
    # Recenter & re-factor via SVD
    X <- U %*% S %*% t(vjoint)
    
    if (center) {
        X <- sweep(X, 2, colMeans(X))
        X <- sweep(X, 1, rowMeans(X))
    }
    if (scale) {
        X <- scale(X, center = FALSE, scale = TRUE)
    }
    svd.res <- svd(X)
    u <- svd.res$u[, seq_len(n.components), drop = FALSE]
    v <- svd.res$v[, seq_len(n.components), drop = FALSE]
    s.eig <- svd.res$d[seq_len(n.components)]
    
    rownames(u) <- train.samples
    rownames(v) <- rownames(vjoint)
    pc.names <- paste0("PC", seq_len(n.components))
    colnames(u) <- colnames(v) <- pc.names
    
    # Build a named per-view features list
    features_list <- lapply(seq_along(tables), function(i) {
        rid <- rownames(tables[[i]])          
        v[rid, , drop = FALSE]                
    })
    names(features_list) <- names(tables)    
    
    prop.exp <- s.eig^2 / sum(s.eig^2)
    ord_res <- .ordination_results(
        method = "rpca",
        eigvals = setNames(s.eig, pc.names),
        samples = u,
        features = features_list,              
        proportion.explained = setNames(prop.exp, pc.names)
    )
    
    # Project test samples
    if (length(test.samples) > 0) {
        test.matrices <- lapply(tables, function(tbl) tbl[, test.samples, drop = FALSE])
        names(test.matrices) <- names(tables)  
        ord_res <- .transform(ord_res, test.matrices, apply.rclr = FALSE)
    }
    
    # Compute distance matrix and CV error summary
    dist.base <- as.matrix(dist(ord_res$samples))
    
    if (!is.null(sample.order)) {
        order_use <- intersect(sample.order, rownames(dist.base))
        dist.mat  <- dist.base[order_use, order_use, drop = FALSE]
    } else {
        dist.mat  <- dist.base
        order_use <- rownames(dist.base)
    }
    
    dist.res <- .distance_matrix(dist.mat, ids = order_use)
    
    cv.dist <- data.frame(t(dists))
    colnames(cv.dist) <- c("mean_CV", "std_CV")
    cv.dist$run <- sprintf("tables_%d.n.components_%d.max.iterations_%d.n.test_%d",
                           length(tables), n.components, max.iterations, length(test.samples))
    cv.dist$iteration <- seq_len(nrow(cv.dist))
    rownames(cv.dist) <- seq_len(nrow(cv.dist))
    
    return(list(ord_res = ord_res, dist = dist.res, cv_stats = cv.dist))
}

#' Apply Projection of New Compositional Tables to Existing Ordination
#'
#' Internal function that transforms and projects new sample tables into an existing Joint RPCA ordination space.
#' It handles feature alignment, optional rCLR preprocessing, padding of missing features,
#' and merges projected samples with existing ones.
#'
#' @param ordination A list containing previous ordination results: `samples`, `features`, and `eigvals`.
#' @param tables A named list of new compositional tables (matrices or data frames)
#'   with features as rows and samples as columns.
#' @param apply.rclr Logical; whether to apply rCLR transformation to new input tables before projection. Default is `TRUE`.
#'
#' @return An updated ordination list with the `samples` matrix extended to include new projected samples.
#' @keywords internal
#' @noRd

.transform <- function(ordination, tables,
                       apply.rclr = TRUE) {
    
    Udf    <- ordination$samples
    Vobj   <- ordination$features
    s.eig  <- ordination$eigvals
    
    # Ensure tables is a list of views
    if (!is.list(tables)) {
        stop("[.transform] 'tables' must be a list of view matrices (features x samples).")
    }
    
    if (is.list(Vobj) && !is.null(names(Vobj))) {
        if (is.null(names(tables))) {
            names(tables) <- names(Vobj)[seq_along(tables)]
        }
        if (is.null(names(tables))) {
            stop("[.transform] 'tables' must be a *named* list of view matrices (features x samples).")
        }
    }
    
    # rCLR if requested 
    prep_view <- function(tab) {
        mat <- as.matrix(tab)
        rown <- rownames(mat)
        coln <- colnames(mat)
        
        if (apply.rclr) {
            mat <- vegan::decostand(mat, method = "rclr", MARGIN = 2)
            
            dimnames(mat) <- list(rown, coln)
        }
        
        storage.mode(mat) <- "double"
        mat[!is.finite(mat)] <- 0
        mat
    }
    tables <- lapply(tables, prep_view)
    
    if (is.matrix(Vobj)) {
        all.features <- rownames(Vobj)
        tables <- lapply(tables, function(mat) {
            miss <- setdiff(all.features, rownames(mat))
            if (length(miss)) {
                pad <- matrix(0, nrow = length(miss), ncol = ncol(mat),
                              dimnames = list(miss, colnames(mat)))
                mat <- rbind(mat, pad)
            }
            mat[all.features, , drop = FALSE]
        })
        proj.mat <- do.call(cbind, tables)
        colnames(proj.mat) <- make.unique(colnames(proj.mat), sep = "_")
        ordination$samples <- .transform_helper(Udf, Vobj, s.eig, proj.mat)
        return(ordination)
    }
    
    # 3+-omic path: V is a named list per view
    if (!is.list(Vobj) || is.null(names(Vobj))) {
        stop("[.transform] ordination$features is neither a matrix nor a named list.")
    }
    
    # Intersect views by name, preserve training order
    views <- intersect(names(Vobj), names(tables))
    if (!length(views)) stop("[.transform] No overlapping view names between ordination and new tables.")
    
    test.matrices <- list()
    for (vw in views) {
        Vvw <- Vobj[[vw]]
        stopifnot(is.matrix(Vvw), !is.null(rownames(Vvw)))
        mat <- tables[[vw]]
        
        train_feats <- rownames(Vvw)
        miss <- setdiff(train_feats, rownames(mat))
        if (length(miss)) {
            pad <- matrix(0, nrow = length(miss), ncol = ncol(mat),
                          dimnames = list(miss, colnames(mat)))
            mat <- rbind(mat, pad)
        }
        mat <- mat[train_feats, , drop = FALSE]
        
        test.matrices[[vw]] <- mat
    }
    
    ordination$samples <- .transform_helper(Udf, Vobj, s.eig, test.matrices)
    return(ordination)
}

#' Project New Data into Existing Ordination Space
#'
#' Internal function to align rCLR-transformed samples to an existing RPCA ordination space.
#' Handles feature alignment, deduplication of sample names, double-centering normalization,
#' and projection into low-rank space using previously learned components.
#'
#' @param Udf Matrix of training sample embeddings (samples × components).
#' @param Vdf Matrix or list of feature loadings (features × components, or per-view list).
#' @param s.eig Singular values from the RPCA decomposition.
#' @param table.rclr.project New rCLR-transformed table(s) for projection (features × samples).
#'   If `Vdf` is a matrix, provide a single matrix; if `Vdf` is a list, provide a named list of matrices per view.
#' @param dedup.samples Logical; whether to merge samples with identical names (e.g. suffixes like `_1`, `_2`)
#'   by averaging their feature values. Default is `TRUE`. Set to `FALSE` to preserve all duplicate sample IDs.
#'
#' @return A combined matrix of training and projected samples (samples × components).
#' @keywords internal
#' @noRd

.transform_helper <- function(Udf, Vdf, s.eig, table.rclr.project,
                              dedup.samples = TRUE) {
    
    # Legacy path (single view)
    if (is.matrix(Vdf)) {
        stopifnot(is.matrix(table.rclr.project))
        # Align rows by name
        common <- intersect(rownames(Vdf), rownames(table.rclr.project))
        if (length(common) < ncol(Udf))
            stop(sprintf("[.transform_helper] Too few matching features: %d", length(common)))
        
        M <- t(as.matrix(table.rclr.project[common, , drop = FALSE]))   
        V <- as.matrix(Vdf[common, , drop = FALSE])                     
        
        # Dedup of sample IDs
        if (dedup.samples) {
            sid <- sub("_\\d+$", "", rownames(M))
            if (any(duplicated(sid))) {
                M <- rowsum(M, group = sid, reorder = FALSE) / as.vector(table(sid))
            } else {
                rownames(M) <- sid
            }
        }
        
        # Projection (match training scaling)
        Uproj <- M %*% V
        # Scale by singular values
        if (length(s.eig)) {
            Sinv <- diag(1 / s.eig, nrow = length(s.eig))
            Uproj <- Uproj %*% Sinv
        }
        
        colnames(Uproj) <- colnames(Udf)
        U.combined <- rbind(Udf[setdiff(rownames(Udf), rownames(Uproj)), , drop = FALSE], Uproj)
        return(U.combined)
    }
    
    # Multi-view path (named lists)
    stopifnot(is.list(Vdf), is.list(table.rclr.project))
    views <- intersect(names(Vdf), names(table.rclr.project))
    if (!length(views)) stop("[.transform_helper] No overlapping views.")
    
    # Project per view, then sum contributions in the shared latent space
    Usum <- NULL
    ncomp <- ncol(Udf)
    for (vw in views) {
        Vvw <- Vdf[[vw]]
        Tvw <- table.rclr.project[[vw]]
        stopifnot(is.matrix(Vvw), is.matrix(Tvw))
        
        common <- intersect(rownames(Vvw), rownames(Tvw))
        if (length(common) < ncomp) {
            stop(sprintf("[.transform_helper] View '%s': too few matching features (%d).", vw, length(common)))
        }
        
        M <- t(as.matrix(Tvw[common, , drop = FALSE]))     
        V <- as.matrix(Vvw[common, , drop = FALSE])        
        
        # Accumulate per-view U
        Uvw <- M %*% V                                     
        if (is.null(Usum)) {
            Usum <- Uvw
        } else {
            # Align rows (samples) by name before summing
            all_s <- union(rownames(Usum), rownames(Uvw))
            Utmp  <- matrix(0, nrow = length(all_s), ncol = ncol(Udf),
                            dimnames = list(all_s, colnames(Udf)))
            Utmp[rownames(Usum), ] <- Usum
            Utmp[rownames(Uvw), ]  <- Utmp[rownames(Uvw), ] + Uvw
            Usum <- Utmp
        }
    }
    
    # Sample dedup (after combining views)
    if (dedup.samples) {
        sid <- sub("_\\d+$", "", rownames(Usum))
        if (any(duplicated(sid))) {
            Usum <- rowsum(Usum, group = sid, reorder = FALSE) / as.vector(table(sid))
        } else {
            rownames(Usum) <- sid
        }
    }
    
    # Scale by S
    if (length(s.eig)) {
        Sinv <- diag(1 / s.eig, nrow = length(s.eig))
        Usum <- Usum %*% Sinv
    }
    colnames(Usum) <- colnames(Udf)
    
    # Merge with training U, avoiding duplicates
    keep_train <- setdiff(rownames(Udf), rownames(Usum))
    rbind(Udf[keep_train, , drop = FALSE], Usum)
}

#' RPCA Table Filtering and Preprocessing
#'
#' Internal function that performs filtering and cleanup on a compositional data table
#' prior to Robust PCA analysis. Removes low-count features/samples, enforces non-zero frequency thresholds,
#' checks for ID duplication, and returns a matrix suitable for transformation and ordination.
#'
#' @param table A matrix or data frame with features as rows and samples as columns.
#' @param min.sample.count Minimum total count required for a sample to be retained. Default is 0.
#' @param min.feature.count Minimum total count required for a feature to be retained. Default is 0.
#' @param min.feature.frequency Minimum percentage (0–100) of samples in which a feature must be non-zero. Default is 0.
#'
#' @return A filtered numeric matrix containing non-empty features and samples.
#' @keywords internal
#' @noRd

.rpca_table_processing <- function(table,
                                   min.sample.count = 0,
                                   min.feature.count = 0,
                                   min.feature.frequency = 0) {
    # Ensure the input is a matrix
    if (is.data.frame(table)) {
        table <- as.matrix(table)
    }
    
    n.features <- nrow(table)
    n.samples  <- ncol(table)
    
    # Filter features by total count
    if (!is.null(min.feature.count)) {
        feature.totals <- rowSums(table, na.rm = TRUE)
        keep.features <- feature.totals > min.feature.count
        table <- table[keep.features, , drop = FALSE]
    }
    
    # Filter features by frequency across samples
    if (!is.null(min.feature.frequency)) {
        freq.threshold <- min.feature.frequency / 100
        feature.freq <- rowMeans(table > 0, na.rm = TRUE)
        keep.features <- feature.freq > freq.threshold
        table <- table[keep.features, , drop = FALSE]
    }
    
    # Filter samples by total count
    if (!is.null(min.sample.count)) {
        sample.totals <- colSums(table, na.rm = TRUE)
        keep.samples <- sample.totals > min.sample.count
        table <- table[, keep.samples, drop = FALSE]
    }
    
    # Check for duplicate IDs
    if (any(duplicated(colnames(table)))) {
        stop("Data table contains duplicate sample (column) IDs.", call. = FALSE)
    }
    if (any(duplicated(rownames(table)))) {
        stop("Data table contains duplicate feature (row) IDs.", call. = FALSE)
    }
    
    # Remove empty rows and columns if sample filtering applied
    if (!is.null(min.sample.count)) {
        nonzero.features <- rowSums(table, na.rm = TRUE) > 0
        nonzero.samples  <- colSums(table, na.rm = TRUE) > 0
        table <- table[nonzero.features, nonzero.samples, drop = FALSE]
    }
    
    return(table)
}

#' Generate a MaskedMatrix from Numeric Input
#'
#' Internal helper that ensures input is a 2D numeric matrix and returns a masked version,
#' replacing non-finite values (e.g., NA, NaN, Inf) with `NA` and recording their positions in a logical mask.
#'
#' @param mat A numeric matrix or vector. If a vector, it will be converted to a 1-row matrix.
#'
#' @return A list of class \code{"MaskedMatrix"} with two elements:
#' \describe{
#'   \item{data}{The original matrix with non-finite values replaced by \code{NA}.}
#'   \item{mask}{Logical matrix indicating non-finite entries (TRUE if missing).}
#' }
#'
#' @keywords internal
#' @noRd

.mask_value_only <- function(mat) {
    # Ensure matrix is at least 2D
    if (is.vector(mat)) {
        mat <- matrix(mat, nrow = 1)
    }
    
    # Ensure matrix is not more than 2D
    if (length(dim(mat)) > 2) {
        stop("Input matrix can only have two dimensions or less")
    }
    
    # Generate logical mask: TRUE where values are missing
    mask <- !is.finite(mat)  
    
    # Create masked matrix
    masked.mat <- mat
    masked.mat[!is.finite(mat)] <- NA
    
    # Return as a masked matrix
    return(structure(list(
        data = masked.mat,
        mask = mask
    ), class = "MaskedMatrix"))
}

#' Internal constructor for ordination results
#' @keywords internal
#' @noRd
.ordination_results <- function(method, eigvals, samples, features,
                                proportion.explained, dist = NULL, metadata = list()) {
    return(structure(list(
        method = method,
        eigvals = eigvals,
        samples = samples,
        features = features,
        proportion.explained = proportion.explained,
        dist = dist,
        metadata = metadata
    ), class = "OrdinationResults"))
}

#' Internal constructor for a distance matrix object
#' @keywords internal
#' @noRd
.distance_matrix <- function(matrix, ids = NULL, method = "euclidean") {
    if (!is.matrix(matrix)) stop("Input must be a matrix.")
    if (!isSymmetric(matrix)) stop("Distance matrix must be symmetric.")
    if (!is.null(ids)) {
        if (length(ids) != nrow(matrix)) stop("Length of 'ids' must match matrix dimensions.")
        rownames(matrix) <- ids
        colnames(matrix) <- ids
    }
    return(structure(list(
        data = matrix,
        ids = rownames(matrix),
        method = method
    ), class = "DistanceMatrix"))
}

#' Extract per-experiment assay tables from a MultiAssayExperiment
#'
#' @param x A MultiAssayExperiment.
#' @param experiments Character vector of experiment names (or NULL for all).
#'
#' @return A list with `tables`, `experiments`, and `assay_names_used`.
#'
#' @keywords internal
#' @noRd
.extract_mae_tables <- function(x, experiments = NULL) {
    exps <- experiments(x)
    
    if (is.null(experiments)) {
        experiments <- names(exps)
    }
    if (length(experiments) == 0L) {
        stop("No experiments found in 'x'.", call. = FALSE)
    }
    
    assay_names_used <- setNames(character(length(experiments)), experiments)
    tables <- vector("list", length(experiments))
    names(tables) <- experiments
    
    for (i in seq_along(experiments)) {
        e <- experiments[[i]]
        exp_se <- exps[[e]]
        if (is.null(exp_se)) {
            stop(sprintf("Experiment '%s' not found in 'x'.", e), call. = FALSE)
        }
        
        anm <- assayNames(exp_se)
        default_assay <- if (length(anm)) anm[[1]] else 1L
        
        assay_names_used[[e]] <- if (is.character(default_assay)) default_assay else as.character(default_assay)
        tables[[e]] <- assay(exp_se, default_assay)
    }
    
    list(
        tables = tables,
        experiments = experiments,
        assay_names_used = assay_names_used
    )
}