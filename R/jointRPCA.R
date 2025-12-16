#' Joint RPCA front-end and orchestration
#'
#' User-facing wrappers and internal orchestration of Joint Robust Principal Component Analysis (Joint-RPCA)
#' for compositional multi-omic data using the OptSpace algorithm.
#'
#' This file contains:
#' \itemize{
#'   \item \code{jointRPCAuniversal()} — universal front-end detecting input type
#'         (\code{matrix}, \code{list}, \code{SummarizedExperiment}, or \code{MultiAssayExperiment}).
#'   \item \code{runJointRPCA()} — high-level wrapper that runs Joint-RPCA and
#'         stores the embedding in \code{reducedDim()}.
#'   \item \code{.jointRPCA()} — internal Joint-RPCA engine handling preprocessing,
#'         shared-sample alignment, factorization, and output assembly.
#'   \item \code{.joint_optspace_helper()} — internal coordination of multi-view OptSpace factorization.
#'   \item \code{.transform()} and \code{.transform_helper()} — projection of new data into existing ordination space.
#'   \item \code{.rpca_table_processing()} — compositional preprocessing (filtering, zero-removal, etc.).
#'   \item \code{.mask_value_only()} — safe numeric matrix masking for missing values.
#'   \item Other small utility functions used across the Joint-RPCA workflow.
#' }
#'
#' These functions depend on the OptSpace backend defined in \code{optspace.R}.
#' Only \code{jointRPCAuniversal()} and \code{runJointRPCA()} are exported; all others are internal helpers.
#'
#' @keywords internal
#' @noRd
NULL

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
#' @param ... Additional arguments passed to \code{.jointRPCA()}.
#'
#' @return Output from \code{.jointRPCA()} with extra fields when \code{x} is a
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
    
    res <- .jointRPCA(tables = tables, ...)
    
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
#'   and from there to the internal \code{.jointRPCA()} engine (e.g. \code{n.components},
#'   \code{min.sample.count}, \code{max.iterations}, etc.).
#'
#' @return The input object \code{x} with a new entry in
#'   \code{reducedDim(x, name)} containing the Joint-RPCA sample embedding.
#'   The full Joint-RPCA result (including distances, cross-validation
#'   statistics and transformed tables) is stored in
#'   \code{metadata(x)$JointRPCA[[name]]}.
#'
#' @export

runJointRPCA <- function(x,
                         experiments = NULL,
                         altexp = NULL,
                         name   = "JointRPCA",
                         ...) {
    #select the object to operate on
    y <- x
    if (!is.null(altexp)) {
        y <- altExp(x, altexp)
    }
    
    #use universal front-end to build tables + run .jointRPCA()
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

#' Joint Robust PCA on Multiple Compositional Tables
#'
#' Internal engine for Joint Robust Principal Component Analysis (RPCA) using
#' OptSpace on multiple compositional tables.
#'
#' This function assumes a list of already extracted tables and is typically
#' called via \code{jointRPCAuniversal()} or \code{runJointRPCA()}.
#'
#' @param tables A list of compositional data tables (matrices or data frames).
#' @param n.test.samples Integer specifying the number of samples to hold out for testing
#'   (only used if \code{sample.metadata} is \code{NULL}). Default is 10.
#' @param sample.metadata Optional data frame containing sample-level metadata.
#' @param train.test.column The name of the column in \code{sample.metadata}
#'   that defines training vs test samples.
#' @param n.components Integer specifying the number of principal components to compute.
#' @param rclr.transform.tables Logical; whether to apply rCLR transformation to each
#'   input table before ordination. Default is \code{TRUE}.
#' @param min.sample.count Minimum total count required for a sample to be retained.
#' @param min.feature.count Minimum total count required for a feature to be retained.
#' @param min.feature.frequency Minimum percentage (0–100) of samples in which a
#'   feature must be non-zero to be retained.
#' @param max.iterations Maximum number of optimization iterations.
#'
#' @return A list with \code{ord.res}, \code{dist}, \code{cv.stats}, and
#'   \code{rclr.tables}.
#'
#' @keywords internal
#' @noRd

.jointRPCA <- function(tables,
                      n.test.samples = 10,
                      sample.metadata = NULL,
                      train.test.column = NULL,
                      n.components = 3,
                      rclr.transform.tables = TRUE,
                      min.sample.count = 0,
                      min.feature.count = 0,
                      min.feature.frequency = 0,
                      max.iterations = 5) {
    
    if (is.null(names(tables))) {
        names(tables) <- paste0("view", seq_along(tables))
    }
    
    if (n.components < 2) stop("n.components must be at least 2.")
    if (max.iterations < 1) stop("max.iterations must be at least 1.")
    
    #filtering (always done, independent of rclr)
    tables <- lapply(tables, function(tbl) {
        .rpca_table_processing(
            tbl,
            min.sample.count      = min.sample.count,
            min.feature.count     = min.feature.count,
            min.feature.frequency = min.feature.frequency
        )
    })
    
    #find shared samples across views
    sample.sets <- lapply(tables, colnames)
    shared.all.samples <- Reduce(intersect, sample.sets)
    if (length(shared.all.samples) == 0) {
        stop("No samples overlap between all tables. If using pre-transformed tables, set rclr.transform.tables = FALSE.")
    }
    unshared.samples <- setdiff(unique(unlist(sample.sets)), shared.all.samples)
    if (length(unshared.samples) > 0) {
        warning(sprintf("Removing %d sample(s) that do not overlap in tables.", length(unshared.samples)))
    }
    
    #restrict each table to the shared sample set
    tables <- lapply(tables, function(tbl) {
        tbl[, shared.all.samples, drop = FALSE]
    })
    shared.all.samples <- Reduce(intersect, lapply(tables, colnames))
    
    #transform tables: rCLR or masking
    rclr.tables <- lapply(tables, function(tbl) {
        mat <- as.matrix(tbl)
        rown <- rownames(mat)
        coln <- colnames(mat)
        
        if (rclr.transform.tables) {
            
            mat[!is.finite(mat)] <- 0
            mat[mat < 0]         <- 0
            
            out <- vegan::decostand(mat, method = "rclr", MARGIN = 2)
            
            dimnames(out) <- list(rown, coln)
            out
        } else {
            .mask_value_only(mat)$data
        }
    })
    names(rclr.tables) <- names(tables)
    
    #determine train/test split
    if (!is.null(sample.metadata) && !is.null(train.test.column)) {
        md <- as.data.frame(sample.metadata)
        md <- md[shared.all.samples, , drop = FALSE]
        train.samples <- rownames(md)[md[[train.test.column]] == "train"]
        test.samples  <- rownames(md)[md[[train.test.column]] == "test"]
    } else {
        ord.tmp <- .optspace_helper(
            rclr.table   = t(rclr.tables[[1]]),
            feature.ids  = rownames(rclr.tables[[1]]),
            subject.ids  = colnames(rclr.tables[[1]]),
            n.components = n.components,
            max.iterations = max.iterations
        )$ord.res
        sorted.ids <- rownames(ord.tmp$samples[order(ord.tmp$samples[, 1]), ])
        idx <- round(seq(1, length(sorted.ids), length.out = n.test.samples))
        test.samples <- sorted.ids[idx]
        train.samples <- setdiff(shared.all.samples, test.samples)
    }
    
    #run joint OptSpace
    result <- .joint_optspace_helper(
        tables        = rclr.tables,
        n.components  = n.components,
        max.iterations = max.iterations,
        test.samples  = test.samples,
        train.samples = train.samples,
        sample.order  = shared.all.samples
    )
    
    list(
        ord.res      = result$ord.res,
        dist         = result$dist,
        cv.stats     = result$cv.stats,
        rclr.tables  = rclr.tables
    )
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
#'   \item{ord.res}{An \code{OrdinationResults} object containing embeddings, loadings, and variance explained.}
#'   \item{dist}{A \code{DistanceMatrix} object for sample embeddings.}
#'   \item{cv.stats}{A data frame summarizing reconstruction error across iterations and tables.}
#' }
#'
#' @keywords internal

.joint_optspace_helper <- function(tables,
                                   n.components,
                                   max.iterations,
                                   test.samples,
                                   train.samples,
                                   sample.order = NULL) {
    
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
        n.components      = n.components,
        max.iter          = max.iterations
    )
    
    U <- opt.result$U
    S <- opt.result$S
    V.list <- opt.result$V.list
    dists <- opt.result$dists
    
    #assign row/column names to loadings
    pc.names <- paste0("PC", seq_len(n.components))
    
    #combine feature loadings with table-derived row names
    vjoint <- do.call(rbind, Map(function(tbl, V) {
        rownames(V) <- rownames(tbl)
        colnames(V) <- pc.names
        V
    }, tables, V.list))
    
    U <- U[seq_along(train.samples), , drop = FALSE]
    rownames(U) <- train.samples
    colnames(U) <- pc.names
    
    #recenter & re-factor via SVD
    X <- U %*% S %*% t(vjoint)
    X <- sweep(X, 2, colMeans(X))
    X <- sweep(X, 1, rowMeans(X))
    svd.res <- svd(X)
    u <- svd.res$u[, seq_len(n.components), drop = FALSE]
    v <- svd.res$v[, seq_len(n.components), drop = FALSE]
    s.eig <- svd.res$d[seq_len(n.components)]
    
    rownames(u) <- train.samples
    rownames(v) <- rownames(vjoint)
    pc.names <- paste0("PC", seq_len(n.components))
    colnames(u) <- colnames(v) <- pc.names
    
    #build a named per-view features list
    features_list <- lapply(seq_along(tables), function(i) {
        rid <- rownames(tables[[i]])          
        v[rid, , drop = FALSE]                
    })
    names(features_list) <- names(tables)    
    
    prop.exp <- s.eig^2 / sum(s.eig^2)
    ord.res <- .OrdinationResults(
        method = "rpca",
        eigvals = setNames(s.eig, pc.names),
        samples = u,
        features = features_list,              
        proportion.explained = setNames(prop.exp, pc.names)
    )
    
    #project test samples
    if (length(test.samples) > 0) {
        test.matrices <- lapply(tables, function(tbl) tbl[, test.samples, drop = FALSE])
        names(test.matrices) <- names(tables)  
        ord.res <- .transform(ord.res, test.matrices, apply.rclr = FALSE)
    }
    
    #compute distance matrix and CV error summary
    dist.base <- as.matrix(dist(ord.res$samples))
    
    if (!is.null(sample.order)) {
        order_use <- intersect(sample.order, rownames(dist.base))
        dist.mat  <- dist.base[order_use, order_use, drop = FALSE]
    } else {
        dist.mat  <- dist.base
        order_use <- rownames(dist.base)
    }
    
    dist.res <- .DistanceMatrix(dist.mat, ids = order_use)
    
    cv.dist <- data.frame(t(dists))
    colnames(cv.dist) <- c("mean_CV", "std_CV")
    cv.dist$run <- sprintf("tables_%d.n.components_%d.max.iterations_%d.n.test_%d",
                           length(tables), n.components, max.iterations, length(test.samples))
    cv.dist$iteration <- seq_len(nrow(cv.dist))
    rownames(cv.dist) <- seq_len(nrow(cv.dist))
    
    list(ord.res = ord.res, dist = dist.res, cv.stats = cv.dist)
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

.transform <- function(ordination, tables,
                       apply.rclr = TRUE) {
    
    Udf    <- ordination$samples
    Vobj   <- ordination$features
    s.eig  <- ordination$eigvals
    
    #ensure tables is a list of views
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
    
    #rCLR if requested 
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
    
    #3+-omic path: V is a named list per view
    if (!is.list(Vobj) || is.null(names(Vobj))) {
        stop("[.transform] ordination$features is neither a matrix nor a named list.")
    }
    
    #intersect views by name, preserve training order
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
    ordination
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

.transform_helper <- function(Udf, Vdf, s.eig, table.rclr.project,
                              dedup.samples = TRUE) {
    
    #legacy path (single view)
    if (is.matrix(Vdf)) {
        stopifnot(is.matrix(table.rclr.project))
        #align rows by name
        common <- intersect(rownames(Vdf), rownames(table.rclr.project))
        if (length(common) < ncol(Udf))
            stop(sprintf("[.transform_helper] Too few matching features: %d", length(common)))
        
        M <- t(as.matrix(table.rclr.project[common, , drop = FALSE]))   
        V <- as.matrix(Vdf[common, , drop = FALSE])                     
        
        #dedup of sample IDs
        if (dedup.samples) {
            sid <- sub("_\\d+$", "", rownames(M))
            if (any(duplicated(sid))) {
                M <- rowsum(M, group = sid, reorder = FALSE) / as.vector(table(sid))
            } else {
                rownames(M) <- sid
            }
        }
        
        #projection (match training scaling)
        Uproj <- M %*% V
        #scale by singular values
        if (length(s.eig)) {
            Sinv <- diag(1 / s.eig, nrow = length(s.eig))
            Uproj <- Uproj %*% Sinv
        }
        
        colnames(Uproj) <- colnames(Udf)
        U.combined <- rbind(Udf[setdiff(rownames(Udf), rownames(Uproj)), , drop = FALSE], Uproj)
        return(U.combined)
    }
    
    #multi-view path (named lists)
    stopifnot(is.list(Vdf), is.list(table.rclr.project))
    views <- intersect(names(Vdf), names(table.rclr.project))
    if (!length(views)) stop("[.transform_helper] No overlapping views.")
    
    #project per view, then sum contributions in the shared latent space
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
        
        #accumulate per-view U
        Uvw <- M %*% V                                     
        if (is.null(Usum)) {
            Usum <- Uvw
        } else {
            #align rows (samples) by name before summing
            all_s <- union(rownames(Usum), rownames(Uvw))
            Utmp  <- matrix(0, nrow = length(all_s), ncol = ncol(Udf),
                            dimnames = list(all_s, colnames(Udf)))
            Utmp[rownames(Usum), ] <- Usum
            Utmp[rownames(Uvw), ]  <- Utmp[rownames(Uvw), ] + Uvw
            Usum <- Utmp
        }
    }
    
    #sample dedup (after combining views)
    if (dedup.samples) {
        sid <- sub("_\\d+$", "", rownames(Usum))
        if (any(duplicated(sid))) {
            Usum <- rowsum(Usum, group = sid, reorder = FALSE) / as.vector(table(sid))
        } else {
            rownames(Usum) <- sid
        }
    }
    
    #scale by S
    if (length(s.eig)) {
        Sinv <- diag(1 / s.eig, nrow = length(s.eig))
        Usum <- Usum %*% Sinv
    }
    colnames(Usum) <- colnames(Udf)
    
    #merge with training U, avoiding duplicates
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

.rpca_table_processing <- function(table,
                                   min.sample.count = 0,
                                   min.feature.count = 0,
                                   min.feature.frequency = 0) {
    #ensure the input is a matrix
    if (is.data.frame(table)) {
        table <- as.matrix(table)
    }
    
    n.features <- nrow(table)
    n.samples  <- ncol(table)
    
    #filter features by total count
    if (!is.null(min.feature.count)) {
        feature.totals <- rowSums(table, na.rm = TRUE)
        keep.features <- feature.totals > min.feature.count
        table <- table[keep.features, , drop = FALSE]
    }
    
    #filter features by frequency across samples
    if (!is.null(min.feature.frequency)) {
        freq.threshold <- min.feature.frequency / 100
        feature.freq <- rowMeans(table > 0, na.rm = TRUE)
        keep.features <- feature.freq > freq.threshold
        table <- table[keep.features, , drop = FALSE]
    }
    
    #filter samples by total count
    if (!is.null(min.sample.count)) {
        sample.totals <- colSums(table, na.rm = TRUE)
        keep.samples <- sample.totals > min.sample.count
        table <- table[, keep.samples, drop = FALSE]
    }
    
    #check for duplicate IDs
    if (any(duplicated(colnames(table)))) {
        stop("Data table contains duplicate sample (column) IDs.")
    }
    if (any(duplicated(rownames(table)))) {
        stop("Data table contains duplicate feature (row) IDs.")
    }
    
    #remove empty rows and columns if sample filtering applied
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

.mask_value_only <- function(mat) {
    #ensure matrix is at least 2D
    if (is.vector(mat)) {
        mat <- matrix(mat, nrow = 1)
    }
    
    #ensure matrix is not more than 2D
    if (length(dim(mat)) > 2) {
        stop("Input matrix can only have two dimensions or less")
    }
    
    #generate logical mask: TRUE where values are missing
    mask <- !is.finite(mat)  
    
    #create masked matrix
    masked.mat <- mat
    masked.mat[!is.finite(mat)] <- NA
    
    #return as a masked matrix
    return(structure(list(
        data = masked.mat,
        mask = mask
    ), class = "MaskedMatrix"))
}