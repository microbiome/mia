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