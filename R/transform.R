#' Apply Projection of New Compositional Tables to Existing Ordination
#'
#' Internal function that transforms and projects new sample tables into an existing Joint RPCA ordination space.
#' It handles feature alignment, optional rCLR preprocessing, padding of missing features,
#' and merges projected samples with existing ones.
#'
#' @param ordination A list containing previous ordination results: `samples`, `features`, and `eigvals`.
#' @param tables A named list of new compositional tables (matrices or data frames)
#'   with features as rows and samples as columns.
#' @param subset.tables Logical; if `TRUE`, removes unshared features in new tables not present in training.
#'   If `FALSE`, stops with an error if features mismatch.
#' @param apply.rclr Logical; whether to apply rCLR transformation to new input tables before projection. Default is `TRUE`.
#'
#' @return An updated ordination list with the `samples` matrix extended to include new projected samples.
#' @keywords internal

.transform <- function(ordination, tables,
                       subset.tables = TRUE,
                       apply.rclr = TRUE) {
  
    Udf    <- ordination$samples                
    Vobj   <- ordination$features               
    s.eig  <- ordination$eigvals
  
    #ensure tables is a named list of views
    if (!is.list(tables)) {
        stop("[.transform] 'tables' must be a list of view matrices (features x samples).")
    }
  
    if (is.null(names(tables)) && is.list(Vobj) && !is.null(names(Vobj))) {
        names(tables) <- names(Vobj)[seq_along(tables)]
    }
  
    if (is.null(names(tables))) {
        stop("[.transform] 'tables' must be a *named* list of view matrices (features x samples).")
    }
  
    #rCLR if requested (must match training)
    prep_view <- function(tab) {
        mat <- as.matrix(tab)
        if (apply.rclr) {
            mat <- vegan::decostand(mat, method = "rclr")
        }
        storage.mode(mat) <- "double"
        mat[!is.finite(mat)] <- 0
        mat
    }
    tables <- lapply(tables, prep_view)
  
    if (is.matrix(Vobj)) {
        all.features <- rownames(Vobj)
        #pad/reorder to training features
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
  
    #3-omic path: V is a named list per view
    if (!is.list(Vobj) || is.null(names(Vobj))) {
        stop("[.transform] ordination$features is neither a matrix nor a named list.")
    }
  
    #intersect views by name, preserve training order
    views <- intersect(names(Vobj), names(tables))
    if (!length(views)) stop("[.transform] No overlapping view names between ordination and new tables.")
  
    #build per-view test matrices aligned to that view's training features
    test.matrices <- list()
    for (vw in views) {
        Vvw <- Vobj[[vw]]
        stopifnot(is.matrix(Vvw), !is.null(rownames(Vvw)))
        mat <- tables[[vw]]
    
        #pad + reorder to training features of this view
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
  
    #project with view-aware helper
    ordination$samples <- .transform_helper(Udf, Vobj, s.eig, test.matrices)
    ordination
}