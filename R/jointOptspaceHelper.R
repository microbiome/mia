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
                                   train.samples) {
    #split and transpose training/test data per table
    tables.split <- lapply(tables, function(tbl) {
        list(t(tbl[, test.samples, drop = FALSE]),
             t(tbl[, train.samples, drop = FALSE]))
    })
  
    #format input for solver
    tables.for.solver <- lapply(tables.split, function(pair) {
        lapply(pair, as.matrix)
    })
  
    #run joint OptSpace solver
    opt.result <- .joint_optspace_solve(tables.for.solver,
                                        n.components = n.components,
                                        max.iter = max.iterations)
      
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
    dist.mat <- as.matrix(dist(ord.res$samples))
    dist.res <- .DistanceMatrix(dist.mat, ids = rownames(ord.res$samples))
  
    cv.dist <- data.frame(t(dists))
    colnames(cv.dist) <- c("mean_CV", "std_CV")
    cv.dist$run <- sprintf("tables_%d.n.components_%d.max.iterations_%d.n.test_%d",
                           length(tables), n.components, max.iterations, length(test.samples))
    cv.dist$iteration <- seq_len(nrow(cv.dist))
    rownames(cv.dist) <- seq_len(nrow(cv.dist))
  
    list(ord.res = ord.res, dist = dist.res, cv.stats = cv.dist)
}