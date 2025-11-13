#' OptSpace back-end (joint & single-view)
#'
#' Internal OptSpace solver + helpers. Prefers \code{vegan::optspace()} when
#' available; otherwise falls back to the internal \code{.optspace}.
#'
#' This file contains:
#' \itemize{
#'   \item \code{.optspace_helper()} — single-view rCLR → OptSpace → biplot
#'   \item \code{.joint_optspace_solve()} — multi-view joint factorization
#'   \item \code{.optspace()} and its auxiliaries — internal fallback solver
#' }
#'
#' @keywords internal
#' @noRd
NULL

#' OptSpace-Based Dimensionality Reduction and RPCA Biplot Generation
#'
#' Internal function that fits an OptSpace model to a rCLR-transformed compositional table,
#' reconstructs the low-rank matrix, applies PCA, and constructs an ordination result
#' capturing sample embeddings, feature loadings, and explained variance. A distance matrix
#' is also generated using Aitchison geometry.
#'
#' @param rclr.table A numeric matrix representing rCLR-transformed compositional data.
#' @param feature.ids Character vector of feature names (used for row labeling of loadings).
#' @param subject.ids Character vector of sample names (used for row labeling of embeddings).
#' @param n.components Integer specifying number of principal components to retain. Default is 3.
#' @param max.iterations Maximum number of iterations to run OptSpace optimization. Default is 5.
#'
#' @return A list with:
#' \describe{
#'   \item{ord.res}{An \code{OrdinationResults} object containing PCA scores, loadings, and metadata.}
#'   \item{dist}{A sample-by-sample \code{DistanceMatrix} object using Aitchison geometry.}
#'   \item{opt.fit}{The raw OptSpace fit result containing matrices \code{X}, \code{Y}, and \code{S}.}
#' }
#'
#' @keywords internal

.optspace_helper <- function(rclr.table,
                             feature.ids,
                             subject.ids,
                             n.components = 3,
                             max.iterations = 5) {
  
  #prefer vegan::optspace() if available; otherwise fall back to internal .optspace
  opt.result <- NULL
  if (requireNamespace("vegan", quietly = TRUE) &&
      "optspace" %in% getNamespaceExports("vegan")) {
    opt.result <- try({
      getExportedValue("vegan", "optspace")(
        rclr.table,
        r      = n.components,       
        niter  = max.iterations,
        tol    = 1e-5,
        center = TRUE,               
        scale  = FALSE
      )
    }, silent = TRUE)
    #ensure expected slots exist; otherwise fall back
    if (inherits(opt.result, "try-error") ||
        is.null(opt.result$X) || is.null(opt.result$S) || is.null(opt.result$Y)) {
      opt.result <- NULL
    }
  }
  if (is.null(opt.result)) {
    opt.result <- .optspace(
      rclr.table,
      ropt    = n.components,
      niter   = max.iterations,
      tol     = 1e-5,
      verbose = FALSE
    )
  }
  
  #update n.components
  n.components <- ncol(opt.result$S)
  
  #reconstruct and re-center matrix
  X.hat <- opt.result$X %*% opt.result$S %*% t(opt.result$Y)
  X.hat <- scale(X.hat, center = TRUE, scale = FALSE)
  X.hat <- t(scale(t(X.hat), center = TRUE, scale = FALSE))
  
  #PCA
  svd.out <- svd(X.hat)
  u <- svd.out$u[, 1:n.components, drop = FALSE]
  s <- svd.out$d[1:n.components]
  v <- svd.out$v[, 1:n.components, drop = FALSE]
  
  #label loadings
  rename.cols <- paste0("PC", seq_len(n.components))
  sample.scores <- data.frame(u, row.names = subject.ids)
  feature.scores <- data.frame(v, row.names = feature.ids)
  colnames(sample.scores) <- rename.cols
  colnames(feature.scores) <- rename.cols
  
  #proportion explained
  prop.var <- s^2 / sum(svd.out$d^2)
  names(prop.var) <- rename.cols
  names(s) <- rename.cols
  
  #add PC3 for 2D case
  if (n.components == 2) {
    sample.scores$PC3 <- 0
    feature.scores$PC3 <- 0
    s <- c(s, PC3 = 0)
    prop.var <- c(prop.var, PC3 = 0)
    rename.cols <- c(rename.cols, "PC3")
  }
  
  #compute distance
  dist.matrix.raw <- as.matrix(dist(u))
  rownames(dist.matrix.raw) <- subject.ids
  colnames(dist.matrix.raw) <- subject.ids
  
  #wrap with DistanceMatrix
  dist.res <- .DistanceMatrix(dist.matrix.raw, ids = subject.ids, method = "aitchison")
  
  #build OrdinationResults object
  ord.res <- .OrdinationResults(
    method = "rpca_biplot",
    eigvals = s,
    samples = sample.scores,
    features = feature.scores,
    proportion.explained = prop.var,
    dist = dist.matrix.raw,  
    metadata = list(
      long.method.name = "(Robust Aitchison) RPCA Biplot",
      run.id = sprintf("optspace_helper_n.components_%d.max.iterations_%d", 
                       n.components, max.iterations)
    )
  )
  
  return(list(
    ord.res = ord.res,
    dist = dist.res,  
    opt.fit = opt.result
  ))
}

#' Joint OptSpace Optimization Across Multiple Train/Test Splits
#'
#' Internal function that performs joint matrix factorization using OptSpace across a set of paired train/test compositional tables.
#' Stacks training matrices horizontally, applies low-rank optimization, splits feature loadings per table,
#' and evaluates projection error on test data via Frobenius norm.
#'
#' @param train.test.pairs A list of paired matrices where each element is a two-item list: \code{[[test, train]]}.
#' @param n.components Integer specifying number of components to retain in the OptSpace model.
#' @param max.iter Maximum number of optimization iterations. Default is 50.
#' @param verbose Logical; whether to print progress messages. Default is \code{TRUE}.
#'
#' @return A list with:
#' \describe{
#'   \item{U}{Shared sample embedding matrix across all input tables.}
#'   \item{S}{Singular values matrix from OptSpace decomposition.}
#'   \item{V.list}{List of per-table feature loading matrices.}
#'   \item{dists}{Matrix of reconstruction errors (rows: error type, columns: tables).}
#' }
#'
#' @keywords internal

.joint_optspace_solve <- function(train.test.pairs, n.components,
                                  max.iter = 50, verbose = TRUE) {
  #prepare lists to hold training matrices and dimensions
  train.matrices <- list()
  test.matrices <- list()
  dims <- list()
  
  for (pair in train.test.pairs) {
    test.mat <- pair[[1]]
    train.mat <- pair[[2]]
    train.matrices <- append(train.matrices, list(train.mat))
    test.matrices <- append(test.matrices, list(test.mat))
    dims <- append(dims, list(dim(train.mat)))
  }
  
  #stack training matrices horizontally 
  train.stacked <- do.call(cbind, train.matrices)
  
  #apply OptSpace to stacked matrix
  if (verbose) message("Running optspace() on stacked training data...")
  #prefer vegan::optspace() if available; otherwise fall back to internal .optspace
  fit <- NULL
  if (requireNamespace("vegan", quietly = TRUE) &&
      "optspace" %in% getNamespaceExports("vegan")) {
    fit <- try({
      getExportedValue("vegan", "optspace")(
        train.stacked,
        r      = n.components,      
        niter  = max.iter,
        tol    = 1e-5,
        center = TRUE,              
        scale  = FALSE
      )
    }, silent = TRUE)
    #if anything odd, drop to internal
    if (inherits(fit, "try-error") || is.null(fit$X) || is.null(fit$S) || is.null(fit$Y)) {
      fit <- NULL
    }
  }
  
  if (is.null(fit)) {
    fit <- .optspace(
      train.stacked,
      ropt  = n.components,
      niter = max.iter,
      tol   = 1e-5,
      verbose = FALSE
    )
  }
  
  #extract sample loadings
  U.shared <- fit$X
  S.shared <- fit$S
  
  #split V back into per-table pieces
  feat.indices <- cumsum(sapply(dims, function(d) d[2]))
  feat.starts <- c(1, head(feat.indices, -1) + 1)
  V.list <- Map(function(start, end) fit$Y[start:end, , drop = FALSE],
                feat.starts, feat.indices)
  
  dists <- matrix(0, nrow = 2, ncol = max.iter)
  for (i in seq_along(test.matrices)) {
    V.k <- V.list[[i]]
    test.mat <- test.matrices[[i]]
    
    #project test samples: U.test = test × V
    U.test <- as.matrix(test.mat) %*% V.k
    U.test <- sweep(U.test, 2, diag(S.shared), "/")
    recon.test <- U.test %*% S.shared %*% t(V.k)
    
    #center for consistency
    recon.test <- scale(recon.test, center = TRUE, scale = FALSE)
    recon.test <- t(scale(t(recon.test), center = TRUE, scale = FALSE))
    
    error <- test.mat - recon.test
    error[is.na(error)] <- 0
    error.val <- norm(error, "F") / sqrt(sum(!is.na(test.mat)))
    
    dists[1, i] <- error.val
    dists[2, i] <- 0  
  }
  
  list(U = U.shared, S = S.shared, V.list = V.list, dists = dists)
}