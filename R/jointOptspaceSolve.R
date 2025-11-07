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
    fit <- .optspace(train.stacked, ropt = n.components, niter = max.iter, tol = 1e-5, verbose = FALSE)
  
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