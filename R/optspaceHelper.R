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
  
    #fit OptSpace
    opt.result <- .optspace(rclr.table, ropt = n.components, niter = max.iterations, tol = 1e-5, verbose = FALSE)
  
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