#FUNCTION FOR JOINT RPCA

#' Joint Robust PCA on Multiple Compositional Tables
#'
#' Performs Joint Robust Principal Component Analysis (RPCA) using OptSpace on multiple compositional tables.
#' Automatically handles shared sample alignment, optional train/test splitting, and compositional preprocessing.
#'
#' @param tables A list of compositional data tables (matrices or data frames).
#'        For automatic handling of Bioconductor objects (e.g., MultiAssayExperiment, 
#'        TreeSummarizedExperiment, SummarizedExperiment), or single matrices, 
#'        consider using `jointRPCAuniversal()` instead.
#' @param n.test.samples Integer specifying the number of samples to hold out for testing (only used if `sample.metadata` is NULL). Default is 10.
#' @param sample.metadata Optional data frame containing sample-level metadata. Must include a column indicating train/test labels.
#' @param train.test.column The name of the column in `sample.metadata` that defines training vs test samples.
#' @param n.components Integer specifying the number of principal components to compute. Must be at least 2.
#' @param rclr.transform.tables Logical; whether to apply rclr transformation to each input table before ordination. Default is TRUE.
#' @param min.sample.count Minimum total count required for a sample to be retained during filtering. Default is 0 (no filtering).
#' @param min.feature.count Minimum total count required for a feature to be retained. Default is 0 (no filtering).
#' @param min.feature.frequency Minimum percentage (0–100) of samples in which a feature must be non-zero to be retained. Default is 0.
#' @param max.iterations Maximum number of optimization iterations. Must be at least 1.
#'
#' @return A list with:
#' \describe{
#'   \item{ord.res}{An \code{OrdinationResults} object containing sample scores, feature loadings, and variance explained.}
#'   \item{U.dist.res}{A distance matrix (samples × samples) based on the learned sample embeddings.}
#'   \item{cv.dist}{A data frame of cross-validation error statistics across iterations.}
#'   \item{rclr.tables}{A list of CLR-transformed tables used as input for the ordination. These are useful for accessing dataset-specific embeddings or for reuse in downstream analyses.}
#' }
#'
#' @examples
#' \dontrun{
#' # See examples/joint_rpca_example.qmd for a full reproducible test case.
#' }
#'
#' @export

jointRPCA <- function(tables,
                       n.test.samples = 10,
                       sample.metadata = NULL,
                       train.test.column = NULL,
                       n.components = 3,
                       rclr.transform.tables = TRUE,
                       min.sample.count = 0,
                       min.feature.count = 0,
                       min.feature.frequency = 0,
                       max.iterations = 5) {
  
    if (n.components < 2) stop("n.components must be at least 2.")
    if (max.iterations < 1) stop("max.iterations must be at least 1.")
  
    #filter each table
    if (rclr.transform.tables) {
        tables <- lapply(tables, function(tbl) {
            .rpca_table_processing(tbl,
                                min.sample.count = min.sample.count,
                                min.feature.count = min.feature.count,
                                min.feature.frequency = min.feature.frequency)
        })
    }
  
    #find shared samples
    sample.sets <- lapply(tables, colnames)
    shared.all.samples <- Reduce(intersect, sample.sets)
    if (length(shared.all.samples) == 0) {
        stop("No samples overlap between all tables. If using pre-transformed tables, set rclr.transform.tables = FALSE.")
    }
    unshared.samples <- setdiff(unique(unlist(sample.sets)), shared.all.samples)
    if (length(unshared.samples) > 0) {
        warning(sprintf("Removing %d sample(s) that do not overlap in tables.", length(unshared.samples)))
    }
  
    #re-filter tables to shared samples
    tables <- lapply(tables, function(tbl) {
        tbl <- tbl[, shared.all.samples, drop = FALSE]
        if (rclr.transform.tables) {
            .rpca_table_processing(tbl,
                                min.sample.count = min.sample.count,
                                min.feature.count = min.feature.count,
                                min.feature.frequency = min.feature.frequency)
        } else {
            tbl
        }
    })
    shared.all.samples <- Reduce(intersect, lapply(tables, colnames))
  
    #transform tables
    rclr.tables <- lapply(tables, function(tbl) {
        mat <- as.matrix(tbl)
        if (rclr.transform.tables) {
            vegan::decostand(mat, method = "rclr", MARGIN = 2)
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
          rclr.table = t(rclr.tables[[1]]),
          feature.ids = rownames(rclr.tables[[1]]),
          subject.ids = colnames(rclr.tables[[1]]),
          n.components = n.components,
          max.iterations = max.iterations
        )$ord.res
        sorted.ids <- rownames(ord.tmp$samples[order(ord.tmp$samples[, 1]), ])
        idx <- round(seq(1, length(sorted.ids), length.out = n.test.samples))
        test.samples <- sorted.ids[idx]
        train.samples <- setdiff(shared.all.samples, test.samples)
    }
  
    #run joint RPCA via helper
    result <- .joint_optspace_helper(
        tables = rclr.tables,
        n.components = n.components,
        max.iterations = max.iterations,
        test.samples = test.samples,
        train.samples = train.samples
    )
  
    return(list(
      ord.res = result$ord.res,
      dist = result$dist,
      cv.stats = result$cv.stats,
      rclr.tables = rclr.tables 
    ))
}
