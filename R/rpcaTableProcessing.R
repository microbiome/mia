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