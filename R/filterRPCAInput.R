#' Filter raw count assays before RPCA or Joint-RPCA
#'
#' @description
#' `filterRPCAInput()` applies Gemelli-style preprocessing filters to raw count
#' assays before robust centered log-ratio transformation and RPCA.
#'
#' The filters are intended for raw count data, not transformed assays such as
#' `rclr`. Samples are filtered by total count. Features are filtered by total
#' count and prevalence across samples.
#'
#' @param x A `SummarizedExperiment` or `MultiAssayExperiment` object.
#' @param assay.type Character scalar. Assay name to filter when `x` is a
#'   `SummarizedExperiment`.
#' @param experiments Character or integer vector. Experiments to filter when
#'   `x` is a `MultiAssayExperiment`.
#' @param assay.types Character vector. Assay names corresponding to
#'   `experiments` when `x` is a `MultiAssayExperiment`.
#' @param min.sample.count Numeric scalar or `NULL`. Samples are retained when
#'   their total count is strictly greater than this value. If `NULL`, this
#'   filter is skipped.
#' @param min.feature.count Numeric scalar or `NULL`. Features are retained when
#'   their total count is strictly greater than this value. If `NULL`, this
#'   filter is skipped.
#' @param min.feature.frequency Numeric scalar or `NULL`. Features are retained
#'   when their percentage prevalence across samples is strictly greater than
#'   this value. If `NULL`, this filter is skipped.
#'
#' @return The filtered object.
#'
#' @details
#' This helper mirrors Gemelli-style RPCA table preprocessing:
#'
#' * sample total count must be `> min.sample.count`
#' * feature total count must be `> min.feature.count`
#' * feature prevalence percentage must be `> min.feature.frequency`
#'
#' The strict `>` comparison is intentional.
#'
#' For exact comparison with the current Gemelli script, use:
#'
#' ```
#' min.sample.count = 0
#' min.feature.count = 0
#' min.feature.frequency = 0
#' ```
#'
#' @export
filterRPCAInput <- function(
    x,
    assay.type = "counts",
    experiments = NULL,
    assay.types = NULL,
    min.sample.count = 0,
    min.feature.count = 0,
    min.feature.frequency = 0) {
  UseMethod("filterRPCAInput")
}


#' @export
filterRPCAInput.SummarizedExperiment <- function(
    x,
    assay.type = "counts",
    experiments = NULL,
    assay.types = NULL,
    min.sample.count = 0,
    min.feature.count = 0,
    min.feature.frequency = 0) {
  .validate_rpca_filter_args(
    min.sample.count = min.sample.count,
    min.feature.count = min.feature.count,
    min.feature.frequency = min.feature.frequency
  )
  
  if (!assay.type %in% SummarizedExperiment::assayNames(x)) {
    stop("'", assay.type, "' is not an assay in 'x'.", call. = FALSE)
  }
  
  mat <- SummarizedExperiment::assay(x, assay.type)
  
  keep <- .get_rpca_filter_indices(
    mat = mat,
    min.sample.count = min.sample.count,
    min.feature.count = min.feature.count,
    min.feature.frequency = min.feature.frequency
  )
  
  x <- x[keep$features, keep$samples, drop = FALSE]
  
  return(x)
}


#' @export
filterRPCAInput.MultiAssayExperiment <- function(
    x,
    assay.type = "counts",
    experiments = NULL,
    assay.types = NULL,
    min.sample.count = 0,
    min.feature.count = 0,
    min.feature.frequency = 0) {
  .validate_rpca_filter_args(
    min.sample.count = min.sample.count,
    min.feature.count = min.feature.count,
    min.feature.frequency = min.feature.frequency
  )
  
  if (is.null(experiments) || is.null(assay.types)) {
    stop(
      "'experiments' and 'assay.types' must be specified for ",
      "MultiAssayExperiment input.",
      call. = FALSE
    )
  }
  
  if (length(experiments) != length(assay.types)) {
    stop(
      "The lengths of 'experiments' and 'assay.types' must match.",
      call. = FALSE
    )
  }
  
  experiment_names <- names(MultiAssayExperiment::experiments(x))
  
  if (is.numeric(experiments)) {
    if (any(experiments < 1L) ||
        any(experiments > length(experiment_names))) {
      stop(
        "'experiments' contains invalid experiment indices.",
        call. = FALSE
      )
    }
    
    experiments <- experiment_names[experiments]
  }
  
  missing_experiments <- setdiff(experiments, experiment_names)
  
  if (length(missing_experiments) > 0L) {
    stop(
      "Experiment(s) not found: ",
      paste(missing_experiments, collapse = ", "),
      call. = FALSE
    )
  }
  
  # First pass:
  # Filter each raw count table independently.
  for (i in seq_along(experiments)) {
    exp_name <- experiments[[i]]
    assay_name <- assay.types[[i]]
    
    se <- x[[exp_name]]
    
    if (!assay_name %in% SummarizedExperiment::assayNames(se)) {
      stop(
        "'", assay_name, "' is not an assay in experiment '",
        exp_name, "'.",
        call. = FALSE
      )
    }
    
    mat <- SummarizedExperiment::assay(se, assay_name)
    
    keep <- .get_rpca_filter_indices(
      mat = mat,
      min.sample.count = min.sample.count,
      min.feature.count = min.feature.count,
      min.feature.frequency = min.feature.frequency
    )
    
    x[[exp_name]] <- se[keep$features, keep$samples, drop = FALSE]
  }
  
  # Gemelli Joint-RPCA keeps only samples shared across all selected tables.
  shared_samples <- Reduce(
    intersect,
    lapply(experiments, function(exp_name) {
      colnames(x[[exp_name]])
    })
  )
  
  if (length(shared_samples) == 0L) {
    stop(
      "No samples overlap between all selected experiments.",
      call. = FALSE
    )
  }
  
  for (exp_name in experiments) {
    se <- x[[exp_name]]
    x[[exp_name]] <- se[, shared_samples, drop = FALSE]
  }
  
  # Second pass:
  # After shared-sample subsetting, filter features again. This mirrors
  # Gemelli's joint_rpca() structure.
  for (i in seq_along(experiments)) {
    exp_name <- experiments[[i]]
    assay_name <- assay.types[[i]]
    
    se <- x[[exp_name]]
    mat <- SummarizedExperiment::assay(se, assay_name)
    
    keep <- .get_rpca_filter_indices(
      mat = mat,
      min.sample.count = min.sample.count,
      min.feature.count = min.feature.count,
      min.feature.frequency = min.feature.frequency
    )
    
    x[[exp_name]] <- se[keep$features, keep$samples, drop = FALSE]
  }
  
  return(x)
}


.validate_rpca_filter_args <- function(
    min.sample.count,
    min.feature.count,
    min.feature.frequency) {
  
  args <- list(
    min.sample.count = min.sample.count,
    min.feature.count = min.feature.count,
    min.feature.frequency = min.feature.frequency
  )
  
  for (arg_name in names(args)) {
    value <- args[[arg_name]]
    
    if (is.null(value)) {
      next
    }
    
    if (!(is.numeric(value) &&
          length(value) == 1L &&
          !is.na(value) &&
          is.finite(value))) {
      stop(
        "'", arg_name,
        "' must be a single finite numeric value or NULL.",
        call. = FALSE
      )
    }
    
    if (value < 0) {
      stop(
        "'", arg_name,
        "' must be non-negative.",
        call. = FALSE
      )
    }
  }
  
  if (!is.null(min.feature.frequency) &&
      min.feature.frequency > 100) {
    stop(
      "'min.feature.frequency' must be between 0 and 100.",
      call. = FALSE
    )
  }
}


.get_rpca_filter_indices <- function(
    mat,
    min.sample.count = 0,
    min.feature.count = 0,
    min.feature.frequency = 0) {
  
  if (!(is.matrix(mat) || is.data.frame(mat))) {
    stop("'mat' must be a matrix or data.frame.", call. = FALSE)
  }
  
  mat <- as.matrix(mat)
  
  if (!is.numeric(mat)) {
    stop("'mat' must contain numeric raw counts.", call. = FALSE)
  }
  
  if (anyNA(mat)) {
    stop("Missing values are not allowed in raw count data.", call. = FALSE)
  }
  
  if (any(!is.finite(mat))) {
    stop("Non-finite values are not allowed.", call. = FALSE)
  }
  
  if (any(mat < 0)) {
    stop("Negative values are not allowed in raw count data.", call. = FALSE)
  }
  
  if (is.null(rownames(mat))) {
    stop("Feature IDs are required as row names.", call. = FALSE)
  }
  
  if (is.null(colnames(mat))) {
    stop("Sample IDs are required as column names.", call. = FALSE)
  }
  
  if (anyDuplicated(rownames(mat))) {
    stop("Data table contains duplicate feature IDs.", call. = FALSE)
  }
  
  if (anyDuplicated(colnames(mat))) {
    stop("Data table contains duplicate sample IDs.", call. = FALSE)
  }
  
  original_features <- rownames(mat)
  original_samples <- colnames(mat)
  
  # Gemelli calculates the frequency denominator before filtering.
  frequency_denominator <- ncol(mat)
  
  # 1. Feature-count filtering
  if (!is.null(min.feature.count)) {
    keep <- rowSums(mat) > min.feature.count
    mat <- mat[keep, , drop = FALSE]
  }
  
  if (nrow(mat) == 0L) {
    stop("No features remain after feature-count filtering.", call. = FALSE)
  }
  
  # 2. Feature-frequency filtering
  if (!is.null(min.feature.frequency)) {
    frequency <- rowSums(mat > 0) / frequency_denominator * 100
    keep <- frequency > min.feature.frequency
    mat <- mat[keep, , drop = FALSE]
  }
  
  if (nrow(mat) == 0L) {
    stop(
      "No features remain after feature-frequency filtering.",
      call. = FALSE
    )
  }
  
  # 3. Sample filtering is calculated from the already-filtered features.
  if (!is.null(min.sample.count)) {
    keep <- colSums(mat) > min.sample.count
    mat <- mat[, keep, drop = FALSE]
    
    if (ncol(mat) == 0L) {
      stop("No samples remain after sample filtering.", call. = FALSE)
    }
    
    # Reproduce Gemelli's remove_empty() behaviour.
    mat <- mat[rowSums(mat) > 0, , drop = FALSE]
    mat <- mat[, colSums(mat) > 0, drop = FALSE]
  }
  
  if (nrow(mat) == 0L) {
    stop("No features remain after RPCA filtering.", call. = FALSE)
  }
  
  if (ncol(mat) == 0L) {
    stop("No samples remain after RPCA filtering.", call. = FALSE)
  }
  
  list(
    features = original_features %in% rownames(mat),
    samples = original_samples %in% colnames(mat)
  )
}