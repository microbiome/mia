context("RPCA input filtering")


# Dataset setup -----------------------------------------------------------

data("ibdmdb", package = "mia")

filter_mae <- ibdmdb

selected_experiments <- c(1L, 2L)
raw_assay_types <- c("mgx", "mtx")

selected_experiment_names <- names(
  MultiAssayExperiment::experiments(filter_mae)
)[selected_experiments]


# SummarizedExperiment filtering -----------------------------------------

test_that("sample-count filtering uses a strict comparison", {
  se <- filter_mae[[1L]]
  
  mat <- SummarizedExperiment::assay(
    se,
    raw_assay_types[[1L]]
  )
  
  sample_totals <- colSums(mat, na.rm = TRUE)
  threshold <- min(sample_totals)
  
  expect_true(any(sample_totals == threshold))
  expect_true(any(sample_totals > threshold))
  
  expected_samples <- colnames(mat)[
    sample_totals > threshold
  ]
  
  filtered <- filterRPCAInput(
    se,
    assay.type = raw_assay_types[[1L]],
    min.sample.count = threshold,
    min.feature.count = NULL,
    min.feature.frequency = NULL
  )
  
  expect_identical(
    colnames(filtered),
    expected_samples
  )
})


test_that("feature-count filtering uses a strict comparison", {
  se <- filter_mae[[1L]]
  
  mat <- SummarizedExperiment::assay(
    se,
    raw_assay_types[[1L]]
  )
  
  feature_totals <- rowSums(mat, na.rm = TRUE)
  threshold <- min(feature_totals)
  
  expect_true(any(feature_totals == threshold))
  expect_true(any(feature_totals > threshold))
  
  expected_features <- rownames(mat)[
    feature_totals > threshold
  ]
  
  filtered <- filterRPCAInput(
    se,
    assay.type = raw_assay_types[[1L]],
    min.sample.count = NULL,
    min.feature.count = threshold,
    min.feature.frequency = NULL
  )
  
  expect_identical(
    rownames(filtered),
    expected_features
  )
})


test_that("feature-frequency filtering uses percentage prevalence", {
  se <- filter_mae[[1L]]
  
  mat <- SummarizedExperiment::assay(
    se,
    raw_assay_types[[1L]]
  )
  
  feature_frequency <- (
    rowSums(mat > 0, na.rm = TRUE) / ncol(mat)
  ) * 100
  
  threshold <- min(feature_frequency)
  
  expect_true(any(feature_frequency == threshold))
  expect_true(any(feature_frequency > threshold))
  
  expected_features <- rownames(mat)[
    feature_frequency > threshold
  ]
  
  filtered <- filterRPCAInput(
    se,
    assay.type = raw_assay_types[[1L]],
    min.sample.count = NULL,
    min.feature.count = NULL,
    min.feature.frequency = threshold
  )
  
  expect_identical(
    rownames(filtered),
    expected_features
  )
})

test_that("sample counts are calculated after feature filtering", {
  mat <- rbind(
    retained_feature = c(s1 = 0, s2 = 6, s3 = 6),
    removed_feature = c(s1 = 100, s2 = 0, s3 = 0)
  )
  
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = mat)
  )
  
  filtered <- filterRPCAInput(
    se,
    assay.type = "counts",
    min.sample.count = 0,
    min.feature.count = NULL,
    min.feature.frequency = 50
  )
  
  # removed_feature occurs in only 1/3 samples and is removed.
  expect_identical(
    rownames(filtered),
    "retained_feature"
  )
  
  # s1 becomes empty after removed_feature is discarded.
  expect_identical(
    colnames(filtered),
    c("s2", "s3")
  )
})

test_that("NULL disables filtering for SummarizedExperiment input", {
  se <- filter_mae[[1L]]
  
  filtered <- filterRPCAInput(
    se,
    assay.type = raw_assay_types[[1L]],
    min.sample.count = NULL,
    min.feature.count = NULL,
    min.feature.frequency = NULL
  )
  
  expect_identical(
    rownames(filtered),
    rownames(se)
  )
  
  expect_identical(
    colnames(filtered),
    colnames(se)
  )
  
  expect_identical(
    SummarizedExperiment::assay(
      filtered,
      raw_assay_types[[1L]]
    ),
    SummarizedExperiment::assay(
      se,
      raw_assay_types[[1L]]
    )
  )
})



test_that("invalid raw counts are rejected", {
  negative <- matrix(
    c(1, -1, 2, 3),
    nrow = 2,
    dimnames = list(
      c("f1", "f2"),
      c("s1", "s2")
    )
  )
  
  missing <- matrix(
    c(1, NA_real_, 2, 3),
    nrow = 2,
    dimnames = list(
      c("f1", "f2"),
      c("s1", "s2")
    )
  )
  
  negative_se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = negative)
  )
  
  missing_se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = missing)
  )
  
  expect_error(
    filterRPCAInput(
      negative_se,
      assay.type = "counts"
    ),
    "Negative values"
  )
  
  expect_error(
    filterRPCAInput(
      missing_se,
      assay.type = "counts"
    ),
    "Missing values"
  )
})

# MultiAssayExperiment filtering -----------------------------------------

test_that("filterRPCAInput applies Gemelli-compatible zero thresholds", {
  filtered <- filterRPCAInput(
    filter_mae,
    experiments = selected_experiments,
    assay.types = raw_assay_types,
    min.sample.count = 0,
    min.feature.count = 0,
    min.feature.frequency = 0
  )
  
  expect_s4_class(
    filtered,
    "MultiAssayExperiment"
  )
  
  expect_identical(
    colnames(filtered[[selected_experiments[[1L]]]]),
    colnames(filtered[[selected_experiments[[2L]]]])
  )
  
  expect_true(
    length(
      colnames(filtered[[selected_experiments[[1L]]]])
    ) > 0L
  )
  
  for (i in seq_along(selected_experiments)) {
    original_se <- filter_mae[[selected_experiments[[i]]]]
    filtered_se <- filtered[[selected_experiments[[i]]]]
    
    expect_lte(
      nrow(filtered_se),
      nrow(original_se)
    )
    
    expect_lte(
      ncol(filtered_se),
      ncol(original_se)
    )
    
    filtered_mat <- SummarizedExperiment::assay(
      filtered_se,
      raw_assay_types[[i]]
    )
    
    expect_true(
      all(
        colSums(
          filtered_mat,
          na.rm = TRUE
        ) > 0
      )
    )
  }
})


test_that("integer and character experiment selections agree", {
  filtered_by_index <- filterRPCAInput(
    filter_mae,
    experiments = selected_experiments,
    assay.types = raw_assay_types,
    min.sample.count = 0,
    min.feature.count = 0,
    min.feature.frequency = 0
  )
  
  filtered_by_name <- filterRPCAInput(
    filter_mae,
    experiments = selected_experiment_names,
    assay.types = raw_assay_types,
    min.sample.count = 0,
    min.feature.count = 0,
    min.feature.frequency = 0
  )
  
  for (i in seq_along(selected_experiments)) {
    index_se <- filtered_by_index[[selected_experiments[[i]]]]
    
    name_se <- filtered_by_name[[selected_experiment_names[[i]]]]
    
    expect_identical(
      rownames(index_se),
      rownames(name_se)
    )
    
    expect_identical(
      colnames(index_se),
      colnames(name_se)
    )
  }
})


test_that("MultiAssayExperiment arguments are validated", {
  expect_error(
    filterRPCAInput(
      filter_mae,
      experiments = NULL,
      assay.types = NULL
    ),
    "must be specified"
  )
  
  expect_error(
    filterRPCAInput(
      filter_mae,
      experiments = selected_experiments,
      assay.types = "mgx"
    ),
    "lengths"
  )
  
  expect_error(
    filterRPCAInput(
      filter_mae,
      experiments = 100L,
      assay.types = "mgx"
    ),
    "invalid experiment"
  )
  
  expect_error(
    filterRPCAInput(
      filter_mae,
      experiments = c(
        selected_experiment_names[[1L]],
        "missing"
      ),
      assay.types = raw_assay_types
    ),
    "not found"
  )
  
  expect_error(
    filterRPCAInput(
      filter_mae,
      experiments = selected_experiments,
      assay.types = c("mgx", "missing")
    ),
    "not an assay"
  )
})


test_that("filtering thresholds are validated with ibdmdb", {
  expect_error(
    filterRPCAInput(
      filter_mae[[1L]],
      assay.type = raw_assay_types[[1L]],
      min.sample.count = -1
    ),
    "non-negative"
  )
  
  expect_error(
    filterRPCAInput(
      filter_mae[[1L]],
      assay.type = raw_assay_types[[1L]],
      min.feature.count = c(1, 2)
    ),
    "single"
  )
  
  expect_error(
    filterRPCAInput(
      filter_mae[[1L]],
      assay.type = raw_assay_types[[1L]],
      min.sample.count = Inf
    ),
    "finite"
  )
  
  expect_error(
    filterRPCAInput(
      filter_mae[[1L]],
      assay.type = raw_assay_types[[1L]],
      min.feature.frequency = 101
    ),
    "between 0 and 100"
  )
})


test_that("MultiAssayExperiment input requires shared samples", {
  shared_samples <- intersect(
    colnames(filter_mae[[1L]]),
    colnames(filter_mae[[2L]])
  )
  
  expect_true(length(shared_samples) >= 4L)
  
  first_se <- filter_mae[[1L]][
    ,
    shared_samples[1:2],
    drop = FALSE
  ]
  
  second_se <- filter_mae[[2L]][
    ,
    shared_samples[3:4],
    drop = FALSE
  ]
  
  no_overlap_mae <- MultiAssayExperiment::MultiAssayExperiment(
    experiments = list(
      first = first_se,
      second = second_se
    )
  )
  
  expect_error(
    filterRPCAInput(
      no_overlap_mae,
      experiments = c("first", "second"),
      assay.types = raw_assay_types,
      min.sample.count = 0,
      min.feature.count = 0,
      min.feature.frequency = 0
    ),
    "No samples overlap"
  )
})
