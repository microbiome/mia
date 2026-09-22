context("filterRPCAInput")

data(GlobalPatterns, package = "mia")
tse <- GlobalPatterns

test_that("filterRPCAInput filters samples and features correctly", {
    count_matrix <- assay(tse, "counts")

    sample_keep <- colSums(count_matrix) > 10000
    feature_keep <- rowSums(count_matrix) > 100
    feature_keep <- feature_keep &
        (rowSums(count_matrix > 0) / ncol(count_matrix) > 0.10)

    filtered_tse <- filterRPCAInput(
        tse,
        min.sample.count = 10000,
        min.feature.count = 100,
        min.feature.frequency = 0.10
    )

    expect_identical(
        rownames(filtered_tse),
        rownames(tse)[feature_keep]
    )

    expect_identical(
        colnames(filtered_tse),
        colnames(tse)[sample_keep]
    )
})

test_that("Each filtering criterion can be disabled with NULL", {
    filtered_tse <- filterRPCAInput(
        tse,
        min.sample.count = NULL,
        min.feature.count = NULL,
        min.feature.frequency = NULL
    )

    expect_identical(dimnames(filtered_tse), dimnames(tse))

    filtered_tse <- filterRPCAInput(
        tse,
        min.sample.count = NULL,
        min.feature.count = 100,
        min.feature.frequency = NULL
    )

    expect_identical(
        rownames(filtered_tse),
        rownames(tse)[rowSums(assay(tse, "counts")) > 100]
    )

    filtered_tse <- filterRPCAInput(
        tse,
        min.sample.count = NULL,
        min.feature.count = NULL,
        min.feature.frequency = 0.10
    )

    expect_identical(
        rownames(filtered_tse),
        rownames(tse)[
            rowSums(assay(tse, "counts") > 0) / ncol(tse) > 0.10
        ]
    )

    filtered_tse <- filterRPCAInput(
        tse,
        min.sample.count = 10000,
        min.feature.count = NULL,
        min.feature.frequency = NULL
    )

    expect_identical(
        colnames(filtered_tse),
        colnames(tse)[colSums(assay(tse, "counts")) > 10000]
    )
})

test_that("Strict greater-than comparison is used", {
    count_matrix <- assay(tse, "counts")

    sample_threshold <- min(colSums(count_matrix))
    feature_threshold <- min(rowSums(count_matrix))
    frequency_threshold <- min(rowSums(count_matrix > 0) / ncol(count_matrix))

    filtered_tse <- filterRPCAInput(
        tse,
        min.sample.count = NULL,
        min.feature.count = feature_threshold,
        min.feature.frequency = frequency_threshold
    )
    expect_false(any(rowSums(assay(filtered_tse, "counts")) == feature_threshold))

    filtered_tse <- filterRPCAInput(
        tse,
        min.sample.count = sample_threshold,
        min.feature.count = NULL,
        min.feature.frequency = NULL
    )
    expect_false(any(colSums(assay(filtered_tse, "counts")) == sample_threshold))
})

test_that("SummarizedExperiment compatibility is preserved", {
    se <- as(tse, "SummarizedExperiment")

    filtered_se <- filterRPCAInput(
        se,
        min.sample.count = 10000,
        min.feature.count = 100,
        min.feature.frequency = 0.10
    )

    expect_s4_class(filtered_se, "SummarizedExperiment")
})

test_that("Boundary values are handled correctly", {
    filtered_tse <- filterRPCAInput(
        tse,
        min.sample.count = NULL,
        min.feature.count = NULL,
        min.feature.frequency = 1
    )

    expect_equal(nrow(filtered_tse), 0)
})

test_that("Invalid assay name is rejected", {
    expect_error(
        filterRPCAInput(
            tse,
            assay.type = "foo"
        )
    )
})

test_that("min.sample.count is validated", {
    expect_error(filterRPCAInput(tse, min.sample.count = -1))
    expect_error(filterRPCAInput(tse, min.sample.count = 1.5))
    expect_error(filterRPCAInput(tse, min.sample.count = TRUE))
    expect_error(filterRPCAInput(tse, min.sample.count = "1"))
    expect_error(filterRPCAInput(tse, min.sample.count = c(1L, 2L)))
})

test_that("min.feature.count is validated", {
    expect_error(filterRPCAInput(tse, min.feature.count = -1))
    expect_error(filterRPCAInput(tse, min.feature.count = 1.5))
    expect_error(filterRPCAInput(tse, min.feature.count = TRUE))
    expect_error(filterRPCAInput(tse, min.feature.count = "1"))
    expect_error(filterRPCAInput(tse, min.feature.count = c(1L, 2L)))
})

test_that("min.feature.frequency is validated", {
    expect_error(filterRPCAInput(tse, min.feature.frequency = -0.1))
    expect_error(filterRPCAInput(tse, min.feature.frequency = 1.1))
    expect_error(filterRPCAInput(tse, min.feature.frequency = TRUE))
    expect_error(filterRPCAInput(tse, min.feature.frequency = "0.1"))
    expect_error(filterRPCAInput(tse, min.feature.frequency = c(0.1, 0.2)))
})

test_that("na.rm is validated", {
    expect_error(filterRPCAInput(tse, na.rm = 1))
    expect_error(filterRPCAInput(tse, na.rm = "TRUE"))
    expect_error(filterRPCAInput(tse, na.rm = NA))
})

test_that("na.rm argument is accepted", {
    filtered_tse <- filterRPCAInput(tse, na.rm = FALSE)

    expect_s4_class(filtered_tse, class(tse))
})

test_that("NA values are handled correctly", {

  tse_na <- tse
  assay(tse_na, "counts")[1, 1] <- NA

  ## With na.rm = FALSE, the affected sample and feature are removed.
  filtered_tse <- filterRPCAInput(
      tse_na,
      na.rm = FALSE
  )

  expect_false(rownames(tse_na)[1] %in% rownames(filtered_tse))
  expect_false(colnames(tse_na)[1] %in% colnames(filtered_tse))

  ## With na.rm = TRUE, they are retained because the NA is ignored.
  filtered_tse <- filterRPCAInput(
      tse_na,
      na.rm = TRUE
  )

  expect_true(rownames(tse_na)[1] %in% rownames(filtered_tse))
  expect_true(colnames(tse_na)[1] %in% colnames(filtered_tse))
})
