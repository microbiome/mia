context("subsetBy")

data(GlobalPatterns, package = "mia")
# Plain SE keeps the tests fast.
tse <- SummarizedExperiment(
    assays = assays(GlobalPatterns), rowData = rowData(GlobalPatterns),
    colData = colData(GlobalPatterns))

test_that("subsetByAbundance", {
    mat <- assay(tse, "counts")
    # Features absent from all samples are removed by default
    res <- subsetByAbundance(tse)
    expect_identical(rownames(res), rownames(tse)[rowSums(mat) > 0])
    expect_identical(colnames(res), colnames(tse))
    # Threshold is strict
    res <- subsetByAbundance(tse, threshold = 100)
    expect_identical(rownames(res), rownames(tse)[rowSums(mat) > 100])
    res <- subsetByAbundance(tse, threshold = min(rowSums(mat)))
    expect_false(any(rowSums(assay(res)) == min(rowSums(mat))))
    expect_equal(nrow(subsetByAbundance(tse, threshold = max(rowSums(mat)))), 0)
    # Samples
    res <- subsetByAbundance(tse, threshold = 1e5, MARGIN = "samples")
    expect_identical(colnames(res), colnames(tse)[colSums(mat) > 1e5])
    expect_identical(rownames(res), rownames(tse))
    expect_identical(
        dimnames(subsetByAbundance(tse, threshold = 1e5, MARGIN = 2L)),
        dimnames(res))
    # Non-integer threshold and other assay
    tse2 <- transformAssay(tse, method = "relabundance")
    res <- subsetByAbundance(
        tse2, assay.type = "relabundance", threshold = 0.01)
    expect_identical(
        rownames(res),
        rownames(tse)[rowSums(assay(tse2, "relabundance")) > 0.01])
    # Other statistics
    rel <- assay(tse2, "relabundance")
    res <- subsetByAbundance(
        tse2, assay.type = "relabundance", method = "mean", threshold = 0.001)
    expect_identical(rownames(res), rownames(tse)[rowMeans(rel) > 0.001])
    res <- subsetByAbundance(
        tse2, assay.type = "relabundance", method = "max", threshold = 0.01)
    expect_identical(
        rownames(res), rownames(tse)[apply(rel, 1, max) > 0.01])
    res <- subsetByAbundance(tse, method = "median", threshold = 1)
    expect_identical(
        rownames(res), rownames(tse)[apply(mat, 1, median) > 1])
    res <- subsetByAbundance(
        tse, method = "mean", threshold = 1e5, MARGIN = "samples")
    expect_identical(colnames(res), colnames(tse)[colMeans(mat) > 1e5])
    # Missing values
    tse_na <- tse
    assay(tse_na, "counts")[1, 1] <- NA
    res <- subsetByAbundance(tse_na)
    expect_false(rownames(tse)[1] %in% rownames(res))
    res <- subsetByAbundance(tse_na, MARGIN = "samples")
    expect_false(colnames(tse)[1] %in% colnames(res))
    res <- subsetByAbundance(tse_na, na.rm = TRUE)
    expect_true(rownames(tse)[1] %in% rownames(res))
    # Works with DelayedArray and TreeSE
    se <- tse[1:50, ]
    assay(se, "counts") <- DelayedArray(assay(se, "counts"))
    expect_identical(
        dimnames(subsetByAbundance(se, threshold = 10)),
        dimnames(subsetByAbundance(tse[1:50, ], threshold = 10)))
    res <- subsetByAbundance(GlobalPatterns[1:50, ], threshold = 10)
    expect_s4_class(res, "TreeSummarizedExperiment")
    # Errors
    expect_error(subsetByAbundance(tse, assay.type = "foo"))
    expect_error(subsetByAbundance(tse, threshold = -1))
    expect_error(subsetByAbundance(tse, threshold = "1"))
    expect_error(subsetByAbundance(tse, threshold = TRUE))
    expect_error(subsetByAbundance(tse, threshold = c(1, 2)))
    expect_error(subsetByAbundance(tse, threshold = NA_real_))
    expect_error(subsetByAbundance(tse, MARGIN = 3))
    expect_error(subsetByAbundance(tse, na.rm = "TRUE"))
    expect_error(subsetByAbundance(tse, na.rm = NA))
    expect_error(subsetByAbundance(tse, method = "foo"))
    expect_error(subsetByAbundance(tse, method = c("sum", "mean")))
})

test_that("subsetByVariance", {
    mat <- matrix(
        c(1, 1, 1, 1,
          1, 2, 3, 4,
          0, 0, 0, 0,
          1, 2, 1, 2),
        nrow = 4, byrow = TRUE,
        dimnames = list(paste0("f", 1:4), paste0("s", 1:4)))
    se <- SummarizedExperiment(assays = list(counts = mat))
    vars <- apply(mat, 1, var)
    # Zero-variance features are removed by default
    res <- subsetByVariance(se)
    expect_identical(rownames(res), rownames(mat)[vars > 0])
    # Threshold
    res <- subsetByVariance(se, threshold = 0.5)
    expect_identical(rownames(res), rownames(mat)[vars > 0.5])
    # Samples
    res <- subsetByVariance(se, MARGIN = "samples")
    expect_identical(colnames(res), colnames(mat)[apply(mat, 2, var) > 0])
    expect_identical(
        dimnames(subsetByVariance(se, MARGIN = 2L)), dimnames(res))
    # Other assay
    assay(se, "test") <- mat * 0
    expect_equal(nrow(subsetByVariance(se, assay.type = "test")), 0L)
    # Missing values
    assay(se, "na") <- mat
    assay(se, "na")[2, 1] <- NA
    res <- subsetByVariance(se, assay.type = "na")
    expect_false("f2" %in% rownames(res))
    res <- subsetByVariance(se, assay.type = "na", na.rm = TRUE)
    expect_true("f2" %in% rownames(res))
    # Works with DelayedArray
    se2 <- se
    assay(se2, "counts") <- DelayedArray(mat)
    expect_identical(
        dimnames(subsetByVariance(se2)), dimnames(subsetByVariance(se)))
    # Errors
    expect_error(subsetByVariance(se, assay.type = "foo"))
    expect_error(subsetByVariance(se, threshold = -1))
    expect_error(subsetByVariance(se, threshold = "a"))
    expect_error(subsetByVariance(se, threshold = c(1, 2)))
    expect_error(subsetByVariance(se, MARGIN = 3))
    expect_error(subsetByVariance(se, na.rm = "TRUE"))
    # Works with TreeSE
    tse <- GlobalPatterns[1:50, 1:5]
    res <- subsetByVariance(tse)
    expect_s4_class(res, "TreeSummarizedExperiment")
    expect_true(nrow(res) < nrow(tse))
})

test_that("subsetBySd", {
    mat <- matrix(
        c(1, 1, 1, 1,
          1, 2, 3, 4,
          0, 0, 0, 0,
          1, 2, 1, 2),
        nrow = 4, byrow = TRUE,
        dimnames = list(paste0("f", 1:4), paste0("s", 1:4)))
    se <- SummarizedExperiment(assays = list(counts = mat))
    sds <- apply(mat, 1, sd)
    # Zero-sd features are removed by default, same as with variance
    res <- subsetBySd(se)
    expect_identical(rownames(res), rownames(mat)[sds > 0])
    expect_identical(dimnames(res), dimnames(subsetByVariance(se)))
    # Threshold is on the sd scale: sd(f4) ~ 0.58 > 0.5 but var(f4) = 1/3 < 0.5
    res <- subsetBySd(se, threshold = 0.5)
    expect_identical(rownames(res), rownames(mat)[sds > 0.5])
    expect_true("f4" %in% rownames(res))
    # Samples
    res <- subsetBySd(se, MARGIN = "samples")
    expect_identical(colnames(res), colnames(mat)[apply(mat, 2, sd) > 0])
    # Missing values
    assay(se, "na") <- mat
    assay(se, "na")[2, 1] <- NA
    expect_false("f2" %in% rownames(subsetBySd(se, assay.type = "na")))
    expect_true(
        "f2" %in% rownames(subsetBySd(se, assay.type = "na", na.rm = TRUE)))
    # Works with DelayedArray and TreeSE
    se2 <- se
    assay(se2, "counts") <- DelayedArray(mat)
    expect_identical(dimnames(subsetBySd(se2)), dimnames(subsetBySd(se)))
    expect_s4_class(
        subsetBySd(GlobalPatterns[1:50, 1:5]), "TreeSummarizedExperiment")
    # Errors
    expect_error(subsetBySd(se, assay.type = "foo"))
    expect_error(subsetBySd(se, threshold = -1))
    expect_error(subsetBySd(se, MARGIN = 3))
    expect_error(subsetBySd(se, na.rm = "TRUE"))
})

test_that("subsetByMissing", {
    mat <- matrix(
        c(1,  2,  3,  4,  5,
          NA, 2,  3,  4,  5,
          NA, NA, 3,  4,  5,
          NA, NA, NA, NA, 5,
          1,  2,  NA, NA, NA),
        nrow = 5, byrow = TRUE,
        dimnames = list(paste0("f", 1:5), paste0("s", 1:5)))
    se <- SummarizedExperiment(
        assays = list(counts = mat),
        colData = DataFrame(group = c("a", "a", "b", "b", "b")))
    prop <- rowMeans(is.na(mat))
    # Default: at most 20% missing
    res <- subsetByMissing(se)
    expect_identical(rownames(res), rownames(mat)[prop <= 0.2])
    expect_identical(rownames(res), c("f1", "f2"))
    # Threshold
    expect_identical(rownames(subsetByMissing(se, threshold = 0)), "f1")
    expect_identical(
        rownames(subsetByMissing(se, threshold = 0.6)),
        rownames(mat)[prop <= 0.6])
    expect_identical(dim(subsetByMissing(se, threshold = 1)), dim(se))
    # Samples
    res <- subsetByMissing(se, threshold = 0.4, MARGIN = "samples")
    expect_identical(
        colnames(res), colnames(mat)[colMeans(is.na(mat)) <= 0.4])
    # Groups: f3 is complete in group b, f5 is complete in group a, f4 is
    # mostly missing in both groups
    res <- subsetByMissing(se, col.var = "group")
    expect_identical(rownames(res), c("f1", "f2", "f3", "f5"))
    expect_identical(colnames(res), colnames(se))
    # Works with DelayedArray and TreeSE
    se2 <- se
    assay(se2, "counts") <- DelayedArray(mat)
    expect_identical(
        dimnames(subsetByMissing(se2, col.var = "group")), dimnames(res))
    expect_s4_class(
        subsetByMissing(GlobalPatterns[1:20, ]), "TreeSummarizedExperiment")
    # Errors
    expect_error(subsetByMissing(se, assay.type = "foo"))
    expect_error(subsetByMissing(se, threshold = 1.5))
    expect_error(subsetByMissing(se, threshold = -1))
    expect_error(subsetByMissing(se, threshold = "0.2"))
    expect_error(subsetByMissing(se, MARGIN = 3))
    expect_error(subsetByMissing(se, col.var = "foo"))
    expect_error(subsetByMissing(se, col.var = "group", MARGIN = "samples"))
    se$group[1] <- NA
    expect_error(subsetByMissing(se, col.var = "group"))
})

test_that("subsetByTop", {
    tse <- GlobalPatterns[1:50, ]
    for( method in c("mean", "sum", "median", "sd", "cv", "prevalence") ){
        ref <- getTop(tse, top = 5, method = method)
        res <- subsetByTop(tse, top = 5, method = method)
        expect_setequal(rownames(res), ref)
        # Original order is preserved
        expect_identical(rownames(res), rownames(tse)[rownames(tse) %in% ref])
    }
    # Arguments are passed to getTop
    ref <- getTop(tse, top = 3, method = "prevalence", detection = 10)
    res <- subsetByTop(tse, top = 3, method = "prevalence", detection = 10)
    expect_setequal(rownames(res), ref)
    # Agglomeration
    ref <- getTop(tse, top = 3, rank = "Family")
    res <- subsetByTop(tse, top = 3, rank = "Family")
    expect_setequal(rownames(res), ref)
    # Duplicated rownames
    se <- SummarizedExperiment(assays = list(
        counts = matrix(c(1, 5, 3, 2, 6, 4), nrow = 3,
            dimnames = list(c("a", "a", "b"), NULL))))
    res <- subsetByTop(se, top = 2)
    expect_equal(assay(res), assay(se)[2:3, ], check.attributes = FALSE)
    # Errors
    expect_error(subsetByTop(tse, top = 1000))
    expect_error(subsetByTop(tse, method = "foo"))
})
