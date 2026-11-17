context("subsetBy")

test_that("subsetBy", {
    data(GlobalPatterns, package = "mia")
    tse <- GlobalPatterns
    mat <- assay(tse, "counts")
    # Rows
    res <- subsetBy(tse, rows = rowSums(counts) > 10)
    expect_s4_class(res, "TreeSummarizedExperiment")
    expect_identical(rownames(res), rownames(tse)[rowSums(mat) > 10])
    expect_identical(ncol(res), ncol(tse))
    # Columns
    res <- subsetBy(tse, cols = colSums(counts) > 1e6)
    expect_identical(colnames(res), colnames(tse)[colSums(mat) > 1e6])
    expect_identical(nrow(res), nrow(tse))
    # rowData and colData variables, and variables from the caller
    th <- "Feces"
    res <- subsetBy(
        tse, rows = Kingdom == "Bacteria", cols = SampleType == th)
    expect_true(all(rowData(res)$Kingdom == "Bacteria"))
    expect_true(all(res$SampleType == "Feces"))
    # Stepwise filtering computes each condition from the subsetted data
    res <- tse |>
        subsetBy(rows = rowSums(counts) > 10) |>
        subsetBy(rows = rowMeans(counts > 0) > 0.25) |>
        subsetBy(cols = colSums(counts) > 100000)
    ref <- tse
    ref <- ref[rowSums(assay(ref, "counts")) > 10, ]
    ref <- ref[rowMeans(assay(ref, "counts") > 0) > 0.25, ]
    ref <- ref[, colSums(assay(ref, "counts")) > 100000]
    expect_identical(dimnames(res), dimnames(ref))
    # cols is evaluated after rows
    res <- subsetBy(
        tse, rows = rowSums(counts) > 10, cols = colSums(counts) > 100000)
    sub <- tse[rowSums(mat) > 10, ]
    expect_identical(
        colnames(res),
        colnames(sub)[colSums(assay(sub, "counts")) > 100000])
    # NULL keeps everything
    expect_identical(dim(subsetBy(tse)), dim(tse))
    # NA is treated as FALSE
    res <- subsetBy(tse, rows = ifelse(rowSums(counts) > 10, TRUE, NA))
    expect_identical(nrow(res), sum(rowSums(mat) > 10))
    # Errors
    expect_error(subsetBy(tse, rows = 1))
    expect_error(subsetBy(tse, rows = c(TRUE, FALSE)))
    expect_error(subsetBy(tse, cols = rowSums(counts) > 0))
})
