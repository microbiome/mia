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

test_that("subsetBy with collapse", {
    data(GlobalPatterns, package = "mia")
    tse <- GlobalPatterns
    mat <- assay(tse, "counts")
    keep <- rowSums(mat) > 1e5
    res <- subsetBy(tse, rows = rowSums(counts) > 1e5, collapse = TRUE)
    # Kept features in the original order, then the merged feature
    expect_identical(rownames(res), c(rownames(tse)[keep], "Other"))
    expect_equal(
        unname(assay(res, "counts")["Other", ]),
        unname(colSums(mat[!keep, ])))
    # Sample totals do not change
    expect_equal(colSums(assay(res, "counts")), colSums(mat))
    # Kept features are unchanged
    expect_equal(
        assay(res, "counts")[rownames(tse)[keep], ], mat[keep, ])
    # Custom name
    res <- subsetBy(
        tse, rows = rowSums(counts) > 1e5, collapse = TRUE,
        collapse.name = "Rest")
    expect_identical(tail(rownames(res), 1), "Rest")
    # Nothing to merge: the object is returned as subsetted
    res <- subsetBy(tse, rows = rowSums(counts) >= 0, collapse = TRUE)
    expect_identical(dim(res), dim(tse))
    # cols is evaluated after collapsing, so totals include the merged feature
    res <- subsetBy(
        tse, rows = rowSums(counts) > 1e5, cols = colSums(counts) > 1e6,
        collapse = TRUE)
    expect_identical(colnames(res), colnames(tse)[colSums(mat) > 1e6])
    # Errors
    expect_error(subsetBy(tse, rows = rowSums(counts) > 0, collapse = 1))
    expect_error(
        subsetBy(tse, rows = rowSums(counts) > 0, collapse = TRUE,
            collapse.name = c("a", "b")))
    tse2 <- tse
    rownames(tse2)[1] <- "Other"
    expect_error(
        subsetBy(tse2, rows = rowSums(counts) > 1e5 | rownames(tse2) ==
            "Other", collapse = TRUE))
})

test_that("subsetBy with stored quality control metrics", {
    data(GlobalPatterns, package = "mia")
    tse <- GlobalPatterns
    tse <- addQC(tse, by = "rows")
    tse <- addQC(tse, by = "cols")
    # rowData and colData metrics
    res <- subsetBy(tse, rows = sd > 0 & detected >= 3, cols = total > 1e6)
    ref <- tse[
        rowData(tse)$sd > 0 & rowData(tse)$detected >= 3, tse$total > 1e6]
    expect_identical(dimnames(res), dimnames(ref))
    # Stored metrics combined with assays and collapse
    keep <- rowData(tse)$sd > 0 & rowSums(assay(tse, "counts") > 0) >= 3
    res <- subsetBy(
        tse, rows = sd > 0 & rowSums(counts > 0) >= 3, collapse = TRUE)
    expect_identical(rownames(res), c(rownames(tse)[keep], "Other"))
    # Stored metrics are not updated after subsetting, unlike metrics computed
    # in the condition
    sub <- subsetBy(tse, rows = rowSums(counts) > 1e5)
    expect_identical(sub$total, tse$total)
    expect_false(
        identical(
            colnames(subsetBy(sub, cols = total > 1e5)),
            colnames(subsetBy(sub, cols = colSums(counts) > 1e5))))
    # Assays take precedence over rowData/colData fields with the same name
    tse2 <- tse
    rowData(tse2)$counts <- rep(-1, nrow(tse2))
    res <- subsetBy(tse2, rows = rowSums(counts) > 1e5)
    expect_identical(
        rownames(res), rownames(tse)[rowSums(assay(tse, "counts")) > 1e5])
})
