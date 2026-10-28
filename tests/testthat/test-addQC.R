test_that("getQC and addQC", {
    data(GlobalPatterns, package = "mia")
    tse <- GlobalPatterns
    mat <- assay(tse, "counts")
    # Samples
    res <- getQC(tse)
    expect_s4_class(res, "DataFrame")
    expect_identical(rownames(res), colnames(tse))
    expect_equal(res[["total"]], unname(colSums(mat)))
    expect_equal(res[["detected"]], unname(colSums(mat > 0)))
    expect_equal(res[["mean"]], unname(colMeans(mat)))
    expect_equal(res[["sd"]], unname(apply(mat, 2, sd)))
    expect_equal(res[["q0"]], unname(apply(mat, 2, min)))
    expect_equal(res[["q50"]], unname(apply(mat, 2, median)))
    expect_equal(res[["q100"]], unname(apply(mat, 2, max)))
    # Features
    res <- getQC(tse, by = "rows")
    expect_identical(rownames(res), rownames(tse))
    expect_equal(res[["total"]], unname(rowSums(mat)))
    # Results are added to colData or rowData
    tse <- addQC(tse)
    expect_equal(tse[["total"]], unname(colSums(mat)))
    tse <- addQC(tse, by = "rows")
    expect_equal(rowData(tse)[["detected"]], unname(rowSums(mat > 0)))
    # Names of the metrics can be changed
    tse <- addQC(tse, name = paste0("qc_", 1:9))
    expect_equal(tse[["qc_1"]], unname(colSums(mat)))
    expect_error(getQC(tse, name = "total"))
    # Metrics are calculated and added for altExp
    altExp(tse, "genus") <- agglomerateByRank(GlobalPatterns, rank = "Genus")
    tse <- addQC(tse, altexp = "genus")
    expect_equal(
        altExp(tse, "genus")[["total"]],
        unname(colSums(assay(altExp(tse, "genus"), "counts"))))
    # Input check
    expect_error(getQC(tse, assay.type = "test"))
    expect_error(getQC(tse, by = 3))
})
