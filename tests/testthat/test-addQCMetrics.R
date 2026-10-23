test_that("getQCMetrics and addQCMetrics", {
    data(GlobalPatterns, package = "mia")
    tse <- GlobalPatterns
    mat <- assay(tse, "counts")
    # Samples
    res <- getQCMetrics(tse)
    expect_s4_class(res, "DataFrame")
    expect_identical(rownames(res), colnames(tse))
    expect_equal(res[["total"]], unname(colSums(mat)))
    expect_equal(res[["detected"]], unname(colSums(mat > 0)))
    expect_equal(res[["mean"]], unname(colMeans(mat)))
    expect_equal(res[["median"]], unname(apply(mat, 2, median)))
    expect_equal(res[["sd"]], unname(apply(mat, 2, sd)))
    # Features
    res <- getQCMetrics(tse, by = "rows")
    expect_identical(rownames(res), rownames(tse))
    expect_equal(res[["total"]], unname(rowSums(mat)))
    # Results are added to colData or rowData
    tse <- addQCMetrics(tse)
    expect_equal(tse[["total"]], unname(colSums(mat)))
    tse <- addQCMetrics(tse, by = "rows")
    expect_equal(rowData(tse)[["detected"]], unname(rowSums(mat > 0)))
    # Names of the metrics can be changed
    tse <- addQCMetrics(tse, name = paste0("qc_", 1:5))
    expect_equal(tse[["qc_1"]], unname(colSums(mat)))
    expect_error(getQCMetrics(tse, name = "total"))
    # Metrics are calculated and added for altExp
    altExp(tse, "genus") <- agglomerateByRank(GlobalPatterns, rank = "Genus")
    tse <- addQCMetrics(tse, altexp = "genus")
    expect_equal(
        altExp(tse, "genus")[["total"]],
        unname(colSums(assay(altExp(tse, "genus"), "counts"))))
    # Input check
    expect_error(getQCMetrics(tse, assay.type = "test"))
    expect_error(getQCMetrics(tse, by = 3))
})
