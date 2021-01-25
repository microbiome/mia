context("relabundance")
test_that("relabundance", {
    mat <- matrix(1:60, nrow = 6)
    df <- DataFrame(n = c(1:6))
    expect_error(relAbundanceCounts(SummarizedExperiment(assays = list(mat = mat),
                                                         rowData = df)),
                 "'abund_values' must be a valid name of assays")
    se <- SummarizedExperiment(assays = list(counts = mat),
                               rowData = df)
    expect_error(relAbundanceCounts(se, name = FALSE),
                 "'name' must be a non-empty single character value")
    actual <- relAbundanceCounts(se)
    expect_named(assays(actual), c("counts", "relabundance"))
    expect_equal(assay(actual,"relabundance")[,1],
                 seq.int(1,6)/21)
    expect_equal(assay(actual,"relabundance"),
                 relabundance(actual))
    rel_mat <- relabundance(actual)
    f <- rev(seq_len(ncol(relabundance(actual))))
    relabundance(actual) <- relabundance(actual)[,f]
    expect_equal(rel_mat,
                 relabundance(actual)[,f])
})
