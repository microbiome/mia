context("Clustering rowData")
test_that("Clustering rowData", {
    data("GlobalPatterns")
    tse <- GlobalPatterns
    tse <- agglomerateByRank(tse, rank="Class")
    # Checking output
    tse <- rowClust(tse, k=6)
    expect_true("clust_id" %in% colnames(rowData(tse)))
    expect_true(is.factor(rowData(tse)$clust_id))
    expect_equal(length(unique(rowData(tse)$clust_id)),6)
    # Checking args
    expect_error(rowClust(tse, abund_values = "assay_not_present"))
    expect_error(rowClust(tse, dissimilarity_FUN = "FUN_not_char"))
    expect_error(rowClust(tse, clustering_method = "non-existing-method"))
    expect_error(rowClust(tse, param = "arg_not_in_cutree"))
})