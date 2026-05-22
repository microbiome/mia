# Check that input checks are working
test_that("addMDS errors", {
    tse <- makeTSE()
    assayNames(tse) <- "counts"
    #
    getMDS(tse, assay.type = "test") |> expect_error()
    addMDS(tse, name = 1) |> expect_error()
    addMDS(tse, name = c("test", "test2")) |> expect_error()
    addMDS(tse, name = NULL) |> expect_error()
})

test_that("Compare addMDS with runMDS", {
    data("GlobalPatterns")
    tse <- GlobalPatterns
    #
    res <- getMDS(tse, assay.type = "counts", method = "bray")
    ref <- calculateMDS(tse, assay.type = "counts", method = "bray", FUN = getDissimilarity)
    expect_equal(res, ref)
    #
    res <- addMDS(tse, assay.type = "counts", method = "bray", name = "test")
    ref <- runMDS(tse, assay.type = "counts", method = "bray", FUN = getDissimilarity, name = "test")
    expect_equal(reducedDim(res, "test"), reducedDim(ref, "test"))
    #
    res <- getMDS(tse, assay.type = "counts", method = "unifrac", name = "test") |> expect_warning()
    ref <- calculateMDS(tse, assay.type = "counts", method = "unifrac", tree = rowTree(tse), FUN = getDissimilarity, name = "test") |> expect_warning()
    expect_equal(res, ref)
})

test_that("Check rarefaction", {
    data("GlobalPatterns")
    tse <- GlobalPatterns
    #
    sample <- assay(tse) |> colSums() |> min()
    sample <- sample+1
    res <- addMDS(tse, assay.type = "counts", method = "bray", niter = 2L, sample = sample) |> expect_warning()
    expect_true( ncol(res) == sum(colSums(assay(tse)) >= sample) )
    #
    res <- addMDS(tse, assay.type = "counts", method = "bray", niter = 2L, sample = sample, subset.result = FALSE) |> expect_warning()
    expect_true( ncol(res) == ncol(tse) )
})

test_that("Test named vector as node labels", {
    data("GlobalPatterns")
    tse <- GlobalPatterns
    #
    named_vector <- rowLinks(tse)[["nodeLab"]]
    names(named_vector) <- rownames(tse)
    #
    tse <- tse[1:100, ]
    vector <- rowLinks(tse)[["nodeLab"]]
    #
    res1 <- getMDS(
        tse,
        assay.type = "counts",
        method = "unifrac",
        tree = rowTree(tse),
        node.label = named_vector
    ) |> expect_warning()
    res2 <- getMDS(
        tse,
        assay.type = "counts",
        method = "unifrac",
        tree = rowTree(tse),
        node.label = named_vector
    ) |> expect_warning()
    expect_equal(res1, res2)
})
