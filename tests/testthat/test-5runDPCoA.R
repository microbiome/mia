context("runDPCoA")
test_that("runDPCoA", {
    data(esophagus)
    #
    esophagus <- runDPCoA(esophagus)
    expect_named(reducedDims(esophagus),"DPCoA")
    expect_true(is.matrix(reducedDim(esophagus,"DPCoA")))
    expect_equal(dim(reducedDim(esophagus,"DPCoA")),c(3,2))
    red <- reducedDim(esophagus,"DPCoA")
    expect_equal(names(attributes(red)),
                 c("dim","dimnames","eig","sample_red","feature_weights",
                   "sample_weights" ))
    expect_equal(dim(attr(red,"sample_red")),c(58,2))
    
    # ERRORs
    expect_error(
        runDPCoA(esophagus, assay.type = "test", tree_name = "phylo", ncomponents = 2, ntop = NULL,
                 subset_row = NULL, scale = FALSE, transposed = FALSE)
    )
    expect_error(
        runDPCoA(esophagus, assay.type = 1, tree_name = "phylo", ncomponents = 2, ntop = NULL,
                 subset_row = NULL, scale = FALSE, transposed = FALSE)
    )
    expect_error(
        runDPCoA(esophagus, assay.type = TRUE, tree_name = "phylo", ncomponents = 2, ntop = NULL,
                 subset_row = NULL, scale = FALSE, transposed = FALSE)
    )
    expect_error(
        runDPCoA(esophagus, assay.type = "counts", tree_name = "test", ncomponents = 2, ntop = NULL,
                 subset_row = NULL, scale = FALSE, transposed = FALSE)
    )
    expect_error(
        runDPCoA(esophagus, assay.type = "counts", tree_name = 1, ncomponents = 2, ntop = NULL,
                 subset_row = NULL, scale = FALSE, transposed = FALSE)
    )
    expect_error(
        runDPCoA(esophagus, assay.type = "counts", tree_name = "phylo", ncomponents = TRUE, ntop = NULL,
                 subset_row = NULL, scale = FALSE, transposed = FALSE)
    )
    expect_error(
        runDPCoA(esophagus, assay.type = "counts", tree_name = "phylo", ncomponents = "test", ntop = NULL,
                 subset_row = NULL, scale = FALSE, transposed = FALSE)
    )
    expect_error(
        runDPCoA(esophagus, assay.type = "counts", tree_name = "phylo", ncomponents = "test", ntop = "test",
                 subset_row = NULL, scale = FALSE, transposed = FALSE)
    )
    expect_error(
        runDPCoA(esophagus, assay.type = "counts", tree_name = "phylo", ncomponents = 1.3, ntop = "test",
                 subset_row = NULL, scale = FALSE, transposed = FALSE)
    )
    expect_error(
        runDPCoA(esophagus, name = c("test", "test2"), assay.type = "counts",
                 tree_name = "phylo", ncomponents = 1.3, ntop = "test",
                 subset_row = NULL, scale = FALSE, transposed = FALSE)
    )
    expect_error(
        runDPCoA(esophagus, name = 1, assay.type = "counts", tree_name = "phylo", ncomponents = 1.3, ntop = "test",
                 subset_row = NULL, scale = FALSE, transposed = FALSE)
    )
    
    data("GlobalPatterns")
    tse <- mergeSEs(esophagus, GlobalPatterns)
    # expect_warning(runDPCoA(tse))
    # expect_warning(runDPCoA(tse, tree_name = "phylo.1"))
})
