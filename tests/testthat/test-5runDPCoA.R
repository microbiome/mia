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
        runDPCoA(esophagus, assay_name = "test", tree_name = "phylo", dimred = NULL, 
                 n_dimred = NULL, ncomponents = 2, ntop = NULL,
                 subset_row = NULL, scale = FALSE, transposed = FALSE)
    )
    expect_error(
        runDPCoA(esophagus, assay_name = 1, tree_name = "phylo", dimred = NULL, 
                 n_dimred = NULL, ncomponents = 2, ntop = NULL,
                 subset_row = NULL, scale = FALSE, transposed = FALSE)
    )
    expect_error(
        runDPCoA(esophagus, assay_name = TRUE, tree_name = "phylo", dimred = NULL, 
                 n_dimred = NULL, ncomponents = 2, ntop = NULL,
                 subset_row = NULL, scale = FALSE, transposed = FALSE)
    )
    expect_error(
        runDPCoA(esophagus, assay_name = "counts", tree_name = "test", dimred = NULL, 
                 n_dimred = NULL, ncomponents = 2, ntop = NULL,
                 subset_row = NULL, scale = FALSE, transposed = FALSE)
    )
    expect_error(
        runDPCoA(esophagus, assay_name = "counts", tree_name = 1, dimred = NULL, 
                 n_dimred = NULL, ncomponents = 2, ntop = NULL,
                 subset_row = NULL, scale = FALSE, transposed = FALSE)
    )
    expect_error(
        runDPCoA(esophagus, assay_name = "counts", tree_name = "phylo", dimred = NULL, 
                 n_dimred = NULL, ncomponents = TRUE, ntop = NULL,
                 subset_row = NULL, scale = FALSE, transposed = FALSE)
    )
    expect_error(
        runDPCoA(esophagus, assay_name = "counts", tree_name = "phylo", dimred = NULL, 
                 n_dimred = NULL, ncomponents = "test", ntop = NULL,
                 subset_row = NULL, scale = FALSE, transposed = FALSE)
    )
    expect_error(
        runDPCoA(esophagus, assay_name = "counts", tree_name = "phylo", dimred = NULL, 
                 n_dimred = NULL, ncomponents = "test", ntop = "test",
                 subset_row = NULL, scale = FALSE, transposed = FALSE)
    )
    expect_error(
        runDPCoA(esophagus, assay_name = "counts", tree_name = "phylo", dimred = NULL, 
                 n_dimred = NULL, ncomponents = 1.3, ntop = "test",
                 subset_row = NULL, scale = FALSE, transposed = FALSE)
    )
    expect_error(
        runDPCoA(esophagus, name = c("test", "test2"), assay_name = "counts", 
                 tree_name = "phylo", dimred = NULL, 
                 n_dimred = NULL, ncomponents = 1.3, ntop = "test",
                 subset_row = NULL, scale = FALSE, transposed = FALSE)
    )
    expect_error(
        runDPCoA(esophagus, name = 1, assay_name = "counts", tree_name = "phylo", dimred = NULL, 
                 n_dimred = NULL, ncomponents = 1.3, ntop = "test",
                 subset_row = NULL, scale = FALSE, transposed = FALSE)
    )
})
