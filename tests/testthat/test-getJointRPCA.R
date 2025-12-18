test_that("getJointRPCA stores embedding in reducedDim and metadata", {
    skip_if_not(requireNamespace("SingleCellExperiment", quietly = TRUE))
    
    set.seed(123)
    X <- abs(matrix(rnorm(20 * 10), nrow = 20, ncol = 10))
    rownames(X) <- paste0("f", 1:20)
    colnames(X) <- paste0("s", 1:10)
    
    sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = X)
    )
    
    sce2 <- mia:::getJointRPCA(
        sce,
        name = "JointRPCA_test",
        n.components = 2,
        max.iterations = 2,
        rclr.transform.tables = FALSE,
        n.test.samples = 3
    )
    
    #reducedDim populated
    emb <- reducedDim(sce2, "JointRPCA_test")
    expect_true(is.matrix(emb))
    expect_equal(nrow(emb), ncol(X))
    expect_equal(rownames(emb), colnames(X))
    
    #metadata populated
    jr <- metadata(sce2)$JointRPCA[["JointRPCA_test"]]
    expect_type(jr, "list")
    expect_true(all(c("ord.res", "dist", "cv.stats", "rclr.tables") %in% names(jr)))
})