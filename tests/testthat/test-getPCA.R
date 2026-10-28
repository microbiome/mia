test_that("getPCA and addPCA", {
    data(GlobalPatterns, package = "mia")
    # The signs of the components from scater::calculatePCA() can differ
    # between identical calls, so components are compared up to their sign
    expect_same_pca <- function(res, ref){
        expect_equal(unname(abs(diag(cor(res, ref)))), rep(1, ncol(ref)))
        expect_equal(attr(res, "percentVar"), attr(ref, "percentVar"))
    }
    tse <- transformAssay(GlobalPatterns, method = "clr", pseudocount = TRUE)
    # All features are used by default
    res <- getPCA(tse, assay.type = "clr", ncomponents = 5)
    ref <- scater::calculatePCA(
        tse, assay.type = "clr", ncomponents = 5, ntop = Inf)
    expect_same_pca(res, ref)
    # ntop is passed to scater
    res <- getPCA(tse, assay.type = "clr", ncomponents = 5, ntop = 500)
    ref <- scater::calculatePCA(
        tse, assay.type = "clr", ncomponents = 5, ntop = 500)
    expect_same_pca(res, ref)
    # Results are added to reducedDim
    tse <- addPCA(tse, assay.type = "clr", ncomponents = 5, name = "test")
    expect_true("test" %in% reducedDimNames(tse))
    expect_equal(dim(reducedDim(tse, "test")), c(ncol(tse), 5L))
    # Input check
    expect_error(getPCA(tse, assay.type = "test"))
    expect_error(addPCA(tse, assay.type = "clr", name = c("a", "b")))
})
