context("DMN")
test_that("DMN", {
    data(dmn_se)
    actual <- calculateDMN(dmn_se)
    expect_true(is.list(actual))
    expect_s4_class(actual[[1L]],"DMN")
    # return a list of DMN objects
    actual <- getDMN(dmn_se)
    expect_true(is.list(actual))
    expect_s4_class(actual[[1L]],"DMN")
    # return, which objects fits best
    expect_equal(bestDMNFit(dmn_se, type = "laplace"), 4)
    # return the model, which fits best
    actual <- getBestDMNFit(dmn_se, type = "laplace")
    expect_s4_class(actual,"DMN")
})
