context("DMN")
test_that("DMN", {
    data(dmn_se)
    actual <- expect_warning(calculateDMN(dmn_se))
    expect_true(is.list(actual))
    expect_s4_class(actual[[1L]],"DMN")
    # return a list of DMN objects
    actual <- expect_warning(getDMN(dmn_se))
    expect_true(is.list(actual))
    expect_s4_class(actual[[1L]],"DMN")
    # return, which objects fits best
    expect_equal(expect_warning(bestDMNFit(dmn_se, type = "laplace")), 4)
    # return the model, which fits best
    actual <- expect_warning(getBestDMNFit(dmn_se, type = "laplace"))
    expect_s4_class(actual,"DMN")
})
