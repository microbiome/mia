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
})
