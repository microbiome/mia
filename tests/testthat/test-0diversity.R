
context("diversity estimates")
test_that("diversity estimates", {

    skip_if_not(requireNamespace("vegan", quietly = TRUE))
    data("esophagus")

    tse <- esophagus
    
    tse <- estimateDiversity(tse, threshold = 0.473)
    cd <- colData(tse)
    expect_equal(unname(round(cd$shannon, 5)), c(2.24937, 2.76239, 2.03249))
    expect_equal(unname(round(cd$simpson_diversity, 6)),
                 c(0.831372, 0.903345, 0.665749))
    expect_equal(unname(round(cd$inv_simpson, 5)),
                 c(5.93021, 10.34606, 2.99177))
    expect_equal(unname(round(cd$coverage, 0)), c(2,3,1))
    expect_equal(estimateDiversity(esophagus,
                                   index = "coverage",
                                   threshold = 0.9),
                 estimateCoverage(esophagus))
})
