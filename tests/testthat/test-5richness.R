
context("richness estimates")
test_that("richness estimates", {
    skip_if_not(requireNamespace("vegan", quietly = TRUE))
    data("esophagus")
    
    esophagus <- estimateRichness(esophagus)
    cd <- colData(esophagus)
    expect_equal(unname(round(cd$observed, 0)), c(28, 33, 38))

    esophagus <- estimateRichness(esophagus, detection = 0)
    cd <- colData(esophagus)
    expect_equal(unname(round(cd$observed, 0)), c(28, 33, 38))

    esophagus <- estimateRichness(esophagus, detection = 1)
    cd <- colData(esophagus)
    expect_equal(unname(round(cd$observed, 0)), c(15, 24, 16))

    # These are unaffected by detection parameter
    expect_equal(unname(round(cd$Chao1, 4)), c(39.1429, 37.5000, 71.0000))
    expect_equal(unname(round(cd$ACE, 4)), c(49.0970, 40.9465, 88.9768))
    expect_equal(unname(round(cd$Hill, 4)), c(9.4817, 15.8376, 7.6331))

})

