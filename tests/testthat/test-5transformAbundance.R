context("transformAbundance")

test_that("transformAbundance", {

    # TSE object
    data(GlobalPatterns)
    tse <- GlobalPatterns

    # No method specified. Should be an error.
    testthat::expect_error(mia::transformAbundance(tse))

    # Method is not provided. Should be an error.
    testthat::expect_error(mia::transformAbundance(tse, method="test"))

    # Counts table should not be changed
    testthat::expect_equal(assays(mia::transformAbundance(tse, method = "pa"))$counts, assays(tse)$counts)

    # Calculates relative abundances. Should be equal.
    testthat::expect_equal(assays(mia::transformAbundance(tse, method = "relabundance"))$relabundance,
                           apply(assays(tse)$counts, 2, FUN=function(x){
                               x/sum(x)
                           }))

    # Calculates log10 transformation with pseudocount = 1. Should be equal.
    testthat::expect_equal(assays(mia::transformAbundance(tse, method = "log10", pseudocount = 1))$log10,
                           apply(assays(tse)$counts, 2, FUN=function(x){
                               log10(x+1)
                           }))

    # Calculates pa transformation. Should be equal.
    testthat::expect_equal(as.vector(assays(mia::transformAbundance(tse, method = "pa"))$pa),
                           as.integer(assays(tse)$counts > 0))


})

