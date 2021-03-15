context("diversity estimates")
test_that("diversity estimates", {

    skip_if_not(requireNamespace("vegan", quietly = TRUE))
    data("esophagus")

    tse <- esophagus
    tse <- relAbundanceCounts(tse)
    
    tse_idx <- estimateDiversity(tse, threshold = 0.473)

    # Checks that the type of output is the same as the type of input.
    expect_true(typeof(tse_idx) == typeof(tse))

    # Check that every index is calculated by checking the column names from
    # colData.
    # Check that the order of indices is right / the same as the order
    # in the input vector.
    expect_named(colData(tse_idx), c("shannon","gini_simpson","inverse_simpson", "coverage", "fisher"))

    lambda <- unname(colSums(assays(tse_idx)$relabundance^2))
    ginisimpson <- 1 - lambda
    invsimpson <- 1 / lambda    

    expect_equal(lambda, .simpson_lambda(assays(tse_idx)$relabundance))
    expect_equal(ginisimpson, colData(tse_idx)$gini_simpson)
    expect_equal(ginisimpson, .calc_gini_simpson(assays(tse_idx)$relabundance))

    expect_equal(invsimpson, colData(tse_idx)$inverse_simpson)    
    expect_equal(invsimpson, .calc_inverse_simpson(assays(tse_idx)$relabundance))
    
    cd <- colData(tse_idx)
    expect_equal(unname(round(cd$shannon, 5)), c(2.24937, 2.76239, 2.03249))
    expect_equal(unname(round(cd$gini_simpson, 6)),
                 c(0.831372, 0.903345, 0.665749))
    expect_equal(unname(round(cd$inverse_simpson, 5)),
                 c(5.93021, 10.34606, 2.99177))
    expect_equal(unname(round(cd$coverage, 0)), c(2,3,1))
    expect_equal(unname(round(cd$fisher, 4)), c(8.8037, 10.0989, 13.2783))    
    expect_equal(estimateDiversity(esophagus,
                                   index = "coverage",
                                   threshold = 0.9),
                 estimateCoverage(esophagus))
})
