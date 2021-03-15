
context("estimateRichness")

test_that("estimateRichness", {

    skip_if_not(requireNamespace("vegan", quietly = TRUE))
    data("esophagus")

    tse <- estimateRichness(esophagus, detection = 1)
    cd <- colData(tse)    
    expect_equal(unname(round(cd$observed, 0)), c(15, 24, 16))
    # These are unaffected by detection parameter
    expect_equal(unname(round(cd$Chao1, 4)), c(39.1429, 37.5000, 71.0000))
    expect_equal(unname(round(cd$ACE, 4)), c(49.0970, 40.9465, 88.9768))
    expect_equal(unname(round(cd$Hill, 4)), c(9.4817, 15.8376, 7.6331))

    test_internal_estimateRichness <- function(tse){

        # Calculate all indices.
        tse_idx <- estimateRichness(tse)

        # Check that the type of output is the same as the type of input.
        expect_true(typeof(tse_idx) == typeof(tse))

        # Check that every index is calculated by checking the column names from
        # colData.
        # Check that the order of indices is right / the same as the order
        # in the input vector.
        expect_named(colData(tse_idx), c("observed", "Chao1", "Chao1_se", "ACE", "ACE_se", "Hill"))

        # .get_observed
        mat <- assay(tse, "counts")
        expect_equal(unname(mia:::.calc_observed(mat, detection = 0)), c(28, 33, 38))
        expect_equal(unname(mia:::.calc_observed(mat, detection = 1)), c(15, 24, 16))	

        s <- mia:::.calc_chao1(mat)
        expect_equal(ncol(s), 2)
        expect_equal(unname(round(s[,1], 4)), c(39.1429, 37.5000, 71.0000))

        s <- mia:::.calc_ACE(mat)
        expect_equal(ncol(s), 2)
        expect_equal(unname(round(s[,1], 4)), c(49.0970, 40.9465, 88.9768))

        s <- mia:::.calc_hill(mat)
	expect_false(is.matrix(s))
        expect_equal(unname(round(s, 4)), c(9.4817, 15.8376, 7.6331))

    }

    # TSE object
    tse <- esophagus

    # Standard tse
    test_internal_estimateRichness(tse)
    
    # DelayedArray version of the assay
    assay(tse,"counts") <- DelayedArray(assay(tse,"counts"))
    test_internal_estimateRichness(tse)

})

