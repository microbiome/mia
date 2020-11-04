context("getAbundanceSample")
test_that("getAbundanceSample", {

    data(GlobalPatterns, package = "MicrobiomeExperiment")

    sam_ab <- getAbundanceSample(GlobalPatterns,
                                 sample_id = "CC1",
                                 abund_values = "counts")

    expect_equal(names(sam_ab)[1], "549322")

    expect_error(
        mia:::.check_sample_ids_assays(GlobalPatterns,
                                       sample_id= "bogus",
                                       abund_values="counts"),
        "Please provide a valid 'sample_id'")

})
