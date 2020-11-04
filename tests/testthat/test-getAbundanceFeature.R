context("getAbundanceFeature")
test_that("getAbundanceFeature", {

    data(GlobalPatterns, package = "MicrobiomeExperiment")

    feature_ab <- getAbundanceFeature(GlobalPatterns,
                                      feature_id = "522457",
                                      abund_values = "counts")

    expect_equal(names(feature_ab)[1], "CL3")

    expect_error(
        mia:::.check_feature_ids_assays(GlobalPatterns,
                                        feature_id="x522457",
                                        abund_values="counts"),
        "Please provide a valid 'feature_id'")

})
