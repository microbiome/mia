context("distance")
test_that("distance", {
    library(SEtup)
    data("enterotype", package = "MicrobiomeExperiment")
    data("esophagus", package = "MicrobiomeExperiment")
    # default
    actual <- calculateDistance(enterotype)
    expect_equal(actual[1], 0.22382452, tolerance = .0000001)
    expect_equal(actual[length(actual)], 0.3316593, tolerance = .0000001)
    # JSD
    actual <- calculateJSD(enterotype)
    expect_equal(actual[1], 0.1282428, tolerance = .0000001)
    expect_equal(actual[length(actual)], 0.06595946, tolerance = .0000001)
    # UniFrac
    actual <- calculateUniFrac(esophagus, weighted = FALSE, normalized = TRUE)
    expect_equal(actual[1], 0.5175550, tolerance = .0000001)
    expect_equal(actual[length(actual)], 0.5422394, tolerance = .0000001)
    actual <- calculateUniFrac(esophagus, weighted = TRUE, normalized = TRUE)
    expect_equal(actual[1], 0.2035424, tolerance = .0000001)
    expect_equal(actual[length(actual)], 0.2477016, tolerance = .0000001)
    actual <- calculateUniFrac(esophagus, weighted = TRUE, normalized = FALSE)
    expect_equal(actual[1], 0.1050480, tolerance = .0000001)
    expect_equal(actual[length(actual)], 0.1422409, tolerance = .0000001)
    #
    enterotype <- runMDS2(enterotype, FUN = calculateJSD, name = "JSD")
    expect_named(reducedDims(enterotype),"JSD")
    esophagus <- runMDS2(esophagus, FUN = calculateUniFrac, name = "UniFrac",
                         tree = rowTree(esophagus))
    expect_named(reducedDims(esophagus),"UniFrac")
})
