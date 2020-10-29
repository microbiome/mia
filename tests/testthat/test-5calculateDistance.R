context("calculateDistance")
test_that("calculateDistance", {
    mat <- matrix(1:60, nrow = 6)
    df <- DataFrame(n = c(1:6))
    expect_error(SEtup:::.merge_rows(),
                 'argument "f" is missing')
    se <- SummarizedExperiment(assays = list(counts = mat),
                               rowData = df)
    # default
    actual <- calculateDistance(se)
    actual2 <- dist(t(assay(se,"counts")))
    expect_true(all(actual == actual2))
    #
    data("enterotype")
    data("esophagus")
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
