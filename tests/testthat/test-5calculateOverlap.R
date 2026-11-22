
context("Overlap")

test_that("Overlap", {

    # Get data
    data(esophagus, package="mia")
    tse <- esophagus
    # Test input
    expect_error(
        getDissimilarity(tse, method = "overlap", assay.type = "relabundance", detection = 0.15))
    expect_error(getDissimilarity(tse, method = "overlap", detection = "TEST"))
    expect_error(getDissimilarity(tse, method = "overlap", detection = TRUE))

    # Calculate overlap
    tse <- transformAssay(tse, method = "relabundance")
    result <- getDissimilarity(tse, method = "overlap", assay.type = "relabundance", detection = 0.15)

    # Test output
    expect_true(class(result) == "dist")
    # Test values: dissimilarity is 1 - overlap
    result <- as.matrix(result)
    overlap <- matrix(c(1, 0.3799961, 0.4111838, 0.3799961, 1,
                        0.3634972, 0.4111838, 0.3634972, 1), nrow=3)
    expect_equal(unname(result), 1 - overlap, tolerance = 1e-6)
    # Test with different detection threshold
    result <- getDissimilarity(tse, method = "overlap", assay.type = "relabundance", detection = 0)
    result <- as.matrix(result)
    overlap <- matrix(c(1, 0.8811552, 0.9038734, 0.8811552, 1,
                        0.8390008, 0.9038734, 0.8390008, 1), nrow=3)
    expect_equal(unname(result), 1 - overlap, tolerance = 1e-6)
    # Identical samples have dissimilarity 0
    tse2 <- cbind(tse, tse[, 1])
    colnames(tse2)[4] <- "copy"
    result <- as.matrix(getDissimilarity(
        tse2, method = "overlap", assay.type = "relabundance"))
    expect_equal(result[1, "copy"], 0)
})
