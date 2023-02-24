
context("calculateOverlap")

test_that("calculateOverlap", {
    
    # Get data
    data("esophagus")
    tse <- esophagus
    # Test input
    expect_error(calculateOverlap(tse, assay_name = "relabundance", detection = 0.15))
    expect_error(calculateOverlap(tse, detection = "TEST"))
    expect_error(calculateOverlap(tse, detection = TRUE))
    
    # Calculate overlap
    tse <- transformCounts(tse, method = "relabundance")
    result <- calculateOverlap(tse, assay_name = "relabundance", detection = 0.15)
    # Test output
    expect_true(class(result) == "dist")
    # Test values
    result <- as.matrix(result)
    # Reference
    # data("esophagus", package = "phyloseq")
    # reference <- microbiome::overlap(esophagus, detection = 0.15)
    reference <- matrix(c(0.0000000, 0.3799961, 0.4111838, 0.3799961, 0.0000000,
                          0.3634972, 0.4111838, 0.3634972, 0.0000000), nrow=3)
    expect_equal(round(unname(result), 7), round(reference), 7)
    # Test with different detection threshold
    result <- calculateOverlap(tse, assay_name = "relabundance", detection = 0)
    result <- as.matrix(result)
    # Reference
    # data("esophagus", package = "phyloseq")
    # reference <- microbiome::overlap(esophagus, detection = 0)
    reference <- matrix(c(0.0000000, 0.8811552, 0.9038734, 0.8811552, 0.0000000, 
                          0.8390008, 0.9038734, 0.8390008, 0.00000000), nrow=3)
    expect_equal(round(unname(result), 7), round(reference, 7))
    # Test runOverlap
    expect_error(runOverlap(tse, name = 1))
    expect_error(runOverlap(tse, name = c("A", "B")))
    expect_error(runOverlap(tse, name = TRUE))
    expect_equal( reducedDim(runOverlap(tse)), as.matrix(calculateOverlap(tse)) )
    expect_equal( reducedDim(runOverlap(tse, assay_name = "relabundance", detection = 0.08)), 
                  as.matrix(calculateOverlap(tse, assay_name = "relabundance", detection = 0.08)) )
})

