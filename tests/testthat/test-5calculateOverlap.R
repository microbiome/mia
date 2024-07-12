
context("getOverlap")

test_that("getOverlap", {
    
    # Get data
    data(esophagus, package="mia")
    tse <- esophagus
    # Test input
    expect_error(getOverlap(tse, assay.type = "relabundance", detection = 0.15))
    expect_error(getOverlap(tse, detection = "TEST"))
    expect_error(getOverlap(tse, detection = TRUE))
    
    # Calculate overlap
    tse <- transformAssay(tse, method = "relabundance")
    result <- getOverlap(tse, assay.type = "relabundance", detection = 0.15)

    # Test output
    expect_true(class(result) == "dist")
    # Test values
    result <- as.matrix(result)
    
    # Reference
    reference <- matrix(c(0.0000000, 0.3799961, 0.4111838, 0.3799961, 0.0000000,
                          0.3634972, 0.4111838, 0.3634972, 0.0000000), nrow=3)
    expect_equal(round(unname(result), 7), round(reference), 7)
    # Test with different detection threshold
    result <- getOverlap(tse, assay.type = "relabundance", detection = 0)
    result <- as.matrix(result)
    
    # Reference
    reference <- matrix(c(0.0000000, 0.8811552, 0.9038734, 0.8811552, 0.0000000, 
                          0.8390008, 0.9038734, 0.8390008, 0.00000000), nrow=3)
    expect_equal(round(unname(result), 7), round(reference, 7))
    
    # Test .calculate_overlap
    expect_error(addOverlap(tse, name = 1))
    expect_error(addOverlap(tse, name = c("A", "B")))
    expect_error(addOverlap(tse, name = TRUE))
    expect_equal( reducedDim(addOverlap(tse)), as.matrix(getOverlap(tse)) )
    expect_equal( reducedDim(addOverlap(tse, assay.type = "relabundance", detection = 0.08)), 
                  as.matrix(getOverlap(tse, assay.type = "relabundance", detection = 0.08)) )
})

