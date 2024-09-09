
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
    # Test values
    result <- as.matrix(result)
    
    # Reference
    reference <- matrix(c(0.0000000, 0.3799961, 0.4111838, 0.3799961, 0.0000000,
                          0.3634972, 0.4111838, 0.3634972, 0.0000000), nrow=3)
    expect_equal(round(unname(result), 7), round(reference), 7)
    # Test with different detection threshold
    result <- getDissimilarity(tse, method = "overlap", assay.type = "relabundance", detection = 0)
    result <- as.matrix(result)
    
    # Reference
    reference <- matrix(c(0.0000000, 0.8811552, 0.9038734, 0.8811552, 0.0000000, 
                          0.8390008, 0.9038734, 0.8390008, 0.00000000), nrow=3)
    expect_equal(round(unname(result), 7), round(reference, 7))
    
    # Test .calculate_overlap
    expect_equal( metadata(addDissimilarity(tse, "overlap"))[["overlap"]], 
                  as.matrix(getDissimilarity(tse, method = "overlap")) )
    expect_equal( metadata(addDissimilarity(tse, method = "overlap", assay.type = "relabundance", detection = 0.08))[["overlap"]], 
                  as.matrix(getDissimilarity(tse, method = "overlap", assay.type = "relabundance", detection = 0.08)) )
})

