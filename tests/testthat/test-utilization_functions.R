
test_that("Test getReducedDimAttribute", {
    data(GlobalPatterns)
    tse <- GlobalPatterns
    #
    res <- getReducedDimAttribute(tse) |> expect_error()
    # Calculate MDS
    tse <- runMDS(tse, assay.type = "counts")
    # Get reference values
    mat <- reducedDim(tse, "MDS")
    ref <- attributes(mat)
    ref <- ref[ !names(ref) %in% c("dim", "dimnames") ]
    # Check correct errors
    res <- getReducedDimAttribute(tse, "test") |> expect_error()
    res <- getReducedDimAttribute(tse, TRUE) |> expect_error()
    res <- getReducedDimAttribute(tse, c("test", "test2")) |> expect_error()
    # Test that values are correct
    res <- getReducedDimAttribute(tse)
    expect_equal(res, ref)
    #
    for( nam in names(ref) ){
        # Get randomly an index to test
        if( rnorm(1)>0 ){
            nam <- which(nam == names(ref))
        }
        res <- getReducedDimAttribute(tse, dimred = "MDS", name = nam)
        expect_equal(res, ref[[nam]])
    }
})
