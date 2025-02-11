
test_that("Test getReducedDimElement", {
    data(GlobalPatterns)
    tse <- GlobalPatterns
    #
    res <- getReducedDimElement(tse) |> expect_error()
    # Calculate MDS
    tse <- runMDS(tse, assay.type = "counts")
    # Get reference values
    mat <- reducedDim(tse, "MDS")
    ref <- attributes(mat)
    ref <- ref[ !names(ref) %in% c("dim", "dimnames") ]
    # Check correct errors
    res <- getReducedDimElement(tse, "test") |> expect_error()
    res <- getReducedDimElement(tse, TRUE) |> expect_error()
    res <- getReducedDimElement(tse, c("test", "test2")) |> expect_error()
    # Test that values are correct
    res <- getReducedDimElement(tse)
    expect_equal(res, ref)
    #
    for( nam in names(ref) ){
        # Get randomly an index to test
        if( rnorm(1)>0 ){
            nam <- which(nam == names(ref))
        }
        res <- getReducedDimElement(tse, dimred = "MDS", name = nam)
        expect_equal(res, ref[[nam]])
    }
})
