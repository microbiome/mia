
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

test_that("Test MAE/TreeSE convertsion functions", {
    data("HintikkaXOData")
    mae1 <- HintikkaXOData

    tse1 <- convertToTreeSE(mae1)
    mae2 <- convertToMAE(tse1)
    tse2 <- convertToTreeSE(mae2)

    # Check that the data is retained correctly
    expect_equal(mae1, mae2, check.attributes = FALSE)
    expect_equal(colData(mae1), colData(mae2), check.attributes = FALSE)
    expect_equal(mae1[[1L]], mae2[[1L]], check.attributes = FALSE)
    expect_equal(tse1, tse2, check.attributes = FALSE)
    expect_equal(assay(tse1), assay(tse2), check.attributes = FALSE)
    expect_equal(rowData(tse1), rowData(tse2), check.attributes = FALSE)
    expect_equal(colData(tse1), colData(tse2), check.attributes = FALSE)
})
