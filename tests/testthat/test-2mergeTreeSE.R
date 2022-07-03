context("mergeTreeSE")
test_that("mergeTreeSE", {
    # Load data
    data("GlobalPatterns")
    data("esophagus")
    data("enterotype")
    
    tse1 <- GlobalPatterns[1:50, ]
    tse2 <- esophagus[1:50, ]
    tse3 <- enterotype[1:50, ]
    
    # Expect errors
    expect_error( mergeTreeSE(tse1) )
    # expect_error( mergeTreeSE(tse1, tse2, missing_values = c(3, 3)) )
    expect_error( mergeTreeSE(tse1, tse2, missing_values = TRUE ) )
    expect_error( mergeTreeSE(tse1, tse2, missing_values = 36846 ) )
    expect_error( mergeTreeSE(tse1, tse2, abund_values = "test")  )
    # Calculate relative transform to test abund_values
    tse1 <- transformSamples(tse1, method = "relabundance")
    expect_error( mergeTreeSE(tse1, tse2, abund_values = "relabundance")  )
    expect_error( mergeTreeSE(tse1, tse2, verbose = "test")  )
    expect_error( mergeTreeSE(tse1, tse2, verbose = 1)  )
    expect_error( mergeTreeSE(tse1, tse2, tse3)  )
    expect_error( mergeTreeSE(tse1)  )
    
    # Test that data match if there is only one element
    tse <- mergeTreeSE(list(tse1), abund_values = "relabundance")
    expect_equal( rowData(tse), rowData(tse1))
    expect_equal( colData(tse), colData(tse1))
    expect_equal( assay(tse, "relabundance"), assay(tse1, "relabundance"))
    expect_equal( rowTree(tse), rowTree(tse1))
    
    # Test that data match if there is only same elements
    tse <- mergeTreeSE(list(tse1, tse1, tse1), abund_values = "relabundance")
    # The order of taxa and samples changes
    tse <- tse[ rownames(tse1), colnames(tse1) ]
    expect_equal( rowData(tse), rowData(tse1))
    expect_equal( colData(tse), colData(tse1))
    expect_equal( assay(tse, "relabundance"), assay(tse1, "relabundance"))
    expect_equal( rowTree(tse), rowTree(tse1))
    
    # Expect that rowTree is preserved if rownames match
    tse <- mergeTreeSE(list(tse1, GlobalPatterns), 
                                         abund_values = "counts",
                                         missing_values = NA)
    expect_equal(rowTree(GlobalPatterns), rowTree(tse))
    # Expect some NAs
    tse <- mergeTreeSE(list(tse1, tse2), 
                                         abund_values = "counts",
                                         missing_values = NA)
    expect_true( any(is.na(assay(tse))) )
    
    # Test that dimensions match
    tse <- mergeTreeSE(tse1, tse2)
    expect_equal( dim(tse), dim(tse1)+dim(tse2) )
    # Expect no NAs in assay
    expect_true( all(!is.na(assay(tse))) )
    
    # Test that dimensions match
    tse <- mergeTreeSE(list(tse1, tse2, tse3), missing_values = "MISSING")
    expect_equal( dim(tse), dim(tse1)+dim(tse2)+dim(tse3) )
    # Expect some "MISSING"s
    expect_true( any( assay(tse) == "MISSING" ) )
    
    # Test rownames
    expect_true( all(rownames(tse1) %in% rownames(tse)) )
    expect_true( all(rownames(tse2) %in% rownames(tse)) )
    expect_true( all(rownames(tse2) %in% rownames(tse)) )
    
    # Test colnames
    expect_true( all(colnames(tse1) %in% colnames(tse)) )
    expect_true( all(colnames(tse2) %in% colnames(tse)) )
    expect_true( all(colnames(tse3) %in% colnames(tse)) )
    
    tse <- mergeTreeSE(list(tse2, tse3, tse1, 
                                              tse1[1:2, ], tse1[1, ]), 
                                         missing_values = NA)
    # Get assay (as.matrix to remove links)
    assay <- assay(tse, "counts")
    assay1 <- as.matrix( assay(tse1, "counts") )
    assay2 <- as.matrix( assay(tse2, "counts") )
    assay3 <- as.matrix( assay(tse3, "counts") )
    
    # Expect that the data can be found without modifications
    colnames <- colnames(assay1)
    rownames <- rownames(assay1)
    expect_equal( assay[rownames, colnames], assay1 )
    colnames <- colnames(assay2)
    rownames <- rownames(assay2)
    expect_equal( assay[rownames, colnames], assay2 )
    colnames <- colnames(assay3)
    rownames <- rownames(assay3)
    expect_equal( assay[rownames, colnames], assay3 )
    
    # Get rowData (as.data.frame to remove links)
    row_data <- as.data.frame( rowData(tse) )
    row_data1 <- as.data.frame(rowData(tse1) )
    row_data2 <- as.data.frame(rowData(tse2) )
    row_data3 <- as.data.frame( rowData(tse3) )
    
    # Expect that the data can be found without modifications
    colnames <- colnames(row_data1)
    rownames <- rownames(row_data1)
    expect_equal( row_data[rownames, colnames], row_data1 )
    colnames <- colnames(row_data2)
    rownames <- rownames(row_data2)
    expect_equal( row_data[rownames, colnames], row_data2 )
    colnames <- colnames(row_data3)
    rownames <- rownames(row_data3)
    expect_equal( row_data[rownames, colnames, drop = FALSE], row_data3 )
    
    # Get colData (as.data.frame to remove links)
    col_data <- as.data.frame( colData(tse) )
    col_data1 <- as.data.frame( colData(tse1) )
    col_data2 <- as.data.frame( colData(tse2) )
    col_data3 <- as.data.frame( colData(tse3) )
    
    # Expect that the data can be found without modifications
    colnames <- colnames(col_data1)
    rownames <- rownames(col_data1)
    expect_equal( col_data[rownames, colnames], col_data1 )
    colnames <- colnames(col_data2)
    rownames <- rownames(col_data2)
    expect_equal( col_data[rownames, colnames], col_data2 )
    colnames <- colnames(col_data3)
    rownames <- rownames(col_data3)
    expect_equal( col_data[rownames, colnames], col_data3 )
    
    
    # TODO: test joining-method, all different joining methods, metadata
})
