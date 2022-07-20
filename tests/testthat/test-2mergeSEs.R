context("mergeSEs")
test_that("mergeSEs", {
    # Load data
    data("GlobalPatterns")
    data("esophagus")
    data("enterotype")
    
    tse1 <- GlobalPatterns[1:50, ]
    tse2 <- esophagus[1:50, ]
    tse3 <- enterotype[1:50, ]
    
    # Expect errors
    expect_error( mergeSEs(tse1) )
    expect_error( mergeSEs(tse1, tse2, join = 1) )
    expect_error( mergeSEs(tse1, tse2, join = TRUE) )
    expect_error( mergeSEs(tse1, tse2, join = NA) )
    expect_error( mergeSEs(tse1, tse2, collapse_samples = NA) )
    expect_error( mergeSEs(tse1, tse2, collapse_samples = 1) )
    expect_error( mergeSEs(tse1, tse2, collapse_samples = "test") )
    expect_error( mergeSEs(tse1, tse2, collapse_samples = NULL) )
    expect_error( mergeSEs(list(tse1, tse2, tse), join = "left") )
    expect_error( mergeSEs(list(tse1, tse2, tse), join = "right") )
    expect_error( mergeSEs(tse1, tse2, missing_values = TRUE ) )
    expect_error( mergeSEs(tse1, tse2, missing_values = 36846 ) )
    expect_error( mergeSEs(tse1, tse2, assay_name = "test")  )
    # Calculate relative transform to test assay_name
    tse1 <- transformSamples(tse1, method = "relabundance")
    expect_error( mergeSEs(tse1, tse2, assay_name = "relabundance")  )
    expect_error( mergeSEs(tse1, tse2, verbose = "test")  )
    expect_error( mergeSEs(tse1, tse2, verbose = 1)  )
    expect_error( mergeSEs(tse1, tse2, tse3)  )
    expect_error( mergeSEs(tse1)  )
    
    # Test that data match if there is only one element
    tse <- mergeSEs(list(tse1), assay_name = "relabundance")
    expect_equal( rowData(tse), rowData(tse1))
    expect_equal( colData(tse), colData(tse1))
    expect_equal( assay(tse, "relabundance"), assay(tse1, "relabundance"))
    expect_equal( rowTree(tse), rowTree(tse1))
    
    # Test that data match if there is only same elements
    tse <- mergeSEs(list(tse1, tse1, tse1), assay_name = "relabundance")
    # The order of taxa and samples changes
    tse <- tse[ rownames(tse1), colnames(tse1) ]
    expect_equal( rowData(tse), rowData(tse1))
    expect_equal( as.data.frame( colData(tse) ), as.data.frame( colData(tse1)) )
    expect_equal( assay(tse, "relabundance"), assay(tse1, "relabundance"))
    expect_equal( rowTree(tse), rowTree(tse1))
    
    # Expect that rowTree is preserved if rownames match
    tse <- mergeSEs(list(tse1, GlobalPatterns), 
                                         assay_name = "counts",
                                         missing_values = NA)
    expect_equal(rowTree(GlobalPatterns), rowTree(tse))
    # Expect some NAs
    tse <- mergeSEs(list(tse1, tse2), assay_name = "counts")
    expect_true( any(is.na(assay(tse))) )
    
    # Test that dimensions match
    tse <- mergeSEs(tse1, tse2, missing_values = 0)
    expect_equal( dim(tse), dim(tse1)+dim(tse2) )
    # Expect no NAs in assay
    expect_true( all(!is.na(assay(tse))) )
    
    # Test that dimensions match
    tse <- mergeSEs(list(tse1, tse2, tse3), missing_values = "MISSING")
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
    
    # CHECK FULL JOIN ###################################################
    tse <- mergeSEs(list(tse2, tse3, tse1, tse1[1:2, ], tse1[1, ]), 
                    missing_values = NA)
    # Get assay (as.matrix to remove links)
    assay <- as.matrix( assay(tse, "counts") )
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
    
    # CHECK INNER JOIN ##############################################
    tse <- mergeSEs(list(tse1[, 1:5], tse1[, 5:10], tse3[1:20, 6:10]), 
                       join = "inner")
    expect_true( nrow(tse) == 0 )
    expect_equal( rowTree(tse), NULL )
    tse <- mergeSEs(list(tse1[, 1:5], tse1[, 5:10], tse1[1:20, 6:10]), 
                       join = "inner", collapse_samples = TRUE)
    expect_true( all(dim(tse) == c(20, 10)) )
    expect_equal( rowTree(tse), rowTree(tse1) )
    # Get assay (as.matrix to remove links)
    assay <- as.matrix( assay(tse, "counts") )
    assay1 <- as.matrix( assay(tse1, "counts") )
    
    # Expect that the data can be found without modifications
    colnames <- colnames(assay)
    rownames <- rownames(assay)
    expect_equal( assay1[rownames, colnames], assay )
    
    # Get rowData (as.data.frame to remove links)
    row_data <- as.data.frame( rowData(tse) )
    row_data1 <- as.data.frame(rowData(tse1) )
    
    # Expect that the data can be found without modifications
    colnames <- colnames(row_data)
    rownames <- rownames(row_data)
    expect_equal( row_data1[rownames, colnames], row_data )
    
    # Get colData (as.data.frame to remove links)
    col_data <- as.data.frame( colData(tse) )
    col_data1 <- as.data.frame( colData(tse1) )
    
    # Expect that the data can be found without modifications
    colnames <- colnames(col_data)
    rownames <- rownames(col_data)
    expect_equal( col_data1[rownames, colnames], col_data )
    
    # CHECK LEFT JOIN ##############################################
    tse <- mergeSEs(list(tse1[11:20, 1:13], tse1[10:50, 7:20]), 
                       join = "left", collapse_samples = TRUE)
    expect_true( all(dim(tse) == c(10, 20)) )
    expect_equal( rowTree(tse), rowTree(tse1) )
    # Get assay (as.matrix to remove links)
    assay <- as.matrix( assay(tse, "counts") )
    assay1 <- as.matrix( assay(tse1, "counts") )
    
    # Expect that the data can be found without modifications
    colnames <- colnames(assay)
    rownames <- rownames(assay)
    expect_equal( assay1[rownames, colnames], assay )
    
    # Get rowData (as.data.frame to remove links)
    row_data <- as.data.frame( rowData(tse) )
    row_data1 <- as.data.frame(rowData(tse1) )
    
    # Expect that the data can be found without modifications
    colnames <- colnames(row_data)
    rownames <- rownames(row_data)
    expect_equal( row_data1[rownames, colnames], row_data )
    
    # Get colData (as.data.frame to remove links)
    col_data <- as.data.frame( colData(tse) )
    col_data1 <- as.data.frame( colData(tse1) )
    
    # Expect that the data can be found without modifications
    colnames <- colnames(col_data)
    rownames <- rownames(col_data)
    expect_equal( col_data1[rownames, colnames], col_data )
    
    # CHECK RIGHT JOIN ##############################################
    tse <- mergeSEs(list(tse1[10:50, 1:13], tse1[1:10, 7:20]), 
                    join = "right", missing_values = NA, 
                    collapse_samples = TRUE)
    expect_true( all(dim(tse) == c(10, 20)) )
    expect_equal( rowTree(tse), rowTree(tse1) )
    # Get assay (as.matrix to remove links)
    assay <- as.matrix( assay(tse, "counts") )
    assay1 <- as.matrix( assay(tse1, "counts") )
    
    # Expect that the data can be found without modifications
    colnames <- colnames(assay)
    rownames <- rownames(assay)
    assay1 <- assay1[rownames, colnames]
    # Add NAs to original
    na <- is.na(assay)
    assay1[na] <- NA
    expect_equal( assay1, assay )
    
    # Get rowData (as.data.frame to remove links)
    row_data <- as.data.frame( rowData(tse) )
    row_data1 <- as.data.frame(rowData(tse1) )
    
    # Expect that the data can be found without modifications
    colnames <- colnames(row_data)
    rownames <- rownames(row_data)
    expect_equal( row_data1[rownames, colnames], row_data )
    
    # Get colData (as.data.frame to remove links)
    col_data <- as.data.frame( colData(tse) )
    col_data1 <- as.data.frame( colData(tse1) )
    
    # Expect that the data can be found without modifications
    colnames <- colnames(col_data)
    rownames <- rownames(col_data)
    expect_equal( col_data1[rownames, colnames], col_data )
    
    
    # Check metadata
    metadata(tse1) <- list(abc = c("abc", 123))
    metadata(tse3) <- list(test = 1)
    metadata(tse) <- list( cd = colData(tse) )
    tse4 <- mergeSEs(list(tse, tse3, tse2, tse1), 
                        join = "inner")
    expect_equal( nrow(tse4), 0 )
    expect_equal( metadata(tse4)[["abc"]], metadata(tse1)[["abc"]] )
    expect_equal( metadata(tse4)[["test"]], metadata(tse3)[["test"]] )
    expect_equal( metadata(tse4)[["cd"]], metadata(tse)[["cd"]] )
    expect_true( all( names(metadata(tse4)) %in% 
                           c(names(metadata(tse1)), names(metadata(tse3)), 
                             names(metadata(tse))) ) )
    expect_equal( length( names(metadata(tse4))), 3) 
    
    # Check correct class
    tse <- mergeSEs(list(as(tse1, "SummarizedExperiment"), 
                            as(tse1, "SummarizedExperiment"),
                            as(tse1, "SummarizedExperiment")), 
                       join = "full")
    expect_true(class(tse) == "SummarizedExperiment")
    suppressWarnings(
    tse <- mergeSEs(list(as(tse1, "SummarizedExperiment"), 
                            as(tse1, "SingleCellExperiment"),
                            as(tse1, "TreeSummarizedExperiment")), 
                       join = "inner")
    )
    expect_true(class(tse) == "SummarizedExperiment")
    suppressWarnings(
    tse <- mergeSEs(list(as(tse1, "SummarizedExperiment"), 
                            as(tse1, "SingleCellExperiment"),
                            as(tse1, "SingleCellExperiment")), 
                       join = "full")
    )
    expect_true(class(tse) == "SummarizedExperiment")
    suppressWarnings(
    tse <- mergeSEs(x = as(tse1, "TreeSummarizedExperiment"), 
                       y = as(tse1, "SingleCellExperiment"), 
                       join = "right")
    )
    expect_warning(mergeSEs(x = as(tse1, "TreeSummarizedExperiment"), 
                            y = as(tse1, "SingleCellExperiment"), 
                            join = "right"))
    expect_true(class(tse) == "SingleCellExperiment")
    tse <- mergeSEs(list(as(tse1, "TreeSummarizedExperiment")), 
                       join = "left")
    expect_true(class(tse) == "TreeSummarizedExperiment")
    
    # Test dplyr-like aliases
    tse_test1 <- mergeSEs(x = tse[1:28, 1:3], 
                         y = tse1[23, 1:5], 
                         join = "full")
    tse_test2 <- full_join(x = tse[1:28, 1:3], 
                           y = tse1[23, 1:5])
    expect_equal(tse_test1, tse_test2)
    tse_test1 <- mergeSEs(x = tse[1:28, 1:3], 
                         y = tse1[23, 1:5], 
                         join = "left")
    tse_test2 <- left_join(x = tse[1:28, 1:3], 
                           y = tse1[23, 1:5])
    expect_equal(tse_test1, tse_test2)
    tse_test1 <- mergeSEs(x = list(tse1[1:28, 1:3], tse1[23, 1:5], tse1[2:4, ]), 
                         join = "inner")
    tse_test2 <- inner_join(x = list(tse1[1:28, 1:3], tse1[23, 1:5], tse1[2:4, ]) )
    expect_equal(tse_test1, tse_test2)
    tse_test1 <- mergeSEs(x = list(tse1[1:28, 1:3], tse1[23, 1:5]), 
                         join = "right")
    tse_test2 <- right_join(x = list(tse1[1:28, 1:3], tse1[23, 1:5]) )
    expect_equal(tse_test1, tse_test2)
    
    # Test collapse_samples
    tse_test <- mergeSEs(x = tse[1:28, 1:3], 
                         y = tse[23, 1:5], 
                         join = "full")
    expect_equal( dim(tse_test), c(28, 8))
    tse_test <- mergeSEs(x = list(tse[1:28, 1:3], tse[23, 1:5], tse[1, 1:10]),
                         join = "full")
    expect_equal( dim(tse_test), c(28, 18))
    expect_true( (all( c( paste0(rep(colnames(tse[, 1:3]), each=3), c("", "_2", "_3")), 
                         paste0(rep(colnames(tse[, 4:5]), each=2), c("", "_3")), 
                         colnames(tse[, 6:10]) ) %in% 
                          colnames(tse_test) ) &&
                     all( colnames(tse_test) %in% 
                     c( paste0(rep(colnames(tse[, 1:3]), each=3), c("", "_2", "_3")), 
                        paste0(rep(colnames(tse[, 4:5]), each=2), c("", "_3")), 
                        colnames(tse[, 6:10]) ) ) )
                 )
    # Test that tree is added after agglomeration
    agg_tse1 <- suppressWarnings( aggTSE(tse1, rowLevel = c(6,4,2)) )
    tse <- mergeSEs(tse1, agg_tse1)
    expect_equal(rowTree(tse), rowTree(tse1))
})
