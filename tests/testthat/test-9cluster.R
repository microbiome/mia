context("rarefyAssay")
test_that("rarefyAssay", {
    library(bluster)
    data(GlobalPatterns, package="mia")
    
    # Parameters validity check
    expect_error(addCluster(GlobalPatterns, 
                        KmeansParam(centers = 3), 
                        assay.type = "error"))
    # Checking wrong by (char)
    expect_error(addCluster(GlobalPatterns, 
                        KmeansParam(centers = 3), 
                        by = "error"))
    # Checking wrong by (number)
    expect_error(addCluster(GlobalPatterns, 
                         KmeansParam(centers = 3), 
                         by = 3))
    tse <- addCluster(GlobalPatterns, 
                   KmeansParam(centers = 3), 
                   name = "custommetadata",
                   full = TRUE,
                   clust.col = "customdataname")
    altExp(tse, "test") <- tse[1:1000,]
    # Checking same name that is already present
    expect_warning(addCluster(tse, 
                         KmeansParam(centers = 3), 
                         name = "custommetadata",
                         full = TRUE))
    # Checking wrong clust.col with already-present name
    expect_warning(addCluster(tse, 
                         KmeansParam(centers = 3), 
                         clust.col = "customdataname"))
    # Checking wrong altexp
    expect_error(addCluster(tse, 
                        KmeansParam(centers = 3), 
                        altexp = "error"))
    
    # Parameters check
    tse <- GlobalPatterns
    altExp(tse, "test") <- tse[1:1000,]
    tse <- addCluster(tse, 
                   KmeansParam(centers = 3), 
                   name = "custommetadata",
                   full = TRUE,
                   clust.col = "customdataname")
    tse <- addCluster(tse, 
                   KmeansParam(centers = 3), 
                   name = "custommetadata",
                   full = TRUE,
                   clust.col = "customdataname",
                   altexp = "test")
    # Checking custom metadata/dataname in main/altExp
    expect_true("custommetadata" %in% names(metadata(tse)))
    expect_true("customdataname" %in% names(rowData(tse)))
    expect_true("custommetadata" %in% names(metadata(altExp(tse, "test"))))
    expect_true("customdataname" %in% names(rowData(altExp(tse, "test"))))
    
    # Checking existing custom metadata/dataname in main/altExp
    expect_warning(addCluster(tse, 
                         KmeansParam(centers = 3), 
                         name = "custommetadata",
                         full = TRUE))
    expect_warning(addCluster(tse, 
                         KmeansParam(centers = 3), 
                         clust.col = "customdataname"))
    expect_warning(addCluster(tse, 
                         KmeansParam(centers = 3), 
                         name = "custommetadata",
                         full = TRUE,
                         clust.col = "customdataname",
                         altexp = "test"))
    # Checking working by
    tse <- GlobalPatterns
    altExp(tse, "test") <- tse[1:1000,]
    tse <- addCluster(tse, 
                    KmeansParam(centers = 3),
                    by = "col")
    tse <- addCluster(tse, 
                    KmeansParam(centers = 3),
                    altexp = "test",
                    by = 2)
    expect_true("clusters" %in% names(colData(tse)))
    expect_true("clusters" %in% names(colData(altExp(tse, "test"))))
    
    # Checking wrapper operational
    tse <- GlobalPatterns
    altExp(tse, "test") <- tse[1:2000,]
    tse <- addCluster(tse, 
                   HclustParam(),
                   by = "col")
    tse <- addCluster(tse, 
                   HclustParam(),
                   by = "row",
                   altexp = "test",
                   full = TRUE)
    expectedCol <- clusterRows(t(assay(tse, "counts")), HclustParam())
    expectedRow <- clusterRows(assay(altExp(tse, "test"), "counts"), 
                               HclustParam(), full = TRUE)
    # Checking same output on cols
    expect_identical(expectedCol, colData(tse)$clusters)
    # Checking same output on rows
    expect_identical(expectedRow$clusters, 
                    rowData(altExp(tse, "test"))$clusters)
    # Checking same metdata output on rows
    expect_identical(expectedRow$objects, 
                    metadata(altExp(tse, "test"))$clusters)
})
