context("addCluster")
test_that("addCluster", {
    library(bluster)
    data(GlobalPatterns, package = "mia")

    # Parameters validity check
    expect_error(addCluster(GlobalPatterns,
        KmeansParam(centers = 3),
        assay.type = "error"
    ))
    # Checking wrong by (char)
    expect_error(addCluster(GlobalPatterns,
        KmeansParam(centers = 3),
        by = "error"
    ))
    # Checking wrong by (number)
    expect_error(addCluster(GlobalPatterns,
        KmeansParam(centers = 3),
        by = 3
    ))
    tse <- addCluster(GlobalPatterns,
        KmeansParam(centers = 3),
        name = "custommetadata",
        assay.type = "counts",
        full = TRUE,
        clust.col = "customdataname"
    )
    altExp(tse, "test") <- tse[1:1000, ]
    # Checking same name that is already present
    expect_warning(addCluster(tse,
        KmeansParam(centers = 3),
        assay.type = "counts",
        name = "custommetadata",
        full = TRUE
    ))
    # Checking wrong clust.col with already-present name
    expect_warning(addCluster(tse,
        KmeansParam(centers = 3),
        assay.type = "counts",
        clust.col = "customdataname"
    ))
    # Checking wrong altexp
    expect_error(addCluster(tse,
        KmeansParam(centers = 3),
        assay.type = "counts",
        altexp = "error"
    ))

    # Parameters check
    tse <- GlobalPatterns
    altExp(tse, "test") <- tse[1:1000, ]
    tse <- addCluster(tse,
        KmeansParam(centers = 3),
        assay.type = "counts",
        name = "custommetadata",
        full = TRUE,
        clust.col = "customdataname"
    )
    tse <- addCluster(tse,
        KmeansParam(centers = 3),
        assay.type = "counts",
        name = "custommetadata",
        full = TRUE,
        clust.col = "customdataname",
        altexp = "test"
    )
    # Checking custom metadata/dataname in main/altExp
    expect_true("custommetadata" %in% names(metadata(tse)))
    expect_true("customdataname" %in% names(rowData(tse)))
    expect_true("custommetadata" %in% names(metadata(altExp(tse, "test"))))
    expect_true("customdataname" %in% names(rowData(altExp(tse, "test"))))

    # Checking existing custom metadata/dataname in main/altExp
    expect_warning(addCluster(tse,
        KmeansParam(centers = 3),
        assay.type = "counts",
        name = "custommetadata",
        full = TRUE
    ))
    expect_warning(addCluster(tse,
        KmeansParam(centers = 3),
        assay.type = "counts",
        clust.col = "customdataname"
    ))
    expect_warning(addCluster(tse,
        KmeansParam(centers = 3),
        assay.type = "counts",
        name = "custommetadata",
        full = TRUE,
        clust.col = "customdataname",
        altexp = "test"
    ))
    # Checking working by
    tse <- GlobalPatterns
    altExp(tse, "test") <- tse[1:1000, ]
    tse <- addCluster(tse,
        KmeansParam(centers = 3),
        assay.type = "counts",
        by = "col"
    )
    tse <- addCluster(tse,
        KmeansParam(centers = 3),
        assay.type = "counts",
        altexp = "test",
        by = 2
    )
    expect_true("cluster" %in% names(colData(tse)))
    expect_true("cluster" %in% names(colData(altExp(tse, "test"))))

    # Checking wrapper operational
    tse <- GlobalPatterns
    altExp(tse, "test") <- tse[1:2000, ]
    tse <- addCluster(tse,
        HclustParam(),
        assay.type = "counts",
        by = "col"
    )
    tse <- addCluster(tse,
        HclustParam(),
        assay.type = "counts",
        by = "row",
        altexp = "test",
        full = TRUE
    )
    expectedCol <- clusterRows(t(assay(tse, "counts")), HclustParam())
    expectedRow <- clusterRows(assay(altExp(tse, "test"), "counts"),
        HclustParam(),
        full = TRUE
    )
    # Checking same output on cols
    expect_identical(expectedCol, colData(tse)$cluster)
    # Checking same output on rows
    expect_identical(
        expectedRow$cluster,
        rowData(altExp(tse, "test"))$cluster
    )
    # Checking same metdata output on rows
    expect_identical(
        expectedRow$objects,
        metadata(altExp(tse, "test"))$cluster
    )

    # Check that dimred works
    reducedDim(tse, "random") <- assay(tse)[sample(rownames(tse), 100), ] |> t()
    getCluster(tse,
        HclustParam(),
        dimred = "random",
        by = "row"
    ) |> expect_error()
    res <- getCluster(tse,
        HclustParam(),
        dimred = "random",
        by = "col",
        name = "pca"
    )
    expected <- clusterRows(reducedDim(tse, "random"), HclustParam())
    expect_equal(res, expected)
})
