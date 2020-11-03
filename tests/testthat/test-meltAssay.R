context("meltAssay")
test_that("meltAssay", {
    # .norm_add_row_data .norm_add_col_data
    expect_error(mia:::.norm_add_row_data(),
                 'argument "add_row_data" is missing')
    expect_error(mia:::.norm_add_row_data(TRUE),
                 'argument "x" is missing')
    expect_error(mia:::.norm_add_col_data(),
                 'argument "add_col_data" is missing')
    expect_error(mia:::.norm_add_col_data(TRUE),
                 'argument "x" is missing')
    data(GlobalPatterns, package = "MicrobiomeExperiment")
    x <- GlobalPatterns
    actual <- mia:::.norm_add_row_data(TRUE, x)
    expect_equal(actual, colnames(rowData(x)))
    actual <- mia:::.norm_add_col_data(TRUE, x)
    expect_equal(actual, colnames(colData(x)))
    expect_error(mia:::.norm_add_row_data("test", x),
                 "Please provide valid column names")
    expect_error(mia:::.norm_add_col_data("test", x),
                 "Please provide valid column names")
    x2 <- x
    rowData(x2)$FeatureID <- rownames(x2)
    colData(x2)$SampleID <- colnames(x2)
    expect_warning(mia:::.norm_add_row_data(TRUE, x2),
                   "'x' contains a column")
    expect_warning(mia:::.norm_add_col_data(TRUE, x2),
                   "'x' contains a column")
    expect_error(mia:::.norm_add_row_data(NA, x),
                 "'add_row_data' contains NA")
    expect_error(mia:::.norm_add_col_data(NA, x),
                 "'add_col_data' contains NA")
    #

    se <- GlobalPatterns
    molten_assay <- meltAssay(se,
                              add_row_data = TRUE,
                              add_col_data = c("X.SampleID", "Primer"),
                              abund_values = "counts")
    expect_s3_class(molten_assay, c("tbl_df","tbl","data.frame"))
    expect_equal(colnames(molten_assay)[c(1:4,11)], c("FeatureID","SampleID","counts","Kingdom","X.SampleID"))
    expect_equal(is.numeric(molten_assay$counts), TRUE)

    only_assay <- meltAssay(se, abund_values = "counts")
    expect_equal(colnames(only_assay)[1:3], c("FeatureID","SampleID","counts"))
    expect_equal(is.numeric(only_assay$counts), TRUE)

    assay_taxa <- mia:::.add_row_data_to_molten_assay(only_assay,
                                                se,
                                                add_row_data = taxonomyRanks(se))

    expect_equal(colnames(assay_taxa)[1:4], c("FeatureID","SampleID","counts","Kingdom"))
    expect_equal(is.numeric(assay_taxa$counts), TRUE)

    assay_taxa_coldata <- mia:::.add_col_data_to_molten_assay(assay_taxa,
                                                        se,
                                                        add_col_data=c("X.SampleID", "Primer"))

    expect_equal(colnames(molten_assay)[c(1:4,11)], c("FeatureID","SampleID","counts","Kingdom","X.SampleID"))
    expect_equal(is.numeric(assay_taxa_coldata$counts), TRUE)
    #
    actual <- meltAssay(x, TRUE, TRUE)
    expect_warning(actual2 <- meltAssay(x2, TRUE, TRUE))
    expect_false("FeatureID_row" %in% colnames(actual))
    expect_true("FeatureID_row" %in% colnames(actual2))
    expect_false("SampleID_col" %in% colnames(actual))
    expect_true("SampleID_col" %in% colnames(actual2))
    x3 <- x2
    rownames(x3) <- NULL
    colnames(x3) <- NULL
    actual3 <- meltAssay(x3, TRUE, TRUE)
    expect_false("FeatureID_row" %in% colnames(actual))
    expect_false("SampleID_col" %in% colnames(actual))
})
