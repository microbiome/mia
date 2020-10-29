context("meltAssay")
test_that("meltAssay", {

    data(GlobalPatterns, package="MicrobiomeExperiment")
    se <- GlobalPatterns
    molten_assay <- meltAssay(se,
                              add_row_data = TRUE,
                              add_col_data = c("X.SampleID", "Primer"),
                              abund_values = "counts")
    expect_s3_class(molten_assay, c("tbl_df","tbl","data.frame"))
    expect_equal(colnames(molten_assay)[c(1:4,11)], c("FeatureID","SampleID","Abundance","Kingdom","X.SampleID"))
    expect_equal(is.numeric(molten_assay$Abundance), TRUE)

    only_assay <- .melt_assay(se, abund_values = "counts")
    expect_equal(colnames(only_assay)[1:3], c("FeatureID","SampleID","Abundance"))
    expect_equal(is.numeric(only_assay$Abundance), TRUE)

    assay_taxa <- .add_row_data_to_molten_assay(only_assay,
                                                se,
                                                add_row_data = taxonomyRanks(se))

    expect_equal(colnames(assay_taxa)[1:4], c("FeatureID","SampleID","Abundance","Kingdom"))
    expect_equal(is.numeric(assay_taxa$Abundance), TRUE)

    assay_taxa_coldata <- .add_col_data_to_molten_assay(assay_taxa,
                                                        se,
                                                        add_col_data=c("X.SampleID", "Primer"))

    expect_equal(colnames(molten_assay)[c(1:4,11)], c("FeatureID","SampleID","Abundance","Kingdom","X.SampleID"))
    expect_equal(is.numeric(assay_taxa_coldata$Abundance), TRUE)

})
