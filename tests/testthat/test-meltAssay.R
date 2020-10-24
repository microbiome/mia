context("meltAssay")
test_that("meltAssay", {

    data(GlobalPatterns, package="MicrobiomeExperiment")
    se <- GlobalPatterns
    molten_assay <- meltAssay(se,
                              add_col_data = c("X.SampleID", "Primer"),
                              abund_values = "counts")
    expect_equal(class(molten_assay)[1], "tbl_df")
    expect_equal(class(molten_assay)[2], "tbl")
    expect_equal(class(molten_assay)[3], "data.frame")
    expect_equal(colnames(molten_assay)[1], "uTaxaID")
    expect_equal(colnames(molten_assay)[2], "uSamId")
    expect_equal(colnames(molten_assay)[3], "Abundance")
    expect_equal(colnames(molten_assay)[11], "X.SampleID")
    expect_equal(is.numeric(molten_assay$Abundance), TRUE)

    only_assay <- .melt_assay(se, abund_values = "counts")
    expect_equal(colnames(only_assay)[1], "uTaxaID")
    expect_equal(colnames(only_assay)[2], "uSamId")
    expect_equal(colnames(only_assay)[3], "Abundance")
    expect_equal(is.numeric(only_assay$Abundance), TRUE)

    assay_taxa <- .add_taxonomic_data_to_molten_assay(se, only_assay)
    expect_equal(colnames(assay_taxa)[1], "uTaxaID")
    expect_equal(colnames(assay_taxa)[2], "uSamId")
    expect_equal(colnames(assay_taxa)[3], "Abundance")
    expect_equal(colnames(assay_taxa)[4], "Kingdom")
    expect_equal(is.numeric(assay_taxa$Abundance), TRUE)

    assay_taxa_coldata <- .add_col_data_to_molten_assay(assay_taxa,
                                                              se,
                                                              add_col_data=c("X.SampleID", "Primer"))
    expect_equal(colnames(assay_taxa_coldata)[1], "uTaxaID")
    expect_equal(colnames(assay_taxa_coldata)[2], "uSamId")
    expect_equal(colnames(assay_taxa_coldata)[3], "Abundance")
    expect_equal(colnames(assay_taxa_coldata)[4], "Kingdom")
    expect_equal(colnames(assay_taxa_coldata)[11], "X.SampleID")
    expect_equal(is.numeric(assay_taxa_coldata$Abundance), TRUE)

})
