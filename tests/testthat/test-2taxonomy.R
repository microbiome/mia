context("taxonomy")
test_that("taxonomy", {
    gr <- GRanges("chr1",rep("1-6",9))
    df <- DataFrame(n = c(1:9))
    mcols(gr) <- df
    grl <- splitAsList(gr,1:9)
    mat <- matrix(1:90, nrow = 9)
    xtse <- TreeSummarizedExperiment(assays = list(mat = mat),
                                     rowRanges = unname(grl))
    tax_data <- DataFrame(Phylum = c(rep("a",3),rep("b",3),rep("c",3)),
                          score = 1:9,
                          Family = c("c",NA,"d","e","f","g","h",NA,"h"),
                          n = 7:15)
    rowData(xtse) <- tax_data
    # .get_tax_cols_logical
    expect_error(MicrobiomeExperiment:::.get_tax_cols_logical(),
                 'argument "x" is missing')
    actual <- MicrobiomeExperiment:::.get_tax_cols_logical(colnames(rowData(xtse)))
    expect_type(actual,"logical")
    expect_equal(actual,c(TRUE,FALSE,TRUE,FALSE))
    # .get_tax_cols
    expect_error(MicrobiomeExperiment:::.get_tax_cols(),
                 'argument "x" is missing')
    actual <- MicrobiomeExperiment:::.get_tax_cols(colnames(rowData(xtse)))
    expect_type(actual,"integer")
    expect_equal(actual,c(1,3))
    # .get_tax_cols_from_se
    expect_error(MicrobiomeExperiment:::.get_tax_cols_from_se(),
                 'argument "x" is missing')
    actual <- MicrobiomeExperiment:::.get_tax_cols_from_se(xtse)
    expect_type(actual,"integer")
    expect_equal(actual,c(1,3))
    # .get_tax_groups
    expect_error(MicrobiomeExperiment:::.get_tax_groups(),
                 'argument "x" is missing')
    expect_error( MicrobiomeExperiment:::.get_tax_groups(xtse),
                  'argument "col" is missing')
    actual <- MicrobiomeExperiment:::.get_tax_groups(xtse,1)
    expect_true(is.factor(actual))
    expect_length(actual,nrow(xtse))
    expect_equal(as.character(actual),tax_data$Phylum)
    actual <- MicrobiomeExperiment:::.get_tax_groups(xtse,2)
    expect_true(is.factor(actual))
    expect_length(actual,nrow(xtse))
    expect_equal(levels(actual),c("a_c","a_NA","a_d","b_e","b_f","b_g","c_h","c_NA"))
    # .get_taxonomic_label
    expect_error(MicrobiomeExperiment:::.get_taxonomic_label(),
                 'argument "x" is missing')
    actual <- MicrobiomeExperiment:::.get_taxonomic_label(xtse, empty.fields = c(NA))
    expect_type(actual,"character")
    expect_length(actual,nrow(xtse))
    expect_equal(actual[1:3],c("Family::c","Phylum::a","Family::d"))
    actual <- MicrobiomeExperiment:::.get_taxonomic_label(xtse[-c(2,8),], empty.fields = c(NA))
    expect_type(actual,"character")
    expect_length(actual,nrow(xtse)-2)
    expect_equal(actual[1:3],c("c","d","e"))
    #
    expect_equal(taxonomyRankEmpty(xtse),
                 rep(FALSE,9))
    expect_true(checkTaxonomy(xtse))
    expect_equal(getTaxonomyLabels(xtse),
                 c("Family::c","Phylum::a","Family::d","Family::e","Family::f",
                   "Family::g","Family::h","Phylum::c","Family::h_1"))
})
