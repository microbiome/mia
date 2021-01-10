context("import qiime2 results as TreeSummarizedExperiment")


test_that("make TSE worked properly while no sample data or taxa data", {
    featureTableFile <- system.file("extdata", "table.qza", package = "mia")
    expect_silent(makeTreeSummarizedExperimentFromqiime2(featureTableFile))
})

test_that("`.parse_q2taxonomy` work with any combination of taxonomic ranks", {
    test_taxa <- matrix(
        c("a", "k__Bacteria; c__Bacteroidia; s__", 0.88,
          "b", "k__Bacteria; c__Clostridia; s__", 0.9),
        ncol = 3,
        byrow = TRUE,
        dimnames = list(c("a", "b"), c("Feature.ID", "Taxon", "Confidence"))
    )
    expect_silent(.parse_q2taxonomy(test_taxa))

    # at certain rank (e.g. species): some taxa can not be determined which
    # species it assigned to (NA)
    test_taxa <- matrix(
        c("a", "k__Bacteria; c__Bacteroidia; s__", 0.88,
          "b", "k__Bacteria; c__Clostridia;", 0.9),
        ncol = 3,
        byrow = TRUE,
        dimnames = list(c("a", "b"), c("Feature.ID", "Taxon", "Confidence"))
    )
    expect_true(is.na(.parse_q2taxonomy(test_taxa)[2, 3]))
})

test_that("`.read_q2sample_meta` remove  the row contained `#q2:types`", {
    sampleMetaFile <- system.file("extdata", "sample-metadata.tsv", package = "mia")
    expect_false(any(as(.read_q2sample_meta(sampleMetaFile), "matrix") == "#q2:types"))
})

test_that('get file extension', {
    expect_identical(.get_ext("a.b.c"), "c")
    expect_identical(.get_ext("abc/a.b.c"), "c")
})

test_that('read qza file', {
    sample_file <- system.file("extdata", "sample-metadata.tsv", package = "mia")
    expect_error(.read_qza("abc"), "does not exist")
    expect_error(.read_qza(sample_file), "must be in `qza` format")
})
