context("import qiime2 results as TreeSummarizedExperiment")

featureTableFile <- system.file("extdata", "table.qza", package = "mia")
taxonomyTableFile <- system.file("extdata", "taxonomy.qza", package = "mia")
refSeqFile <- system.file("extdata", "refseq.qza", package = "mia")

test_that("make TSE worked properly while no sample or taxa data", {
    ## no sample data or taxa data
    expect_silent(tse <- makeTreeSummarizedExperimentFromqiime2(featureTableFile))
    expect_s4_class(tse, "TreeSummarizedExperiment")
    expect_equal(dim(tse), c(770,34))
})

test_that("reference sequences of TSE", {
    # 1. fasta file of refseq
    tse <- makeTreeSummarizedExperimentFromqiime2(
        featureTableFile,
        refSeqFile = refSeqFile
    )
    tse2 <-  makeTreeSummarizedExperimentFromqiime2(
        featureTableFile,
        refSeqFile = refSeqFile,
        featureNamesAsRefseq = FALSE
    )
    expect_identical(tse@referenceSeq, .read_qza(refSeqFile))
    expect_identical(tse2@referenceSeq, .read_qza(refSeqFile))

    # 2. row.names of feature table as refseq
    # 2.1 element of row.names of feature table is not DNA sequence
    tse <- makeTreeSummarizedExperimentFromqiime2(
        featureTableFile,
        featureNamesAsRefseq = TRUE
    )
    expect_null(tse@referenceSeq)

    # 2.2 element of row.names of feature table is a DNA sequence
    #
    # codes used for create sample data (donot run)
    if (FALSE) {
        .require_package("biomformat")
        feature_tab <- .read_qza(featureTableFile)
        n_feature <- nrow(feature_tab)
        random_seq <- sapply(
            rep(20, n_feature),
            function(size) {
                paste(
                    sample(Biostrings::DNA_BASES, size, replace = TRUE),
                    collapse=""
                )
            }
        )
        row.names(feature_tab) <- random_seq
        obj <- biomformat::make_biom(feature_tab)
        biomformat::write_biom(obj, "data/table-rownamesInSeq.biom")
        # create qza file using qiime2
        # qiime tools import \
        #     --input-path data/table-rownamesInSeq.biom
        #     --type 'FeatureTable[Frequency]'
        #     --input-format BIOMV100Format
        #     --output-path inst/extdata/table-rownamesInSeq.qza
        unlink("data/table-rownamesInSeq.biom")
    }
    featureTableFile2 <- system.file(
        "extdata",
        "table-rownamesInSeq.qza",
        package = "mia"
    )

    # featureNamesAsRefseq is TRUE, refSeqFile is NULL, set row.names of
    # feature table as reference sequences
    tse <- makeTreeSummarizedExperimentFromqiime2(
        featureTableFile2,
        featureNamesAsRefseq = TRUE
    )
    feature_tab <- .read_qza(featureTableFile2)
    names_seq <- Biostrings::DNAStringSet(row.names(feature_tab))
    names(names_seq) <- paste0("seq_", seq_along(names_seq))
    expect_identical(tse@referenceSeq, names_seq)

    # refSeqFile is not NULL, featureNamesAsRefseq is TRUE,
    # set the sequences from refSeqFile as reference sequences
    tse <- makeTreeSummarizedExperimentFromqiime2(
        featureTableFile2,
        featureNamesAsRefseq = TRUE,
        refSeqFile = refSeqFile
    )
    expect_identical(tse@referenceSeq, .read_qza(refSeqFile))

    # 3. refSeqFile = NULL, featureNamesAsRefseq = FALSE
    tse <- makeTreeSummarizedExperimentFromqiime2(
        featureTableFile,
        refSeqFile = NULL,
        featureNamesAsRefseq = FALSE
    )
    expect_null(tse@referenceSeq)
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

test_that("Confidence of taxa is numberic", {
    tse <- makeTreeSummarizedExperimentFromqiime2(
        featureTableFile,
        taxonomyTableFile = taxonomyTableFile
    )
    expect_true(is.numeric(S4Vectors::mcols(tse)$Confidence))
})
