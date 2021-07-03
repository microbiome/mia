
test_that("Importing biom files yield SummarizedExperiment objects", {
    skip_if_not_installed("biomformat")
    rich_dense_file  = system.file("extdata", "rich_dense_otu_table.biom",
                                   package = "biomformat")
    me <- loadFromBiom(rich_dense_file)
    expect_s4_class(me, "SummarizedExperiment")
    # load from object
    x1 <- biomformat::read_biom(rich_dense_file)
    me2 <- makeSummarizedExperimentFromBiom(x1)
    expect_s4_class(me2, "SummarizedExperiment")
    expect_equal(dim(me), dim(me2))
    expect_equal(rowData(me), rowData(me2))
})

test_that("Importing phyloseq objects yield TreeSummarizedExperiment objects", {
    skip_if_not_installed("phyloseq")
    data(GlobalPatterns, package="phyloseq")
    me <- makeTreeSummarizedExperimentFromphyloseq(GlobalPatterns)
    expect_s4_class(me, "TreeSummarizedExperiment")
    expect_equal(dim(me),c(19216,26))
    data(enterotype, package="phyloseq")
    me <- makeTreeSummarizedExperimentFromphyloseq(enterotype)
    expect_s4_class(me, "TreeSummarizedExperiment")
    expect_equal(dim(me),c(553,280))
    data(esophagus, package="phyloseq")
    me <- makeTreeSummarizedExperimentFromphyloseq(esophagus)
    expect_s4_class(me, "TreeSummarizedExperiment")
    expect_equal(dim(me),c(58,3))
})

test_that("Importing dada2 objects yield TreeSummarizedExperiment objects", {
    skip_if_not_installed("dada2")
    fnF <- system.file("extdata", "sam1F.fastq.gz", package="dada2")
    fnR = system.file("extdata", "sam1R.fastq.gz", package="dada2")
    dadaF <- dada2::dada(fnF, selfConsist=TRUE)
    dadaR <- dada2::dada(fnR, selfConsist=TRUE)

    me <- makeTreeSummarizedExperimentFromDADA2(dadaF, fnF, dadaR, fnR)
    expect_s4_class(me, "TreeSummarizedExperiment")
})

test_that("Importing Mothur files yield SummarizedExperiment objects", {
    
    counts <- system.file("extdata", "mothur_example.shared", package = "mia")
    taxa <- system.file("extdata", "mothur_example.cons.taxonomy", package = "mia")
    taxa2 <- system.file("extdata", "mothur_example.taxonomy", package = "mia")
    meta <- system.file("extdata", "mothur_example.design", package = "mia")
    se <- loadFromMothur(counts)
    se2 <- loadFromMothur(counts)
    expect_s4_class(se, "SummarizedExperiment")
    expect_s4_class(se2, "SummarizedExperiment")
    se <- loadFromMothur(counts, taxa)
    se2 <- loadFromMothur(counts, taxa2)
    expect_s4_class(se, "SummarizedExperiment")
    expect_s4_class(se2, "SummarizedExperiment")
    expect_error(loadFromMothur(counts, meta))
    expect_error(loadFromMothur(counts, meta))
    se <- loadFromMothur(counts, designFile = meta)
    se2 <- loadFromMothur(counts, designFile = meta)
    expect_s4_class(se, "SummarizedExperiment")
    expect_s4_class(se2, "SummarizedExperiment")
    se <- loadFromMothur(counts, taxa, meta)
    se2 <- loadFromMothur(counts, taxa2, meta)
    expect_s4_class(se, "SummarizedExperiment")
    expect_s4_class(se2, "SummarizedExperiment")
    
    # Checks dimensions, rownames, and colnames of assay
    expect_equal(nrow(assays(se)$counts), 100)
    expect_equal(rownames(assays(se)$counts)[1:10],
                          c("Otu001", "Otu002",  "Otu003",  "Otu004",  "Otu005",
                            "Otu006", "Otu007", "Otu008", "Otu009", "Otu010"))
    expect_equal(ncol(assays(se)$counts), 100)
    expect_equal(colnames(assays(se)$counts)[1:10],
                          c("Sample1", "Sample2", "Sample3", "Sample4", "Sample5", 
                            "Sample6", "Sample7", "Sample8", "Sample9", "Sample10"))
    expect_equal(nrow(assays(se2)$counts), 100)
    expect_equal(rownames(assays(se2)$counts)[1:10],
                 c("Otu001", "Otu002",  "Otu003",  "Otu004",  "Otu005",
                   "Otu006", "Otu007", "Otu008", "Otu009", "Otu010"))
    expect_equal(ncol(assays(se2)$counts), 100)
    expect_equal(colnames(assays(se)$counts)[1:10],
                          c("Sample1", "Sample2", "Sample3", "Sample4", "Sample5", 
                            "Sample6", "Sample7", "Sample8", "Sample9", "Sample10"))
    
    # Checks that rowData has right dimensions, rownames, and colnames
    expect_equal(nrow(rowData(se)), 100)
    expect_equal(rownames(rowData(se))[1:10],
                          c("Otu001", "Otu002",  "Otu003",  "Otu004",  "Otu005",
                            "Otu006", "Otu007", "Otu008", "Otu009", "Otu010"))
    expect_equal(colnames(rowData(se)),
                 c("OTU", "Size", "Kingdom", "Phylum", "Order", "Class", "Family", "Genus"))
    expect_equal(nrow(rowData(se2)), 100)
    expect_equal(rownames(rowData(se2))[1:10],
                          c("Otu001", "Otu002",  "Otu003",  "Otu004",  "Otu005",
                            "Otu006", "Otu007", "Otu008", "Otu009", "Otu010"))
    expect_equal(colnames(rowData(se2)),
                 c("OTU", "Kingdom", "Phylum", "Order", "Class", "Family", "Genus"))
    
    expect_equal(nrow(colData(se)), 100)
    expect_equal(rownames(colData(se))[1:10],
                          c("Sample1", "Sample2", "Sample3", "Sample4", "Sample5", 
                            "Sample6", "Sample7", "Sample8", "Sample9", "Sample10"))
    
    # Checks colData's dimensions and names of columns and rows
    expect_equal(colnames(colData(se)),
                          c("group", "sex", "age", "drug", "label", "numOtus", "Group"))
    expect_equal(nrow(colData(se2)), 100)
    expect_equal(rownames(colData(se2))[1:10],
                          c("Sample1", "Sample2", "Sample3", "Sample4", "Sample5", 
                            "Sample6", "Sample7", "Sample8", "Sample9", "Sample10"))
    expect_equal(colnames(colData(se2)),
                          c("group", "sex", "age", "drug", "label", "numOtus", "Group"))
})

featureTableFile <- system.file("extdata", "table.qza", package = "mia")
taxonomyTableFile <- system.file("extdata", "taxonomy.qza", package = "mia")
refSeqFile <- system.file("extdata", "refseq.qza", package = "mia")

test_that("make TSE worked properly while no sample or taxa data", {
    ## no sample data or taxa data
    expect_silent(tse <- loadFromQIIME2(featureTableFile))
    expect_s4_class(tse, "TreeSummarizedExperiment")
    expect_equal(dim(tse), c(770,34))
})

test_that("reference sequences of TSE", {
    # 1. fasta file of refseq
    tse <- loadFromQIIME2(
        featureTableFile,
        refSeqFile = refSeqFile
    )
    tse2 <-  loadFromQIIME2(
        featureTableFile,
        refSeqFile = refSeqFile,
        featureNamesAsRefseq = FALSE
    )
    expect_identical(tse@referenceSeq, .read_qza(refSeqFile))
    expect_identical(tse2@referenceSeq, .read_qza(refSeqFile))

    # 2. row.names of feature table as refseq
    # 2.1 element of row.names of feature table is not DNA sequence
    tse <- loadFromQIIME2(
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
        # create qza file using QIIME2
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
    tse <- loadFromQIIME2(
        featureTableFile2,
        featureNamesAsRefseq = TRUE
    )
    feature_tab <- .read_qza(featureTableFile2)
    names_seq <- Biostrings::DNAStringSet(row.names(feature_tab))
    names(names_seq) <- paste0("seq_", seq_along(names_seq))
    expect_identical(tse@referenceSeq, names_seq)

    # refSeqFile is not NULL, featureNamesAsRefseq is TRUE,
    # set the sequences from refSeqFile as reference sequences
    tse <- loadFromQIIME2(
        featureTableFile2,
        featureNamesAsRefseq = TRUE,
        refSeqFile = refSeqFile
    )
    expect_identical(tse@referenceSeq, .read_qza(refSeqFile))

    # 3. refSeqFile = NULL, featureNamesAsRefseq = FALSE
    tse <- loadFromQIIME2(
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
    expect_silent(mia:::.parse_q2taxonomy(test_taxa))

    # at certain rank (e.g. species): some taxa can not be determined which
    # species it assigned to (NA)
    test_taxa <- matrix(
        c("a", "k__Bacteria; c__Bacteroidia; s__", 0.88,
          "b", "k__Bacteria; c__Clostridia;", 0.9),
        ncol = 3,
        byrow = TRUE,
        dimnames = list(c("a", "b"), c("Feature.ID", "Taxon", "Confidence"))
    )
    expect_true(is.na(mia:::.parse_q2taxonomy(test_taxa)[2,"Species"]))

    # if the expexted order is not present it will return a correct result
    test_taxa <- matrix(
        c("a", "k__Bacteria; s__test; c__Bacteroidia", 0.88,
          "b", "k__Bacteria; c__Clostridia;", 0.9),
        ncol = 3,
        byrow = TRUE,
        dimnames = list(c("a", "b"), c("Feature.ID", "Taxon", "Confidence"))
    )
    expect_equal(mia:::.parse_q2taxonomy(test_taxa)[,"Species"],c("s__test",NA))
    expect_equal(mia:::.parse_q2taxonomy(test_taxa, removeTaxaPrefixes = TRUE)[,"Species"],
                 c("test",NA))
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
    tse <- loadFromQIIME2(
        featureTableFile,
        taxonomyTableFile = taxonomyTableFile
    )
    expect_true(is.numeric(S4Vectors::mcols(tse)$Confidence))
})

test_that("makePhyloseqFromTreeSummarizedExperiment", {

    skip_if_not_installed("phyloseq")

    # TSE object
    data(GlobalPatterns)
    tse <- GlobalPatterns

    phy <- makePhyloseqFromTreeSummarizedExperiment(GlobalPatterns)

    # Test that assay is in otu_table
    expect_equal(as.data.frame(phyloseq::otu_table(phy)@.Data), as.data.frame(assays(tse)$counts))

    # Test that rowData is in tax_table
    expect_equal(as.data.frame(phyloseq::tax_table(phy)@.Data), as.data.frame(rowData(tse)))

    # Test that colData is in sample_table
    expect_equal(phyloseq::sample_data(phy),
                 phyloseq::sample_data(data.frame(colData(tse))))

    # Test that rowTree is in phy_tree
    expect_equal(phyloseq::phy_tree(phy), rowTree(tse))

    # Test that referenceSeq is in refseq. Expect error, because there should not be
    # reference sequences.
    expect_error(phyloseq::refseq(phy))
    
    # Test with agglomeration that an error occurs when rowTree is incorrect
    test1 <- agglomerateByRank(tse, rank = "Phylum")
    test2 <- expect_warning(agglomerateByRank(tse, rank = "Phylum", agglomerateTree = TRUE))
    test2_phy <- makePhyloseqFromTreeSummarizedExperiment(test2)
    
    expect_error(makePhyloseqFromTreeSummarizedExperiment(test1))
    expect_equal(phyloseq::phy_tree(test2_phy), rowTree(test2))

    # TSE object
    data(esophagus)
    tse <- esophagus

    phy <- makePhyloseqFromTreeSummarizedExperiment(esophagus)

    # Test that assay is in otu_table
    expect_equal(as.data.frame(phyloseq::otu_table(phy)@.Data), as.data.frame(assays(tse)$counts))


    # Test that rowTree is in phy_tree
    expect_equal(phyloseq::phy_tree(phy), rowTree(tse))

})
