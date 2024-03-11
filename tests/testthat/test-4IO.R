
test_that("Importing biom files yield SummarizedExperiment objects", {
    skip_if_not(require("biomformat", quietly = TRUE))
    rich_dense_file  = system.file("extdata", "rich_dense_otu_table.biom",
                                   package = "biomformat")
    me <- loadFromBiom(rich_dense_file)
    expect_s4_class(me, "SummarizedExperiment")
    # load from object
    x1 <- biomformat::read_biom(rich_dense_file)
    me2 <- makeTreeSEFromBiom(x1)
    expect_s4_class(me2, "SummarizedExperiment")
    expect_equal(dim(me), dim(me2))
    expect_equal(rowData(me), rowData(me2))
    
    biom_object <- biomformat::read_biom(
        system.file("extdata/testdata/Aggregated_humanization2.biom",
                    package="mia")
    )
    tse <- makeTreeSEFromBiom(biom_object,
                              removeTaxaPrefixes = FALSE,
                              rankFromPrefix = FALSE,
                              remove.artifacts = TRUE,
                              pattern = "\"")
    # Testing no prefixes removed
    expect_true(rowData(tse) %>%
                    apply(2,grepl,pattern="sk__|([dkpcofgs]+)__") %>%
                    all())
    # Testing no taxonomy ranks parsed
    expect_false(
        sapply(tolower(colnames(rowData(tse))),
               function(x) x %in% TAXONOMY_RANKS) %>% all())
    # Testing the remove.artifacts, since the original artifact in the biom file 
    # is '\"'
    expect_false(apply(rowData(tse), 2, grepl, pattern="^\"") %>% all())
    
    # Testing prefixes removed
    tse <- makeTreeSEFromBiom(biom_object,
                              removeTaxaPrefixes=TRUE,
                              rankFromPrefix=TRUE,
                              remove.artifacts = TRUE,
                              pattern = "\"")
    expect_false(rowData(tse) %>%
                     apply(2,grepl,pattern="sk__|([dkpcofgs]+)__") %>%
                     all())
    
    # Testing parsing taxonomy ranks from prefixes
    tse <- makeTreeSEFromBiom(biom_object,
                              removeTaxaPrefixes=FALSE,
                              rankFromPrefix=TRUE,
                              remove.artifacts = TRUE,
                              pattern = "\"")
    expect_true(
        sapply(tolower(colnames(rowData(tse))),
               function(x) x %in% TAXONOMY_RANKS) %>% all())
    
    # Testing the remove.artifacts, the original artifact in the biom file 
    # is '\"', as a test we rather try remove a non existing pattern.
    tse <- makeTreeSEFromBiom(biom_object,
                              removeTaxaPrefixes=FALSE,
                              rankFromPrefix=FALSE,
                              remove.artifacts = TRUE,
                              pattern = "\\*|\\?")
    # with wrong pattern artifact not cleaned
    expect_true(apply(rowData(tse), 2, grepl, pattern="\"") %>% any())
    # Testing the remove.artifacts, with the value 'auto' to automatically 
    # detect the artifact and remove it (in our case the artifact is '\"').
    tse <- makeTreeSEFromBiom(biom_object,
                              removeTaxaPrefixes=FALSE,
                              rankFromPrefix=FALSE,
                              remove.artifacts = TRUE)
    # Checking if 'auto' has detected and cleaned the artifact
    expect_false(apply(rowData(tse), 2, grepl, pattern="\"") %>% any())
    # Testing the remove.artifacts, with the value NULL to not detect or clean 
    # anything.
    tse <- makeTreeSEFromBiom(biom_object,
                              removeTaxaPrefixes=FALSE,
                              rankFromPrefix=FALSE,
                              remove.artifacts = FALSE)
    # Checking if the '\"' artifact still exists.
    expect_true(apply(rowData(tse), 2, grepl, pattern="\"") %>% any())
    
    # General final test
    tse <- makeTreeSEFromBiom(biom_object,
                              removeTaxaPrefixes=TRUE,
                              rankFromPrefix=TRUE,
                              remove.artifacts = TRUE)
    # check if '\"' cleaned
    expect_false(apply(rowData(tse), 2, grepl, pattern="\"") %>% any())
    # check if taxa prefixes removed
    expect_false(rowData(tse) %>%
                     apply(2,grepl,pattern="sk__|([dkpcofgs]+)__") %>%
                     all())
    # Check if rank names were parsed correctly
    expect_true(
        sapply(tolower(colnames(rowData(tse))),
               function(x) x %in% TAXONOMY_RANKS) %>% all())
    
    # General final test with another biom file
    biom_object <- biomformat::read_biom(
        system.file("extdata", "rich_dense_otu_table.biom",
                    package = "biomformat")
    )
    tse <- makeTreeSEFromBiom(biom_object,
                              removeTaxaPrefixes=TRUE,
                              rankFromPrefix=TRUE,
                              remove.artifacts = TRUE)
    # check if taxa prefixes removed
    expect_false(rowData(tse) %>%
                     apply(2,grepl,pattern="sk__|([dkpcofgs]+)__") %>%
                     all())
    # Check if rank names were parsed correctly
    expect_true(
        sapply(tolower(colnames(rowData(tse))),
               function(x) x %in% TAXONOMY_RANKS) %>% all())
})

test_that("Importing phyloseq objects yield TreeSummarizedExperiment objects", {
    skip_if_not_installed("phyloseq")
    data(GlobalPatterns, package="phyloseq")
    me <- makeTreeSEFromPhyloseq(GlobalPatterns)
    expect_s4_class(me, "TreeSummarizedExperiment")
    expect_equal(dim(me),c(19216,26))
    data(enterotype, package="phyloseq")
    me <- makeTreeSEFromPhyloseq(enterotype)
    expect_s4_class(me, "TreeSummarizedExperiment")
    expect_equal(dim(me),c(553,280))
    data(esophagus, package="phyloseq")
    me <- makeTreeSEFromPhyloseq(esophagus)
    expect_s4_class(me, "TreeSummarizedExperiment")
    expect_equal(dim(me),c(58,3))
    esophagus2 <- esophagus
    phyloseq::otu_table(esophagus2) <- t(phyloseq::otu_table(esophagus))
    me2 <- makeTreeSEFromPhyloseq(esophagus2)
    expect_equal(me, me2)
})

test_that("Importing dada2 objects yield TreeSummarizedExperiment objects", {
    skip_if_not_installed("dada2")
    fnF <- system.file("extdata", "sam1F.fastq.gz", package="dada2")
    fnR = system.file("extdata", "sam1R.fastq.gz", package="dada2")
    dadaF <- dada2::dada(fnF, selfConsist=TRUE)
    dadaR <- dada2::dada(fnR, selfConsist=TRUE)

    me <- makeTreeSEFromDADA2(dadaF, fnF, dadaR, fnR)
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
sampleMetaFile <- system.file("extdata", "sample-metadata.tsv", package = "mia")
refSeqFile <- system.file("extdata", "refseq.qza", package = "mia")

test_that("make TSE worked properly while no sample or taxa data", {
    skip_if_not(require("biomformat", quietly = TRUE))
    ## no sample data or taxa data
    expect_silent(tse <- loadFromQIIME2(featureTableFile))
    expect_s4_class(tse, "TreeSummarizedExperiment")
    expect_equal(dim(tse), c(770,34))
})

test_that("reference sequences of TSE", {
    skip_if_not(require("biomformat", quietly = TRUE))
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
    expect_identical(tse@referenceSeq, readQZA(refSeqFile))
    expect_identical(tse2@referenceSeq, readQZA(refSeqFile))

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
        feature_tab <- readQZA(featureTableFile)
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
    feature_tab <- readQZA(featureTableFile2)
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
    expect_identical(tse@referenceSeq, readQZA(refSeqFile))

    # 3. refSeqFile = NULL, featureNamesAsRefseq = FALSE
    tse <- loadFromQIIME2(
        featureTableFile,
        refSeqFile = NULL,
        featureNamesAsRefseq = FALSE
    )
    expect_null(tse@referenceSeq)
})


test_that("`.parse_taxonomy` work with any combination of taxonomic ranks", {
    test_taxa <- matrix(
        c("a", "k__Bacteria; c__Bacteroidia; s__", 0.88,
          "b", "k__Bacteria; c__Clostridia; s__", 0.9),
        ncol = 3,
        byrow = TRUE,
        dimnames = list(c("a", "b"), c("Feature.ID", "Taxon", "Confidence"))
    )
    expect_silent(mia:::.parse_taxonomy(test_taxa))

    # at certain rank (e.g. species): some taxa can not be determined which
    # species it assigned to (NA)
    test_taxa <- matrix(
        c("a", "k__Bacteria; c__Bacteroidia; s__", 0.88,
          "b", "k__Bacteria; c__Clostridia;", 0.9),
        ncol = 3,
        byrow = TRUE,
        dimnames = list(c("a", "b"), c("Feature.ID", "Taxon", "Confidence"))
    )
    expect_true(is.na(mia:::.parse_taxonomy(test_taxa)[2,"Species"]))

    # if the expexted order is not present it will return a correct result
    test_taxa <- matrix(
        c("a", "k__Bacteria; s__test; c__Bacteroidia", 0.88,
          "b", "k__Bacteria; c__Clostridia;", 0.9),
        ncol = 3,
        byrow = TRUE,
        dimnames = list(c("a", "b"), c("Feature.ID", "Taxon", "Confidence"))
    )
    expect_equal(mia:::.parse_taxonomy(test_taxa)[,"Species"],c("s__test",NA))
    expect_equal(mia:::.parse_taxonomy(test_taxa, removeTaxaPrefixes = TRUE)[,"Species"],
                 c("test",NA))
})

test_that("`.read_q2sample_meta` remove  the row contained `#q2:types`", {
    expect_false(any(as(.read_q2sample_meta(sampleMetaFile), "matrix") == "#q2:types"))
})

test_that('get file extension', {
    expect_identical(.get_ext("a.b.c"), "c")
    expect_identical(.get_ext("abc/a.b.c"), "c")
})

test_that('read qza file', {
    expect_error(readQZA("abc"), "does not exist")
    expect_error(readQZA(sampleMetaFile), "must be in `qza` format")
})

test_that("Confidence of taxa is numberic", {
    skip_if_not(require("biomformat", quietly = TRUE))
    tse <- loadFromQIIME2(
        featureTableFile,
        taxonomyTableFile = taxonomyTableFile
    )
    expect_true(is.numeric(S4Vectors::mcols(tse)$Confidence))
})

test_that("dimnames of feature table is identicle with meta data", {
   skip_if_not(require("biomformat", quietly = TRUE))
   feature_tab <- readQZA(featureTableFile)
   
   sample_meta <- .read_q2sample_meta(sampleMetaFile)
   taxa_meta <- readQZA(taxonomyTableFile)
   taxa_meta <- .subset_taxa_in_feature(taxa_meta, feature_tab)
   new_feature_tab <- .set_feature_tab_dimnames(
       feature_tab, 
       sample_meta, 
       taxa_meta
   )
   expect_identical(rownames(new_feature_tab), rownames(taxa_meta))
   expect_identical(colnames(new_feature_tab), rownames(sample_meta))
   
   # sample_meta or feature meta is NULL
   sample_meta2 <- S4Vectors::make_zero_col_DFrame(ncol(feature_tab))
   rownames(sample_meta2) <- colnames(feature_tab)
   taxa_meta2 <- S4Vectors::make_zero_col_DFrame(nrow(feature_tab))
   rownames(taxa_meta2) <- rownames(feature_tab)
   expect_silent(.set_feature_tab_dimnames(feature_tab, sample_meta2, taxa_meta))
   
   # sample meta or feature meta without any information, only contains sample/feature
   # ID in its rownames
   feature_tab3 <- S4Vectors::DataFrame(
       sample1 = 1:3,
       sample2 = 4:6,
       sample3 = 7:9,
       row.names = paste0("feature", 1:3)
   )
   sample_meta3 <- S4Vectors::DataFrame(row.names = paste0("sample", 3:1))
   feature_meta3 <- S4Vectors::DataFrame(row.names = paste0("feature", c(2, 3, 1)))
   new_feature_tab3 <- .set_feature_tab_dimnames(
       feature_tab3, 
       sample_meta3, 
       feature_meta3
    )
   expect_identical(row.names(new_feature_tab3), paste0("feature", c(2, 3, 1)))
   expect_identical(colnames(new_feature_tab3), paste0("sample", 3:1))
})


test_that("makePhyloseqFromTreeSE", {

    skip_if_not_installed("phyloseq")

    # TSE object
    data(GlobalPatterns, package="mia")
    tse <- GlobalPatterns

    phy <- makePhyloseqFromTreeSE(GlobalPatterns)

    # Test that assay is in otu_table
    expect_equal(as.data.frame(phyloseq::otu_table(phy)@.Data), as.data.frame(assays(tse)$counts))

    # Test that rowData is in tax_table
    expect_equal(as.data.frame(phyloseq::tax_table(phy)@.Data), as.data.frame(rowData(tse)))

    # Test that colData is in sample_table
    expect_equal(phyloseq::sample_data(phy),
                 phyloseq::sample_data(data.frame(colData(tse))))

    # Test that rowTree is in phy_tree
    expect_identical(phyloseq::phy_tree(phy), rowTree(tse))

    # Test that referenceSeq is in refseq. Expect error, because there should not be
    # reference sequences.
    expect_error(phyloseq::refseq(phy))
    
    # Test with agglomeration that that pruning is done internally
    test1 <- agglomerateByRank(tse, rank = "Phylum")
    test2 <- agglomerateByRank(tse, rank = "Phylum", agglomerate.tree = TRUE)
    test1_phy <- expect_warning(makePhyloseqFromTreeSE(test1))
    test2_phy <- makePhyloseqFromTreeSE(test2)
    
    expect_equal(length(phyloseq::phy_tree(test1_phy)$node), 
                 length(ape::keep.tip(rowTree(test1), rowLinks(test1)$nodeLab)$node))
    expect_equal(phyloseq::phy_tree(test1_phy)$tip.label, rownames(test2))
    # The tip labels do not match because of renaming
    expect_identical(phyloseq::phy_tree(test2_phy)$edge, rowTree(test2)$edge)
    
    # Check that everything works also with agglomerated data
    for (level in colnames(rowData(tse)) ){
        temp <- agglomerateByRank(tse, rank = level)
        expect_warning(makePhyloseqFromTreeSE(temp))
    }
    
    tse2 <- tse
    # Concerts data frame to factors
    rowData(tse2) <- DataFrame(lapply(rowData(tse2), as.factor))
    phy <- makePhyloseqFromTreeSE(tse)
    phy2 <- makePhyloseqFromTreeSE(tse2)
    expect_equal(phyloseq::tax_table(phy2), phyloseq::tax_table(phy))
    
    # TSE object
    data(esophagus, package="mia")
    tse <- esophagus

    phy <- makePhyloseqFromTreeSE(tse, assay.type="counts")

    # Test that assay is in otu_table
    expect_equal(as.data.frame(phyloseq::otu_table(phy)@.Data), as.data.frame(assays(tse)$counts))

    # Test that rowTree is in phy_tree
    expect_identical(phyloseq::phy_tree(phy), rowTree(tse))
    
    # Test that merging objects lead to correct phyloseq
    tse <- mergeSEs(GlobalPatterns, esophagus, assay.type="counts", missing_values = 0)
    pseq <- makePhyloseqFromTreeSE(tse, assay.type="counts")
    
    tse_compare <- tse[ rownames(GlobalPatterns), ]
    pseq_compare <- makePhyloseqFromTreeSE(tse_compare, assay.type="counts")
    
    expect_equal(phyloseq::otu_table(pseq), phyloseq::otu_table(pseq_compare))
})

test_that("Import HUMAnN file", {
    file_path <- system.file("extdata", "humann_output.tsv", package = "mia")
    tse <- loadFromHumann(file_path)
    #
    expect_true( length(taxonomyRanks(tse)) == 7 )
    expect_true( ncol(rowData(tse)) == 9 )
    expect_true( ncol(tse) == 3 )
    expect_true( nrow(tse) == 12 )
    expect_true( all(!is.na(assay(tse))) )
    #
    expect_error(loadFromHumann("test"))
    expect_error(loadFromHumann(file_path, remove.suffix="test"))
    expect_error(loadFromHumann(file_path, remove.suffix=1))
    expect_error(loadFromHumann(file_path, remove.suffix=c(FALSE, TRUE)))
    #
    tse2 <- loadFromHumann(file_path, remove.suffix=TRUE)
    # There is no suffix, should be equal
    expect_equal(tse, tse2)
    #
    cd <- colData(tse)
    tse2 <- loadFromHumann(file_path, colData = cd)
    # colData should not change the result
    expect_equal(tse, tse2)
})
