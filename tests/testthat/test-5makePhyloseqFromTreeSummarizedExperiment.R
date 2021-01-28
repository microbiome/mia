context("makePhyloseqFromTreeSummarizedExperiment")

test_that("makePhyloseqFromTreeSummarizedExperiment", {

    skip_if_not_installed("phyloseq")

    # TSE object
    data(GlobalPatterns)
    tse <- GlobalPatterns

    phy <- makePhyloseqFromTreeSummarizedExperiment(GlobalPatterns)

    # Test that assay is in otu_table
    testthat::expect_equal(as.data.frame(phyloseq::otu_table(phy)@.Data), as.data.frame(assays(tse)$counts))

    # Test that rowData is in tax_table
    testthat::expect_equal(as.data.frame(phyloseq::tax_table(phy)@.Data), as.data.frame(rowData(tse)))

    # Test that colData is in sample_table
    testthat::expect_equal(phyloseq::sample_data(phy),
                           phyloseq::sample_data(data.frame(colData(tse))))

    # Test that rowTree is in phy_tree
    testthat::expect_equal(phyloseq::phy_tree(phy), rowTree(tse))

    # Test that referenceSeq is in refseq. There is no referenceseq in TSE, so
    # referenceseq slot is not made in the phyloseq object
    #testthat::expect_equal(phyloseq::refseq(phy), referenceSeq(tse))

    # TSE object
    data("esophagus")
    tse <- esophagus

    phy <- makePhyloseqFromTreeSummarizedExperiment(esophagus)

    # Test that assay is in otu_table
    testthat::expect_equal(as.data.frame(phyloseq::otu_table(phy)@.Data), as.data.frame(assays(tse)$counts))


    # Test that rowTree is in phy_tree
    testthat::expect_equal(phyloseq::phy_tree(phy), rowTree(tse))



})


