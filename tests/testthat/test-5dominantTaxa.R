context("dominantTaxa")

test_that("dominantTaxa", {

    # TSE object
    data(GlobalPatterns)
    tse <- GlobalPatterns

    # Phyloseq object
    data("GlobalPatterns", package = "phyloseq")
    phy <- GlobalPatterns

    # Counts table should not be changed
    expect_equal(assays(mia::dominantTaxa(tse))$counts, assays(tse)$counts)

    # Should get same result
    expect_equal(colData(mia::dominantTaxa(tse, name="dominant"))$dominant, microbiome::dominant(phy))

    expect_equal(colData(mia::dominantTaxa(tse, rank = "Species"))$dominant_taxa, microbiome::dominant(phy, level = "Species"))

    expect_equal(colData(mia::dominantTaxa(tse, rank = "Genus"))$dominant_taxa, microbiome::dominant(phy, level = "Genus"))

    expect_equal(colData(mia::dominantTaxa(tse, rank = "Family"))$dominant_taxa, microbiome::dominant(phy, level = "Family"))

    expect_equal(colData(mia::dominantTaxa(tse, rank = "Order"))$dominant_taxa, microbiome::dominant(phy, level = "Order"))

    expect_equal(colData(mia::dominantTaxa(tse, rank = "Class"))$dominant_taxa, microbiome::dominant(phy, level = "Class"))

    expect_equal(colData(mia::dominantTaxa(tse, rank = "Phylum"))$dominant_taxa, microbiome::dominant(phy, level = "Phylum"))

    expect_equal(colData(mia::dominantTaxa(tse, rank = "Kingdom"))$dominant_taxa, microbiome::dominant(phy, level = "Kingdom"))

})


test_that("getDominantTaxa", {

    # TSE object
    tse <- x <- microbiomeDataSets::dietswap()

    # Phyloseq object
    data("dietswap", package = "microbiome")
    phy <- dietswap

    expect_equal(mia::getDominantTaxa(tse), microbiomeutilities::dominant_taxa(phy)$dominant_overview)

    expect_equal(mia::getDominantTaxa(tse, rank = "Family"), microbiomeutilities::dominant_taxa(phy, level = "Family")$dominant_overview)

    expect_equal(mia::getDominantTaxa(tse, rank = "Family", group = "nationality"),
                 microbiomeutilities::dominant_taxa(phy, level = "Family", group = "nationality")$dominant_overview)






})
