context("Unifrac beta diversity")
test_that("Unifrac beta diversity", {
    skip_if_not(require("phyloseq", quietly = TRUE))
    
    data(esophagus, package="mia")
    tse <- esophagus
    tse <- transformAssay(tse, assay.type="counts", method="relabundance")
    
    expect_error(
        calculateUnifrac(tse, assay.type = "test", tree_name = "phylo",
                         weighted = FALSE)
    )
    expect_error(
        calculateUnifrac(tse, assay.type = 2, tree_name = "phylo",
                         weighted = FALSE)
    )
    expect_error(
        calculateUnifrac(tse, assay.type = TRUE, tree_name = "phylo",
                         weighted = FALSE)
    )
    expect_error(
        calculateUnifrac(tse, assay.type = "counts", tree_name = "test",
                         weighted = FALSE)
    )
    expect_error(
        calculateUnifrac(tse, assay.type = "counts", tree_name = 1,
                         weighted = FALSE)
    )
    expect_error(
        calculateUnifrac(tse, assay.type = "counts", tree_name = TRUE,
                         weighted = "FALSE")
    )
    expect_error(
        calculateUnifrac(tse, assay.type = "counts", tree_name = "phylo",
                         weighted = 1)
    )
    
    data(GlobalPatterns, package="mia")
    tse <- GlobalPatterns
    # Calculate unifrac
    unifrac_mia <- as.matrix(calculateUnifrac(tse, weighted = FALSE))
    unifrac_rbiom <- as.matrix(rbiom::unifrac(assay(tse), weighted = FALSE,
                                              rowTree(tse)))
    expect_equal(unifrac_mia, unifrac_rbiom)
    # Calculate unifrac
    unifrac_mia <- as.matrix(calculateUnifrac(tse, weighted = TRUE))
    unifrac_rbiom <- as.matrix(rbiom::unifrac(assay(tse), weighted = TRUE,
                                              rowTree(tse)))
    expect_equal(unifrac_mia, unifrac_rbiom)
    
    # Test with merged object with multiple trees
    tse <- mergeSEs(GlobalPatterns, esophagus, assay.type="counts", missing_values = 0)
    
    # Calculate unifrac
    unifrac_mia <- as.matrix(calculateUnifrac(tse, weighted = FALSE))
    unifrac_rbiom <- as.matrix(rbiom::unifrac(assay(tse), weighted = FALSE,
                                              rowTree(tse)))
    expect_equal(unifrac_mia, unifrac_rbiom)
    # Calculate unifrac
    unifrac_mia <- as.matrix(calculateUnifrac(tse, weighted = TRUE))
    unifrac_rbiom <- as.matrix(rbiom::unifrac(assay(tse), weighted = TRUE,
                                              rowTree(tse)))
    expect_equal(unifrac_mia, unifrac_rbiom)
    
    # Test the function with agglomerated data
    tse <- agglomerateByRank(tse, rank = "Phylum")
    # Calculate unifrac
    unifrac_mia <- as.matrix(calculateUnifrac(tse, weighted = FALSE))
    unifrac_rbiom <- as.matrix(rbiom::unifrac(assay(tse), weighted = FALSE,
                                              rowTree(tse)))
    expect_equal(unifrac_mia, unifrac_rbiom)
})
