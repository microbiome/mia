context("Unifrac beta diversity")
test_that("Unifrac beta diversity", {
    skip_if_not(require("phyloseq", quietly = TRUE))
    
    data(esophagus, package="mia")
    tse <- esophagus
    tse <- transformCounts(tse, assay.type="counts", method="relabundance")
    
    expect_error(
        calculateUnifrac(tse, assay.type = "test", tree_name = "phylo",
                         weighted = FALSE, normalized = TRUE,
                         BPPARAM = SerialParam())
    )
    expect_error(
        calculateUnifrac(tse, assay.type = 2, tree_name = "phylo",
                         weighted = FALSE, normalized = TRUE,
                         BPPARAM = SerialParam())
    )
    expect_error(
        calculateUnifrac(tse, assay.type = TRUE, tree_name = "phylo",
                         weighted = FALSE, normalized = TRUE,
                         BPPARAM = SerialParam())
    )
    expect_error(
        calculateUnifrac(tse, assay.type = "counts", tree_name = "test",
                         weighted = FALSE, normalized = TRUE,
                         BPPARAM = SerialParam())
    )
    expect_error(
        calculateUnifrac(tse, assay.type = "counts", tree_name = 1,
                         weighted = FALSE, normalized = TRUE,
                         BPPARAM = SerialParam())
    )
    expect_error(
        calculateUnifrac(tse, assay.type = "counts", tree_name = TRUE,
                         weighted = "FALSE", normalized = TRUE,
                         BPPARAM = SerialParam())
    )
    expect_error(
        calculateUnifrac(tse, assay.type = "counts", tree_name = "phylo",
                         weighted = 1, normalized = TRUE,
                         BPPARAM = SerialParam())
    )
    expect_error(
        calculateUnifrac(tse, assay.type = "counts", tree_name = "phylo",
                         weighted = FALSE, normalized = "TRUE",
                         BPPARAM = SerialParam())
    )
    expect_error(
        calculateUnifrac(tse, assay.type = "counts", tree_name = "phylo",
                         weighted = FALSE, normalized = 1,
                         BPPARAM = SerialParam())
    )
    
    data(GlobalPatterns, package="mia")
    tse <- GlobalPatterns
    # Compare to phyloseq function
    .require_package("phyloseq")
    # Convert data into phyloseq
    pseq <- makePhyloseqFromTreeSE(tse)
    # Calculate unifrac
    unifrac_tse <- as.matrix(calculateUnifrac(tse))
    unifrac_pseq <- as.matrix(phyloseq::UniFrac(pseq))
    expect_equal(unifrac_tse, unifrac_pseq)
    # Calculate unifrac
    unifrac_tse <- as.matrix(calculateUnifrac(tse, weighted = TRUE, normalized = FALSE))
    unifrac_pseq <- as.matrix(phyloseq::UniFrac(pseq, weighted = TRUE, normalized = FALSE))
    expect_equal(unifrac_tse, unifrac_pseq)
    
    # Test with merged object with multiple trees
    tse <- mergeSEs(GlobalPatterns, esophagus, assay.type="counts", missing_values = 0)
    
    # Convert data into phyloseq
    pseq <- makePhyloseqFromTreeSE(tse, assay.type="counts")
    # Convert back to TreeSE (pseq has pruned tree)
    tse <- makeTreeSEFromPhyloseq(pseq)
    # Calculate unifrac
    unifrac_tse <- as.matrix(calculateUnifrac(tse))
    unifrac_pseq <- as.matrix(phyloseq::UniFrac(pseq))
    expect_equal(unifrac_tse, unifrac_pseq)
    # Calculate unifrac
    unifrac_tse <- as.matrix(calculateUnifrac(tse, weighted = TRUE, normalized = FALSE, 
                                              tree_name = "phylo"))
    unifrac_pseq <- as.matrix(phyloseq::UniFrac(pseq, weighted = TRUE, normalized = FALSE))
    expect_equal(unifrac_tse, unifrac_pseq)
    
    # Test the function with agglomerated data
    tse <- agglomerateByRank(tse, rank = "Phylum")
    # Convert data into phyloseq
    suppressWarnings( pseq <- makePhyloseqFromTreeSE(tse) )
    # Convert back to TreeSE (pseq has pruned tree)
    tse <- makeTreeSEFromPhyloseq(pseq)
    # Calculate unifrac
    unifrac_tse <- as.matrix(calculateUnifrac(tse, normalized = TRUE))
    unifrac_pseq <- as.matrix(phyloseq::UniFrac(pseq, normalized = TRUE))
    expect_equal(unifrac_tse, unifrac_pseq)
    
    # Detach phyloseq package
    unloadNamespace("phyloseq")
})
