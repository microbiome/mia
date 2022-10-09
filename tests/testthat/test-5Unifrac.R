context("Unifrac beta diversity")
test_that("Unifrac beta diversity", {
    
    data(esophagus)
    
    tse <- esophagus
    tse <- relAbundanceCounts(tse)
    
    expect_error(
        calculateUnifrac(tse, assay_name = "test", tree_name = "phylo",
                         weighted = FALSE, normalized = TRUE,
                         BPPARAM = SerialParam())
    )
    expect_error(
        calculateUnifrac(tse, assay_name = 2, tree_name = "phylo",
                         weighted = FALSE, normalized = TRUE,
                         BPPARAM = SerialParam())
    )
    expect_error(
        calculateUnifrac(tse, assay_name = TRUE, tree_name = "phylo",
                         weighted = FALSE, normalized = TRUE,
                         BPPARAM = SerialParam())
    )
    expect_error(
        calculateUnifrac(tse, assay_name = "counts", tree_name = "test",
                         weighted = FALSE, normalized = TRUE,
                         BPPARAM = SerialParam())
    )
    expect_error(
        calculateUnifrac(tse, assay_name = "counts", tree_name = 1,
                         weighted = FALSE, normalized = TRUE,
                         BPPARAM = SerialParam())
    )
    expect_error(
        calculateUnifrac(tse, assay_name = "counts", tree_name = TRUE,
                         weighted = "FALSE", normalized = TRUE,
                         BPPARAM = SerialParam())
    )
    expect_error(
        calculateUnifrac(tse, assay_name = "counts", tree_name = "phylo",
                         weighted = 1, normalized = TRUE,
                         BPPARAM = SerialParam())
    )
    expect_error(
        calculateUnifrac(tse, assay_name = "counts", tree_name = "phylo",
                         weighted = FALSE, normalized = "TRUE",
                         BPPARAM = SerialParam())
    )
    expect_error(
        calculateUnifrac(tse, assay_name = "counts", tree_name = "phylo",
                         weighted = FALSE, normalized = 1,
                         BPPARAM = SerialParam())
    )
    
    data("GlobalPatterns")
    tse <- GlobalPatterns
    # Compare to phyloseq function if it is installed
    if( !require("phyloseq") ){
        # Convert data into phyloseq
        pseq <- makePhyloseqFromTreeSE(tse)
        # Calculate unifrac
        unifrac_tse <- as.matrix(calculateUnifrac(tse))
        unifrac_pseq <- as.matrix(UniFrac(pseq))
        expect_equal(unifrac_tse, unifrac_pseq)
        # Calculate unifrac
        unifrac_tse <- as.matrix(calculateUnifrac(tse, weighted = TRUE, normalized = FALSE))
        unifrac_pseq <- as.matrix(UniFrac(pseq, weighted = TRUE, normalized = FALSE))
        expect_equal(unifrac_tse, unifrac_pseq)
        
        # Test with merged object with multiple trees
        tse <- mergeSEs(GlobalPatterns, esophagus)
        
        # Convert data into phyloseq
        pseq <- makePhyloseqFromTreeSE(tse)
        # Calculate unifrac
        unifrac_tse <- as.matrix(calculateUnifrac(tse))
        unifrac_pseq <- as.matrix(UniFrac(pseq))
        expect_equal(unifrac_tse, unifrac_pseq)
        # Calculate unifrac
        unifrac_tse <- as.matrix(calculateUnifrac(tse, weighted = TRUE, normalized = FALSE, 
                                                  tree_name = "phylo.1"))
        unifrac_pseq <- as.matrix(UniFrac(pseq, weighted = TRUE, normalized = FALSE))
        expect_equal(unifrac_tse, unifrac_pseq)
    }
})
