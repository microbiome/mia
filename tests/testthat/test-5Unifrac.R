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
})
