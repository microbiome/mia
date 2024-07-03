context("rarefyAssay")
test_that("rarefyAssay", {
    seed = 1938
    set.seed(seed)
    data(GlobalPatterns, package="mia")
    
    tse.subsampled <- rarefyAssay(
        GlobalPatterns, 
        sample = 60000,
        name = "subsampled",
        replace = TRUE)
    # check class 
    expect_s4_class(tse.subsampled, "TreeSummarizedExperiment")
    expect_equal(nrow(tse.subsampled), 12403)
    expect_equal(ncol(tse.subsampled), 25)
    # check number of features removed is correct
    expnFeaturesRemoved <- 6813
    obsnFeaturesRemoved <- nrow(GlobalPatterns) - nrow(tse.subsampled)
    expect_equal(obsnFeaturesRemoved, expnFeaturesRemoved)
    
    # check if same Features removed
    obsFeaturesRemoved <- rownames(GlobalPatterns)[!rownames(GlobalPatterns) %in% rownames(tse.subsampled)]
    
    expFeaturesRemoved <- c("951","244423","586076","246140","143239",
                            "244960","144887","141782","215972","31759")
    
    expect_equal(obsFeaturesRemoved[1:10], expFeaturesRemoved)
    
    # check which sample is removed
    expSampleRemoved <- "TRRsed1"
    obsSampleRemoved <- colnames(GlobalPatterns)[!colnames(GlobalPatterns) %in% colnames(tse.subsampled)]
    expect_equal(obsSampleRemoved, expSampleRemoved)
    
    # check if all samples subsampled to even depth
    expColSums <- rep(60000, 25)
    expect_equal(unname(colSums2(assay(tse.subsampled, "subsampled"))), expColSums)
    
    # When replace = FALSE
    seed = 1938
    set.seed(seed)
    tse.subsampled.rp <- rarefyAssay(
        GlobalPatterns, 
        sample = 60000, 
        name = "subsampled",
        replace = TRUE)
    
    # check number of features removed is correct
    expnFeaturesRemovedRp <- 6731
    obsnFeaturesRemovedRp <- nrow(GlobalPatterns) - nrow(tse.subsampled.rp)
    expect_equal(obsnFeaturesRemovedRp, expnFeaturesRemovedRp)
    
    # check if all samples subsampled to even depth
    expColSumsRp <- rep(60000, 25)
    expect_equal(unname(colSums2(assay(tse.subsampled.rp, "subsampled"))), expColSumsRp)
    
    # check if same Features removed
    obsFeaturesRemovedRp <- rownames(GlobalPatterns)[!rownames(GlobalPatterns) %in% rownames(tse.subsampled.rp)]
    
    expFeaturesRemovedRP <- c(
        "951","244423","586076","246140","143239", "31759","30678","138353",
        "406058","1126")
    
    expect_equal(obsFeaturesRemovedRp[1:10], expFeaturesRemovedRP)
})
