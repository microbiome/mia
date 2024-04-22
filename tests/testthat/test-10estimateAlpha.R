test_that("Estimate Alpha Diversity Indices with Rarefaction", {
    data(GlobalPatterns, package="mia")
    tse <- GlobalPatterns
    ## Testing diversity
    # Calculate the default Shannon index with no rarefaction
    tse <- addAlpha(tse, assay.type = "counts", index = "shannon")
    expect_true(any(grepl("shannon", colnames(colData(tse)))))
    tse <- addAlpha(tse, assay.type = "counts", index = "shannon_diversity")
    expect_true(any(grepl("shannon_diversity", colnames(colData(tse)))))
    # Calculate same index with 10 rarefaction rounds
    tse <- addAlpha(
        tse, assay.type = "counts", index = "shannon",
        rarefaction.depth = min(colSums(assay(tse, "counts")), na.rm = TRUE),
        n.iter = 10, name = "shannon_10")
    expect_true(any(grepl("shannon_10", colnames(colData(tse)))))
    # comparing the estimates
    expect_false( all(tse$shannon_diversity == tse$shannon_10) )
    
    ## Testing Dominance
    # Calculate the default gini_dominance index with no rarefaction
    tse <- addAlpha(tse, assay.type = "counts", index = "gini_dominance")
    expect_true( any(grepl("gini_dominance", colnames(colData(tse)))) )
    # Calculate same index with 10 rarefaction rounds
    tse <- addAlpha(
        tse, assay.type = "counts", index = "gini_dominance",
        rarefaction.depth = min(colSums(assay(tse, "counts")), na.rm = TRUE),
        n.iter = 10, name = "gini_dominance_10")
    expect_true(any(grepl("gini_dominance_10", colnames(colData(tse)))))
    # comparing the estimates
    expect_false(all(tse$gini_dominance==tse$gini_dominance_10))
    
    ## Testing Evenness
    # Calculate the default pielou index with no rarefaction
    tse <- addAlpha(tse, assay.type = "counts", index = "pielou")
    expect_true(any(grepl("pielou", colnames(colData(tse)))))
    # Calculate same index with 10 rarefaction rounds
    tse <- addAlpha(
        tse, assay.type = "counts", index = "pielou",
        rarefaction.depth = min(colSums(assay(tse, "counts")), na.rm = TRUE),
        n.iter = 10, name = "pielou_10")
    expect_true(any(grepl("pielou_10", colnames(colData(tse)))))
    # comparing the estimates
    expect_false(all(tse$pielou==tse$pielou_10))
    
    ## Testing Richness
    # Calculate the default chao1 index with no rarefaction
    tse <- addAlpha(tse, assay.type = "counts", index = "chao1")
    expect_true(any(grepl("chao1", colnames(colData(tse)))))
    # Calculate same index with 10 rarefaction rounds
    tse <- addAlpha(
        tse, assay.type = "counts", index = "chao1",
        rarefaction.depth=0.1*mean(colSums(assay(tse, "counts")), na.rm = TRUE),
        n.iter = 10, name = "chao1_10")
                         
    expect_true(any(grepl("pielou_10", colnames(colData(tse)))))
    # comparing the estimates
    expect_false(all(tse$chao1==tse$chao1_10))
    
    # test non existing index
    expect_error(addAlpha(tse, assay.type = "counts", index = "ödsaliufg"))
    
    # comparing 10 iter with 20 iters estimates
    tse <- addAlpha(
        tse, assay.type = "counts", index = "shannon",
        rarefaction.depth = min(colSums(assay(tse, "counts")), na.rm = TRUE),
        n.iter=20, name="shannon_20")
    # comparing the estimates
    expect_false(all(tse$shannon_20==tse$shannon_10))
    
    # Testing with multiple indices
    tse <- addAlpha(
        tse, assay.type = "counts",
        index = c("coverage","absolute", "camargo", "ace"))
    expect_true(any(grepl("coverage", colnames(colData(tse)))))
    expect_true(any(grepl("absolute", colnames(colData(tse)))))
    expect_true(any(grepl("camargo", colnames(colData(tse)))))
    expect_true(any(grepl("ace", colnames(colData(tse)))))
    
    # Testing with multiple indices with rarefaction
    tse <- addAlpha(
        tse, assay.type = "counts", 
        rarefaction.depth=min(colSums(assay(tse, "counts")), na.rm = TRUE),
        n.iter = 10, 
        index = c("coverage","absolute", "camargo", "ace"),
        name = c("coverage_10","absolute_10", "camargo_10", "ace_10"))
    expect_true(any(grepl("coverage_10", colnames(colData(tse)))))
    expect_true(any(grepl("absolute_10", colnames(colData(tse)))))
    expect_true(any(grepl("camargo_10", colnames(colData(tse)))))
    expect_true(any(grepl("ace_10", colnames(colData(tse)))))
    expect_false(all(tse$coverage_==tse$coverage_10))
    expect_false(all(tse$absolute==tse$absolute_10))
    expect_false(all(tse$camargo==tse$camargo_10))
    expect_false(all(tse$ace==tse$ace_10))
})
