test_that("Estimate Alpha Diversity Indices with Rarefaction", {
    data(GlobalPatterns, package="mia")
    tse <- GlobalPatterns
    ## Testing diversity
    # Calculate the default Shannon index with no rarefaction
    tse <- estimateAlpha(tse, assay.type = "counts", index = "shannon")
    expect_true(any(grepl("shannon_diversity", colnames(colData(tse)))))
    # Calculate same index with 10 rarefaction rounds
    tse <- estimateAlpha(tse, assay.type = "counts", index = "shannon",
                         rarify=TRUE, nrounds=10, name="shannon_10")
    expect_true(any(grepl("shannon_10", colnames(colData(tse)))))
    # comparing the estimates
    expect_false(any(colData(tse)$shannon_diversity==colData(tse)$shannon_10))
    
    ## Testing Dominance
    # Calculate the default gini_dominance index with no rarefaction
    tse <- estimateAlpha(tse, assay.type = "counts", index = "gini_dominance")
    expect_true(any(grepl("gini_dominance", colnames(colData(tse)))))
    # Calculate same index with 10 rarefaction rounds
    tse <- estimateAlpha(tse, assay.type = "counts", index = "gini_dominance",
                         rarify=TRUE, nrounds=10, name="gini_dominance_10")
    expect_true(any(grepl("gini_dominance_10", colnames(colData(tse)))))
    # comparing the estimates
    expect_false(any(colData(tse)$gini_dominance==colData(tse)$gini_dominance_10))
    
    ## Testing Evenness
    # Calculate the default pielou index with no rarefaction
    tse <- estimateAlpha(tse, assay.type = "counts", index = "pielou")
    expect_true(any(grepl("pielou_evenness", colnames(colData(tse)))))
    # Calculate same index with 10 rarefaction rounds
    tse <- estimateAlpha(tse, assay.type = "counts", index = "pielou",
                         rarify=TRUE, nrounds=10, name="pielou_10")
    expect_true(any(grepl("pielou_10", colnames(colData(tse)))))
    # comparing the estimates
    expect_false(any(colData(tse)$pielou_evenness==colData(tse)$pielou_10))
    
    ## Testing Richness
    # Calculate the default chao1 index with no rarefaction
    tse <- estimateAlpha(tse, assay.type = "counts", index = "chao1")
    expect_true(any(grepl("chao1_richness", colnames(colData(tse)))))
    # Calculate same index with 10 rarefaction rounds
    tse <- estimateAlpha(tse, assay.type = "counts", index = "chao1",
                         rarify=TRUE, nrounds=10, name="chao1_10")
    expect_true(any(grepl("pielou_10", colnames(colData(tse)))))
    # comparing the estimates
    expect_false(any(colData(tse)$chao1_richness==colData(tse)$chao1_10))
    
    # test non existing index
    expect_error(estimateAlpha(tse, assay.type = "counts", index = "ödsaliufg"))

})