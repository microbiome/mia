context("subsampleCounts")
test_that("subsampleCounts", {
  
  data(GlobalPatterns)
  
  tse.subsampled <- subsampleCounts(GlobalPatterns, 
                                    min_size = 60000, 
                                    name = "subsampled",
                                    replace = TRUE,
                                    seed = 1938,
                                    return_type = "TreeSummarizedExperiment")
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
  expect_equal(colSums2(assay(tse.subsampled, "subsampled")), expColSums)
  
  # When replace = FALSE
  tse.subsampled.rp <- subsampleCounts(GlobalPatterns, 
                                       min_size = 60000, 
                                       name = "subsampled",
                                       replace = FALSE,
                                       seed = 1938,
                                       return_type = "TreeSummarizedExperiment")
  
  # check number of features removed is correct
  expnFeaturesRemovedRp <- 6731
  obsnFeaturesRemovedRp <- nrow(GlobalPatterns) - nrow(tse.subsampled.rp)
  expect_equal(obsnFeaturesRemovedRp, expnFeaturesRemovedRp)
  
  # check if all samples subsampled to even depth
  expColSumsRp <- rep(60000, 25)
  expect_equal(colSums2(assay(tse.subsampled.rp, "subsampled")), expColSumsRp)
  
  # check if same Features removed
  obsFeaturesRemovedRp <- rownames(GlobalPatterns)[!rownames(GlobalPatterns) %in% rownames(tse.subsampled.rp)]
  
  expFeaturesRemovedRP <- c("522457","951","586076","244960","215972",
                            "31759","30678","138353","406058","1126")
  
  expect_equal(obsFeaturesRemovedRp[1:10], expFeaturesRemovedRP)
  
  # Check if altExp 
  tse.altExp <- subsampleCounts(GlobalPatterns, 
                                min_size = 5000, 
                                name = "subsampled",
                                replace = TRUE,
                                seed = 1938,
                                return_type = "TreeSummarizedExperiment")
  expect_equal(altExpNames(tse.altExp),
               "subsampled")
  expect_equal(dim(altExp(tse.altExp,"subsampled")),c(6730,26))
  
  # MultiAssay
  mae <- subsampleCounts(GlobalPatterns, 
                         min_size = 60000, 
                         name = "subsampled",
                         replace = TRUE,
                         seed = 1938,
                         return_type = "MultiAssayExperiment")
  expect_s4_class(mae, "MultiAssayExperiment")
  expect_equal(names(mae), c("inputTreeSE","subsampledTreeSE"))
  expect_equal(dim(experiments(mae)$subsampledTreeSE),c(12403,25))
  # check if original treeSE is stored correctly.
  expect_equal(dim(experiments(mae)$inputTreeSE),dim(GlobalPatterns))
  
})
