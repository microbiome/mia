context("Dissimilarity calculation")
test_that("Dissimilarity calculation", {
  data(GlobalPatterns, package="mia")
  tse <- GlobalPatterns
  # Test if Unifrac dissimilarity works with external tree
  expect_equal(getDissimilarity(tse, method = "unifrac"),
               getDissimilarity(tse, tree = rowTree(tse), method = "unifrac"))
  # Test if results from vegan::vegdist are equal to getDissimilarity results
  mat <- assay(tse, "counts")
  expect_equal(as.matrix(vegan::vegdist(t(mat), "euclidean")),
               as.matrix(getDissimilarity(tse, method = "euclidean")))
})
