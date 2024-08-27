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
  tse <- GlobalPatterns
  # Test if Unifrac dissimilarity works with external tree
  expect_equal(getDissimilarity(tse, method = "unifrac"),
               getDissimilarity(tse, tree = rowTree(tse), method = "unifrac"))
  # Test rarefaction
  ntop <- 5
  tse_sub <- tse[head(rev(order(rowSds(assay(tse, "counts")))), ntop), ]
  mat <- assay(tse_sub, "counts")
  clr <- function (x) {
    vegan::decostand(x, method="clr", pseudocount=1)
  }
  set.seed(123)
  res1 <- vegan::avgdist(t(mat), distfun = vegdist, dmethod = "euclidean",
                         sample = min(colSums2(mat)), 
                         iterations = 10, transf = clr)
  set.seed(123)
  res2 <- getDissimilarity(tse_sub, method = "euclidean", niter = 10, 
                           transf = clr)
  expect_equal(as.matrix(res1), as.matrix(res2))
})
