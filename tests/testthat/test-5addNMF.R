context("addNMF")
test_that("addNMF", {
  skip_if_not_installed("NMF")
  data(GlobalPatterns, package="mia")
  #
  tse <- GlobalPatterns
  tse <- agglomerateByPrevalence(tse, rank = "Phylum")
  tse <- addNMF(tse)
  expect_named(reducedDims(tse),"NMF")
  expect_true(is.matrix(reducedDim(tse,"NMF")))
  expect_equal(dim(reducedDim(tse,"NMF")),c(26,2))
  red <- reducedDim(tse,"NMF")
  expect_equal(names(attributes(red)),
               c("dim","dimnames","loadings", "model"))
  expect_equal(dim(attr(red,"loadings")),c(35,2))
  # Check if ordination matrix returned by NMF::nmf is the same as
  # getNMF and addNMF ones
  mat <- t(assay(tse, "counts"))
  library("NMF")
  nmf_model <- NMF::nmf(mat, 2)
  loadings <- nmf_model@fit@H
  # Compare NMF::nmf and addNMF
  expect_equal(loadings, t(attr(red, "loadings")), tolerance = 10**-4)
  scores2 <- getNMF(tse)
  # Compare NMF::nmf and getNMF
  expect_equal(loadings, attr(scores2, "loadings"), tolerance = 10**-4)
  # ERRORs
  expect_error(
    addNMF(GlobalPatterns, k = "test", assay.type = "counts", name = "NMF")
  )
  expect_error(
    addNMF(GlobalPatterns, k = 1.5, assay.type = "counts", name = "NMF")
  )
  expect_error(
    addNMF(GlobalPatterns, k = TRUE, assay.type = "counts", name = "NMF")
  )
  expect_error(
    addNMF(GlobalPatterns, k = 2, assay.type = "test", name = "NMF")
  )
  expect_error(
    addNMF(GlobalPatterns, k = 2, assay.type = 1, name = "NMF")
  )
  expect_error(
    addNMF(GlobalPatterns, k = 2, assay.type = TRUE, name = "NMF")
  )
  expect_error(
    addNMF(GlobalPatterns, k = 2, assay.type = "counts", name = 1)
  )
  expect_error(
    addNMF(GlobalPatterns, k = 2, assay.type = "counts", name = TRUE)
  )
})
