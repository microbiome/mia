context("addLDA")
test_that("addLDA", {
  skip_if_not(require("topicmodels", quietly = TRUE))
  data(GlobalPatterns, package="mia")
  #
  tse <- addLDA(GlobalPatterns)
  expect_named(reducedDims(tse),"LDA")
  expect_true(is.matrix(reducedDim(tse,"LDA")))
  expect_equal(dim(reducedDim(tse,"LDA")),c(26,2))
  red <- reducedDim(tse,"LDA")
  expect_equal(names(attributes(red)),
               c("dim","dimnames","loadings"))
  expect_equal(dim(attr(red,"loadings")),c(19216,2))
  
  # ERRORs
  expect_error(
    addLDA(GlobalPatterns, k = "test", assay.type = "counts", name = "LDA")
  )
  expect_error(
    addLDA(GlobalPatterns, k = 1.5, assay.type = "counts", name = "LDA")
  )
  expect_error(
    addLDA(GlobalPatterns, k = TRUE, assay.type = "counts", name = "LDA")
  )
  expect_error(
    addLDA(GlobalPatterns, k = 2, assay.type = "test", name = "LDA")
  )
  expect_error(
    addLDA(GlobalPatterns, k = 2, assay.type = 1, name = "LDA")
  )
  expect_error(
    addLDA(GlobalPatterns, k = 2, assay.type = TRUE, name = "LDA")
  )
  expect_error(
    addLDA(GlobalPatterns, k = 2, assay.type = "counts", name = 1)
  )
  expect_error(
    addLDA(GlobalPatterns, k = 2, assay.type = "counts", name = TRUE)
  )
})
