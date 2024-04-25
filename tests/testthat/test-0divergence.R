context("divergence estimates")
test_that("divergence estimates", {
  
  skip_if_not(requireNamespace("vegan", quietly = TRUE))
  
  data(esophagus, package="mia")
  tse <- esophagus
  
  tse_divergence <- addDivergence(tse)
  
  # Checks that the type of output is the same as the type of input.
  expect_true(typeof(tse_divergence) == typeof(tse))
  
  # Check that colData includes divergence results
  expect_true( "divergence" %in% colnames(colData(tse_divergence)))
  
  # Expect errors when input is wrong
  expect_error(addDivergence(tse, name = 123, reference = "median", 
                                   FUN = vegan::vegdist, method = "euclidean") )
  expect_error(addDivergence(tse, name = "test", reference = "test", 
                                   FUN = vegan::vegdist, method = "euclidean") )
  expect_error(addDivergence(tse, name = "test", reference = "median", 
                                   FUN = "test", method = "euclidean") )
  expect_error(addDivergence(tse, name = "test", reference = "median", 
                                   FUN = vegan::vegdist, method = "test") )
  expect_error(addDivergence(tse, reference = rep(0, nrow(tse)), 
                                  FUN = "test",
                                  method = "euclidean"))
  expect_error(addDivergence(tse, 
                                  reference = rep(0, nrow(tse)),
                                  FUN = stats::dist,
                                  method = "test"))
  expect_error(addDivergence(tse, 
                                  reference = "test",
                                  FUN = stats::dist,
                                  method = "euclidean"))
  expect_error(addDivergence(tse, 
                                  reference = rep("test", nrow(tse)),
                                  FUN = stats::dist,
                                  method = "euclidean"))
  expect_error(addDivergence(tse, 
                                  reference = rep(0, nrow(tse)-1),
                                  FUN = stats::dist,
                                  method = "euclidean"))
  expect_error(addDivergence(tse, 
                                  reference = rep("test", nrow(tse)-1),
                                  FUN = stats::dist,
                                  method = "euclidean"))
  expect_error(addDivergence(tse, 
                                  reference = rep(TRUE, nrow(tse)),
                                  FUN = stats::dist,
                                  method = "euclidean"))
  expect_error(addDivergence(tse, 
                                  reference = FALSE,
                                  FUN = stats::dist,
                                  method = "euclidean"))
  
  # Reference values from microbiome pkg's divergence function
  expect_equal(unname(round(colData(
      addDivergence(tse, 
                       reference = "mean",
                       FUN = stats::dist,
                       method = "euclidean"))$divergence, 6)),
    round(c(35.35534, 42.16634, 59.44746)),6)
  
  expect_equal(unname(round(colData(
      addDivergence(tse, 
                       reference = assay(tse, "counts")[,3],
                       FUN = stats::dist,
                       method = "manhattan"))$divergence, 6)),
    round(c(210, 280, 0)),6)
  
  expect_equal(unname(round(colData(
      addDivergence(tse, 
                       reference = assay(tse, "counts")[,1],
                       FUN = vegan::vegdist,
                       method = "chao"))$divergence, 6)),
    round(c(0.00000000, 0.10115766, 0.08239422)),6)

})
