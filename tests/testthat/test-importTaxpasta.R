skip_if_not_installed("biomformat")

test_that("a complete BIOM file can be loaded without error", {
  expect_no_error(
    suppressWarnings({
      importTaxpasta(test_path("data", "complete.biom"))
    })
  )
})

test_that("loading a simple BIOM file without taxonomy information generates a warning", {
  expect_warning(importTaxpasta(test_path("data", "simple.biom")))
})

test_that("loading a BIOM file with less ranks than taxonomy columns generates an error", {
  expect_error(importTaxpasta(test_path("data", "less_ranks.biom")))
})

test_that("loading a BIOM file with more ranks than taxonomy columns generates an error", {
  expect_error(importTaxpasta(test_path("data", "more_ranks.biom")))
})

result <- suppressWarnings({
  importTaxpasta(test_path("data", "complete.biom"))
})

test_that("loading a taxpasta BIOM file returns a TreeSummarizedExperiment", {
  expect_s4_class(result, "TreeSummarizedExperiment")
})

requireNamespace("ape", quietly = TRUE)
requireNamespace("TreeSummarizedExperiment", quietly = TRUE)

test_that("the TreeSummarizedExperiment has expected dimensions", {
  expect_identical(dim(result), c(6L, 2L))
})

test_that("the TreeSummarizedExperiment has expected taxonomic ranks", {
  expect_identical(taxonomyRanks(result), c("Superkingdom", "Clade", "Class", "Order", "Family"))
})

test_that("the TreeSummarizedExperiment has one alternative experiment per taxonomic rank", {
  expect_identical(SingleCellExperiment::altExpNames(result), c("Superkingdom", "Clade", "Class", "Order", "Family"))
})

test_that("the TreeSummarizedExperiment has one tree with two leaves", {
  tree <- TreeSummarizedExperiment::rowTree(result)
  expect_identical(ape::Ntip(tree), 2L)
  expect_setequal(tree$tip.label, c("Family:Lachnospiraceae", "Family:Vallitaleaceae"))
})

test_that("the TreeSummarizedExperiment has one tree with four inner nodes", {
  tree <- TreeSummarizedExperiment::rowTree(result)
  expect_identical(ape::Nnode(tree), 4L)
})

test_that("the TreeSummarizedExperiment has one tree with five edges", {
  tree <- TreeSummarizedExperiment::rowTree(result)
  expect_identical(ape::Nedge(tree), 5L)
})
