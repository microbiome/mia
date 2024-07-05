skip_if_not_installed("biomformat")

# Get default taxonomy ranks (importTaxpasta sets ranks based on data)
default_ranks <- getTaxonomyRanks()

test_that("a complete BIOM file can be loaded without error", {
  expect_no_error(
    suppressWarnings({
      importTaxpasta(system.file("extdata/complete.biom", package = "mia", mustWork = TRUE))
    })
  )
})

test_that("loading a simple BIOM file without taxonomy information generates a warning", {
  expect_warning(importTaxpasta(system.file("extdata/simple.biom", package = "mia", mustWork = TRUE)))
})

test_that("loading a BIOM file with less ranks than taxonomy columns generates an error", {
  expect_error(importTaxpasta(system.file("extdata/less_ranks.biom", package = "mia", mustWork = TRUE)))
})

test_that("loading a BIOM file with more ranks than taxonomy columns generates an error", {
  expect_error(importTaxpasta(system.file("extdata/more_ranks.biom", package = "mia", mustWork = TRUE)))
})

result <- suppressWarnings({
  importTaxpasta(system.file("extdata/complete.biom", package = "mia", mustWork = TRUE))
})

test_that("loading a taxpasta BIOM file returns a TreeSummarizedExperiment", {
  expect_s4_class(result, "TreeSummarizedExperiment")
})

test_that("the TreeSummarizedExperiment has expected dimensions", {
  expect_identical(dim(result), c(6L, 2L))
})

test_that("the TreeSummarizedExperiment has expected taxonomic ranks", {
  expect_identical(taxonomyRanks(result), c("Superkingdom", "Clade", "Class", "Order", "Family"))
})

test_that("the TreeSummarizedExperiment has one alternative experiment per taxonomic rank", {
  expect_identical(altExpNames(result), c("Superkingdom", "Clade", "Class", "Order", "Family"))
})

test_that("the TreeSummarizedExperiment has one tree with two leaves", {
  tree <- rowTree(result)
  expect_identical(ape::Ntip(tree), 2L)
  expect_setequal(tree$tip.label, c("Family:Lachnospiraceae", "Family:Vallitaleaceae"))
})

test_that("the TreeSummarizedExperiment has one tree with four inner nodes", {
  tree <- rowTree(result)
  expect_identical(ape::Nnode(tree), 4L)
})

test_that("the TreeSummarizedExperiment has one tree with five edges", {
  tree <- rowTree(result)
  expect_identical(ape::Nedge(tree), 5L)
})

# Set taxonomy ranks to default (importTaxpasta sets ranks based on data)
setTaxonomyRanks(default_ranks)
