context("addNMDS")
test_that("addNMDS", {
    mat <- matrix(1:60, nrow = 6)
    df <- DataFrame(n = c(1:6))
    se <- SummarizedExperiment(assays = list(counts = mat),
                               rowData = df)
    #
    actual <- getNMDS(se)
    expect_true(is.matrix(actual))
    expect_equal(dim(actual),c(10,2))
    actual2 <- getNMDS(se,nmds="monoMDS",pc=FALSE,scaling=FALSE)
    expect_true(is.matrix(actual))
    expect_equal(dim(actual),c(10,2))
    expect_true(sum(actual2 - actual) < 0.00001)
    #
    data(esophagus, package="mia")
    esophagus <- addNMDS(esophagus, distFUN = vegan::vegdist, name = "BC")
    esophagus <- addNMDS(esophagus, distFUN = vegan::vegdist, name = "euclidean",
                         method = "euclidean")
    expect_named(reducedDims(esophagus),c("BC","euclidean"))
    expect_true(is.matrix(reducedDim(esophagus,"BC")))
    expect_equal(dim(reducedDim(esophagus,"BC")),c(3,2))
    expect_true(is.matrix(reducedDim(esophagus,"BC")))
    expect_equal(dim(reducedDim(esophagus,"BC")),c(3,2))
})
