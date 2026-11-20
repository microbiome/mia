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
    actual2 <- getNMDS(se,nmds.fun="monoMDS",pc=FALSE,scaling=FALSE)
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

test_that("getNMDS uses getDissimilarity by default", {
    skip_if_not(requireNamespace("vegan", quietly = TRUE))
    data(esophagus, package = "mia")
    tse <- esophagus
    # The default gives the same result as vegdist for Bray-Curtis
    set.seed(1)
    res_default <- getNMDS(tse, method = "bray")
    set.seed(1)
    res_vegdist <- getNMDS(tse, method = "bray", FUN = vegan::vegdist)
    expect_equal(res_default, res_vegdist)
    # Function names are accepted
    set.seed(1)
    res_name <- getNMDS(tse, method = "bray", FUN = "getDissimilarity")
    expect_equal(res_default, res_name)
    # Bray-Curtis is used when method is not given, as with vegdist before
    set.seed(1)
    expect_equal(getNMDS(tse), res_default)
    # Dissimilarities that only getDissimilarity supports work by default
    res_unifrac <- getNMDS(
        tse, method = "unifrac", tree = rowTree(tse), keep.dist = TRUE)
    expect_equal(
        unname(as.matrix(attr(res_unifrac, "dist"))),
        unname(as.matrix(getDissimilarity(tse, method = "unifrac"))))
})
