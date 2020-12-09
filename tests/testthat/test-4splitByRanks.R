context("splitByRanks")
test_that("splitByRanks", {
    data(GlobalPatterns)
    x <- GlobalPatterns

    # splitByRanks
    altExps(x) <- splitByRanks(x)
    expect_equal(altExpNames(x),
                 taxonomyRanks(x))
    altExp(x,"Kingdom")
    expect_equal(dim(altExp(x,"Kingdom")),c(2,26))
    expect_equal(dim(altExp(x,"Species")),c(944,26))

    # unsplitByRanks
    x2 <- unsplitByRanks(x)
    expect_equal(dim(altExp(x,"Species")),c(944,26))
    expect_equal(levels(rowData(x2)$taxonomicLevel), taxonomyRanks(x))
})
