context("agglomerateByRanks")
test_that("agglomerateByRanks", {
    data(GlobalPatterns, package="mia")
    x <- GlobalPatterns

    # splitByRanks
    x <- agglomerateByRanks(x)
    expect_equal(altExpNames(x),
                 taxonomyRanks(x))
    altExp(x,"Kingdom")
    expect_equal(dim(altExp(x,"Kingdom")),c(2,26))
    expect_equal(dim(altExp(x,"Species")),c(944,26))
    
    # Chceck that agglomerateByRanks and splitByRanks results to same output
    res1 <- altExps(x)
    res2 <- splitByRanks(x)
    expect_equal( names(res1), names(res2) )
    expect_equal( res1[[1]], res2[[1]] )
    expect_equal( res1[[3]], res2[[3]] )

    # unsplitByRanks
    x2 <- unsplitByRanks(x)
    expect_equal(dim(altExp(x,"Species")),c(944,26))
    expect_equal(levels(rowData(x2)$taxonomicLevel), taxonomyRanks(x))
})
