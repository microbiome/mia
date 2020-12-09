context("prevalence")

test_that("getPrevalence", {

    data(GlobalPatterns)
    expect_error(getPrevalence(GlobalPatterns, detection="test"),
                 "'detection' must be a single numeric value or coercible to one")
    expect_error(getPrevalence(GlobalPatterns, include_lowest="test"),
                 "'include_lowest' must be TRUE or FALSE")
    expect_error(getPrevalence(GlobalPatterns, sort="test"),
                 "'sort' must be TRUE or FALSE")
    expect_error(getPrevalence(GlobalPatterns, as_relative="test"),
                 "'as_relative' must be TRUE or FALSE")
    expect_error(getPrevalence(GlobalPatterns, abund_values="test"),
                 "'abund_values' must be a valid name")
    # Output should be always a frequency between 0 to 1
    pr <- getPrevalence(GlobalPatterns, detection=0.1/100, as_relative=TRUE)
    expect_true(min(pr) >= 0 && max(pr) <= 1)
    pr <- getPrevalence(GlobalPatterns, detection=0.1/100, as_relative=FALSE)
    expect_true(min(pr) >= 0 && max(pr) <= 1)

    # Same prevalences should be returned for as_relative T/F in certain cases.
    pr1 <- getPrevalence(GlobalPatterns, detection=1, include_lowest=TRUE, as_relative=FALSE)
    pr2 <- getPrevalence(GlobalPatterns, detection=0/100, include_lowest=FALSE, as_relative=TRUE)
    expect_true(all(pr1 == pr2))

    # Same prevalences should be returned for as_relative T/F in certain cases.
    pr1 <- getPrevalence(GlobalPatterns, detection=1, include_lowest=TRUE, as_relative=FALSE)
    pr2 <- getPrevalence(GlobalPatterns, detection=0, include_lowest=FALSE, as_relative=FALSE)
    expect_true(all(pr1 == pr2))

    # Different ways to use relative abundance should yield the same output
    pr2 <- getPrevalence(GlobalPatterns, as_relative=TRUE, abund_values = "counts")
    GlobalPatterns <- relAbundanceCounts(GlobalPatterns)
    pr1 <- getPrevalence(GlobalPatterns, as_relative=FALSE, abund_values = "relabundance")
    expect_true(all(pr1 == pr2))

    # Sorting should put the top values first
    pr <- getPrevalence(GlobalPatterns, sort=TRUE, detection = 0.1/100)
    expect_equal(as.vector(which.max(pr)), 1)
    pr <- names(head(getPrevalence(GlobalPatterns, sort=TRUE,  include_lowest = TRUE), 5L))
    actual <- getTopTaxa(GlobalPatterns,
                         method="prevalence",
                         top=5,
                         abund_values="counts")
    expect_equal(pr, actual)

})


test_that("getPrevalentTaxa", {

    data(GlobalPatterns)
    expect_error(getPrevalentTaxa(GlobalPatterns, prevalence="test"),
                 "'prevalence' must be a single numeric value or coercible to one")
    # Results compatible with getPrevalence
    pr1 <- getPrevalentTaxa(GlobalPatterns, detection=0.1/100, as_relative=TRUE, sort=TRUE)
    pr2 <- names(getPrevalence(GlobalPatterns, rank = "Kingdom", detection=0.1/100, as_relative=TRUE, sort=TRUE))
    expect_true(all(pr1 == pr2))

    # Same sorting for toptaxa obtained in different ways
    pr1 <- getPrevalentTaxa(GlobalPatterns, detection=0.1/100, as_relative=TRUE, sort=TRUE)
    pr2 <- names(getPrevalence(GlobalPatterns, rank = "Kingdom", detection=0.1/100, as_relative=TRUE, sort=TRUE))
    expect_true(all(pr1 == pr2))

    # Retrieved taxa are the same for counts and relative abundances
    pr1 <- getPrevalentTaxa(GlobalPatterns, detection=0.1/100, as_relative=TRUE)
    pr2 <- getPrevalentTaxa(GlobalPatterns, detection=0.1/100, as_relative=FALSE)
    expect_true(all(pr1 == pr2))

    # Prevalence and detection threshold at 0 has the same impact on counts and relative abundances
    pr1 <- getPrevalentTaxa(GlobalPatterns, detection=0, prevalence=0, as_relative=TRUE)
    pr2 <- getPrevalentTaxa(GlobalPatterns, detection=0, prevalence=0, as_relative=FALSE)
    expect_true(all(pr1 == pr2))

})
