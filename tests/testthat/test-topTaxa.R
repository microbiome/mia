context("topTaxa")
test_that("topTaxa", {

    data(GlobalPatterns, package = "MicrobiomeExperiment")
    mean.taxa <- c("549656", "331820", "279599", "360229", "317182")
    sum.taxa <- c("549656", "331820", "279599", "360229", "317182")
    median.taxa <- c("549656", "331820", "317182", "94166",  "279599")
    top_mean <- topTaxa(GlobalPatterns, method="mean", top=5, abund_values="counts")
    top_sum <- topTaxa(GlobalPatterns, method="sum", top=5, abund_values="counts")
    top_median <- topTaxa(GlobalPatterns, method="median", top=5, abund_values="counts")

    expect_equal(top_mean, mean.taxa)
    expect_equal(top_sum, sum.taxa)
    expect_equal(top_median, median.taxa)

    expect_error(mia:::.check_max_taxa(GlobalPatterns, 100000000, "counts"),
                 "top value is greater than maximum number of taxa in assay")

})
