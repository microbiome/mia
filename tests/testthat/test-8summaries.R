context("summary")

test_that("summary", {
    
    data(GlobalPatterns, package="mia")
    sumdf <- summary(GlobalPatterns, assay.type="counts")
    samples.sum <- sumdf$samples
    
    # check samples
    colnames.samples <- c("total_counts", "min_counts",
                          "max_counts", "median_counts",
                          "mean_counts" ,"stdev_counts")
    expect_equal(colnames(samples.sum),
                 colnames.samples)
    sample.exp.vals <- c(28216678.0, 58688.0, 2357181.0,
                         1106849.0, 1085256.8, 650145.3)
    expect_equal(ceiling(as.numeric(samples.sum[1,])),
                 ceiling(sample.exp.vals))
    # check features
    feature.sum <- sumdf$features
    colnames.feature <- c("total", "singletons",
                          "per_sample_avg")
    expect_equal(colnames(feature.sum),
                 colnames.feature)
    feature.exp.vals <- c(19216.000,2134.000, 4022.231)
    expect_equal(ceiling(as.numeric(feature.sum[1,])),
                 ceiling(feature.exp.vals))
})


context("getUnique")

test_that("getUnique", {
    
    data(GlobalPatterns, package="mia")
    exp.phy <- c("Crenarchaeota","Euryarchaeota",
                 "Actinobacteria","Spirochaetes","MVP-15")
    
    expect_equal(getUnique(GlobalPatterns, "Phylum")[1:5],
                 exp.phy)
})

context("summaries")

test_that("summaries", {
    
    data(GlobalPatterns, package="mia")
    expect_equal( getTop(GlobalPatterns, 
                             method = "mean",
                             top = 5,
                             assay.type = "counts"), 
                  getTop(GlobalPatterns, 
                                 method = "mean",
                                 top = 5,
                                 assay.type = "counts") )
    expect_equal( getUnique(GlobalPatterns, "Phylum", sort = TRUE),
                  getUnique(GlobalPatterns, "Phylum", sort = TRUE) )
    expect_equal( countDominantTaxa(GlobalPatterns),
                  countDominantFeatures(GlobalPatterns))
    
    # Test with multiple equal dominant taxa in one sample
    assay(GlobalPatterns)[1, 1] <- max(assay(GlobalPatterns)[, 1])
    expect_warning(countDominantFeatures(GlobalPatterns, complete = FALSE))
})
