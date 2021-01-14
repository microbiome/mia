context("estimateDominance")

test_that("estimateDominance", {

    # .simpson_dominance
    #Tests function with input value 1, output should be 1.
    x <- 1
    expect_equal(mia:::.simpson_dominance(x), 1)
    #Tests function with input vector, output should be 0.2048374906.
    x <- c(1,5,2,3,0,10,3,20,3,16)
    expect_equal(mia:::.simpson_dominance(x), 0.2048374906)

    #.gini_dominance
    #reldist package needed. If it is not installed, skips tests tat require it.
    skip_if_not_installed("reldist")
    #Tests function with vector that has value 1 1000 times, output should be 0.
    x <- c(rep(1,1000))
    expect_equal(mia:::.gini_dominance(x), reldist:::gini(x))
    #Tests function with vector that has value 9 one time and value 0 999 times, output should be 0.999.
    x <- c(9,rep(0,999))
    expect_equal(mia:::.gini_dominance(x), reldist:::gini(x))
    #Tests function with vector that has value 1 500 times and value 0 500 times, output should be 0.5.
    x <- c(rep(0,500),rep(1,500))
    expect_equal(mia:::.gini_dominance(x), reldist:::gini(x))
    #Tests function with vector that has values 1,2,3,4,5,6,7,8,9, output should be 0.2962963.
    x <- c(1:9)
    expect_equal(mia:::.gini_dominance(x), reldist:::gini(x))
    #Tests function with vector that has values 1,0,6,1000,2,4739,26,16,10,35,5,28, output should be 0.8738355.
    x <- c(1,0,6,1000,2,4739,26,16,10,35,5,28)
    expect_equal(mia:::.gini_dominance(x), reldist:::gini(x))

    #.get_core_dominance
    #Tests function with microbiome package's core_abundance function
    data("esophagus")
    tse <- esophagus
    data("esophagus", package = "phyloseq")
    phy <- esophagus
    expect_equal(mia:::.get_core_dominance(tse), microbiome:::core_abundance(phy))

    #.get_dominance
    #Tests function with microbiome's dominance_help function. Test all 4 indices with
    #vector of 1000 random decimal numbers.
    x <- runif(1000)
    expect_equal(mia:::.get_dominance(x, index="absolute"), microbiome:::dominance_help(x, index="absolute"))
    x <- runif(1000)
    expect_equal(mia:::.get_dominance(x, index="relative"), microbiome:::dominance_help(x, index="relative"))
    x <- runif(1000)
    expect_equal(mia:::.get_dominance(x, index="DBP"), microbiome:::dominance_help(x, index="dbp"))
    x <- runif(1000)
    expect_equal(mia:::.get_dominance(x, index="DMN"), microbiome:::dominance_help(x, index="dmn"))

    #estimateDominance
    #Calculates all indices.
    tse_idx <- estimateDominance(tse)
    #Checks that the type of output is the same as the type of input.
    expect_true(typeof(tse_idx) == typeof(tse))
    #Checks that every index is calculated by checking the column names from colData.
    #Checks also, that the order of indices is right / the same as the order in the input vector.
    expect_named(colData(tse_idx), c("DBP", "DMN", "absolute", "relative", "simpson", "core_abundance", "gini"))

})
