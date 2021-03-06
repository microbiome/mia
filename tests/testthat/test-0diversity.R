context("diversity estimates")
test_that("diversity estimates", {

    skip_if_not(requireNamespace("vegan", quietly = TRUE))
    data(esophagus)

    tse <- esophagus
    tse <- relAbundanceCounts(tse)

    tse_idx <- estimateDiversity(tse, threshold = 0.473)

    # Checks that the type of output is the same as the type of input.
    expect_true(typeof(tse_idx) == typeof(tse))

    # Check that every index is calculated by checking the column names from
    # colData.
    # Check that the order of indices is right / the same as the order
    # in the input vector.
    expect_named(colData(tse_idx), c("coverage", "fisher", "gini_simpson",
                                     "inverse_simpson", "log_modulo_skewness",
                                     "shannon", "faith"))

    lambda <- unname(colSums(assays(tse_idx)$relabundance^2))
    ginisimpson <- 1 - lambda
    invsimpson <- 1 / lambda

    expect_equal(lambda, .simpson_lambda(assays(tse_idx)$relabundance))
    expect_equal(ginisimpson, colData(tse_idx)$gini_simpson)
    expect_equal(ginisimpson, .calc_gini_simpson(assays(tse_idx)$relabundance))

    expect_equal(invsimpson, colData(tse_idx)$inverse_simpson)
    expect_equal(invsimpson, .calc_inverse_simpson(assays(tse_idx)$relabundance))

    cd <- colData(tse_idx)
    expect_equal(unname(round(cd$shannon, 5)), c(2.24937, 2.76239, 2.03249))
    expect_equal(unname(round(cd$gini_simpson, 6)),
                 c(0.831372, 0.903345, 0.665749))
    expect_equal(unname(round(cd$inverse_simpson, 5)),
                 c(5.93021, 10.34606, 2.99177))
    expect_equal(unname(round(cd$coverage, 0)), c(2,3,1))
    expect_equal(unname(round(cd$fisher, 4)), c(8.8037, 10.0989, 13.2783))
    expect_equal(unname(round(cd$log_modulo_skewness, 6)), c(2.013610, 1.827198, 2.013695))
    
    # Tests that 'quantile' and 'num_of_classes' are working
    expect_equal(unname(round(colData(estimateDiversity(tse,index="log_modulo_skewness",
                                                        quantile=0.855,
                                                        num_of_classes=32)
                                      )$log_modulo_skewness, 
                              6)), c(1.814770, 1.756495, 1.842704))
    
    # Tests that .calc_skewness returns right value
    mat <- assay(tse, "counts")
    num_of_classes <- 61
    quantile <- 0.35
    
    quantile_point <- quantile(max(mat), quantile)
    cutpoints <- c(seq(0, quantile_point, length=num_of_classes), Inf)
    
    freq_table <- table(cut(mat, cutpoints), col(mat))
    test1 <- mia:::.calc_skewness(freq_table)
    
    test2 <- mia:::.calc_skewness(apply(mat, 2, function(x) {
        table(cut(x, cutpoints))
        }))
    
    expect_equal(test1, test2)
    expect_equal(round(test1, 6), c(7.256706, 6.098354, 7.278894))

    # Tests faith index with esophagus data
    for( i in c(1:(length(colnames(tse_idx)))) ){
        # Gets those taxa that are present/absent in the sample
        present <- rownames(tse)[assays(tse)$counts[,i] > 0]
        absent <- rownames(tse)[assays(tse)$counts[,i] == 0]

        # Absent taxa are dropped from the tree
        sub_tree <- ape::drop.tip(rowTree(tse_idx), absent)
        # Faith is now calculated based on the sub tree
        faith <- sum(sub_tree$edge.length)

        expect_equal(cd$faith[i], faith)
    }

    ########## Check that estimateFaith works correctly ##########
    ########## with different SE object types ##########
    
    # Creates SE from TSE by dropping, e.g., rowTree
    se <- as(tse, "SummarizedExperiment")
    
    # Add rownames because they are not included when SE is created from TSE
    rownames(se) <- rownames(tse)
    
    # Calculates "faith" TSE
    tse_only <- estimateFaith(tse)
    
    # tse_only should be TSE object
    expect_true(class(tse_only)== "TreeSummarizedExperiment")
    # tse_only should include "faith"
    expect_equal(colnames(colData(tse_only)), c(colnames(colData(tse)), "faith"))
    
    # Calculates "faith" TSE + TREE
    tse_tree <- estimateFaith(tse, tree = rowTree(tse))
    
    # tse_tree should be TSE object
    expect_true(class(tse_tree)== "TreeSummarizedExperiment")
    # tse_tree should include "faith"
    expect_equal(colnames(colData(tse_tree)), c(colnames(colData(tse)), "faith"))
    
    
    # Calculates "faith" SE + TREE
    se_tree <- estimateFaith(se, tree = rowTree(tse))
    
    # se_tree should be SE object
    expect_true(class(se_tree)== "SummarizedExperiment")
    # se_tree should include "faith"
    expect_equal(colnames(colData(se_tree)), c(colnames(colData(se)), "faith"))
    
    
    
})
