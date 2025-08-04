context("PairwiseTest")

test_that("PairwiseTest", {
    library(dplyr)
    # Data setup - reused throughout all tests
    data(GlobalPatterns, package="mia")
    tse <- GlobalPatterns
    tse <- transformAssay(tse, method = "relabundance")
    tse <- agglomerateByRank(tse, "Phylum")
    
    ## Test 1: Basic functionality with assay data
    result_tse <- addPairwiseTest(
        tse,
        assay.type = "relabundance",
        group = "SampleType",
        name = "basic_test"
    )
    
    expect_true("basic_test" %in% names(metadata(result_tse)))
    result <- metadata(result_tse)$basic_test
    expect_s3_class(result, "data.frame")
    
    # Check required columns
    expected_cols <- c("group1", "group2", "n1", "n2", "p", "p.adj", 
                       "statistic", "mean_group1", "mean_group2", "log2FC")
    expect_true(all(expected_cols %in% colnames(result)))
    
    # Check p-values are valid
    expect_true(is.numeric(result$p))
    expect_true(all(result$p >= 0 & result$p <= 1))
    expect_true(is.numeric(result$p.adj))
    expect_true(all(result$p.adj >= 0 & result$p.adj <= 1))
    
    ## Test 2: Different statistical methods
    # Wilcoxon method
    result_wilcox <- addPairwiseTest(
        tse,
        assay.type = "relabundance",
        group = "SampleType",
        significance.method = "wilcoxon",
        name = "wilcox_test"
    )
    expect_true("wilcox_test" %in% names(metadata(result_wilcox)))
    
    # t-test method
    result_ttest <- addPairwiseTest(
        tse,
        assay.type = "relabundance", 
        group = "SampleType",
        significance.method = "t.test",
        name = "t_test"
    )
    expect_true("t_test" %in% names(metadata(result_ttest)))
    
    # Kruskal-Wallis method
    result_kruskal <- addPairwiseTest(
        tse,
        assay.type = "relabundance",
        group = "SampleType", 
        significance.method = "kruskal",
        name = "kruskal_test"
    )
    expect_true("kruskal_test" %in% names(metadata(result_kruskal)))
    
    # Kruskal should return list with global and pairwise results
    kruskal_results <- metadata(result_kruskal)$kruskal_test
    expect_true("global" %in% names(attributes(kruskal_results)))
    expect_true("effect" %in% names(attributes(kruskal_results)))
    
    ## Test 3: P-value adjustment methods
    p_adjust_methods <- c("fdr", "bonferroni", "holm", "none")
    for(method in p_adjust_methods) {
        result_p_adj <- addPairwiseTest(
            tse,
            assay.type = "relabundance",
            group = "SampleType",
            p.adjust.method = method,
            name = paste0("padj_", method)
        )
        expect_true(paste0("padj_", method) %in% names(metadata(result_p_adj)))
        result_p <- metadata(result_p_adj)[[paste0("padj_", method)]]
        expect_true("p.adj" %in% colnames(result_p))
        expect_true(is.numeric(result_p$p.adj))
    }
    
    ## TODO: Test 4: rowData variables

    
    ## Test 5: colData variables
    tse <- addAlpha(tse, assay.type = "counts", index = "shannon")
    result_col <- addPairwiseTest(
        tse,
        col.var = "shannon",
        group = "SampleType", 
        name = "col_test"
    )
    expect_true("col_test" %in% names(metadata(result_col)))
    result_col_data <- metadata(result_col)$col_test
    expect_s3_class(result_col_data, "data.frame")
    
    ## Test 6: Effect size calculations
    # With effect sizes
    result_with_effect <- addPairwiseTest(
        tse,
        assay.type = "relabundance",
        group = "SampleType",
        include.effect = TRUE,
        name = "with_effect"
    )
    result_eff <- metadata(result_with_effect)$with_effect
    expect_true("effsize" %in% colnames(result_eff))
    expect_true("magnitude" %in% colnames(result_eff))
    
    # Without effect sizes
    result_no_effect <- addPairwiseTest(
        tse,
        assay.type = "relabundance",
        group = "SampleType",
        include.effect = FALSE,
        name = "no_effect"
    )
    result_no_eff <- metadata(result_no_effect)$no_effect
    expect_false("effsize" %in% colnames(result_no_eff))
    expect_false("magnitude" %in% colnames(result_no_eff))
    
    ## Test 7: Paired analysis
    # result_paired <- addPairwiseTest(
    #     tse,
    #     assay.type = "relabundance",
    #     group = "Timepoint",
    #     pair.by = "Subject",
    #     paired = TRUE,
    #     name = "paired_test"
    # )
    # expect_true("paired_test" %in% names(metadata(result_paired)))
    # result_paired_data <- metadata(result_paired)$paired_test
    # expect_s3_class(result_paired_data, "data.frame")
    
    ## Test 8: Features parameter
    selected_features <- rownames(tse)[1:2]
    result_features <- addPairwiseTest(
        tse,
        assay.type = "relabundance",
        group = "SampleType",
        features = selected_features,
        name = "features_test"
    )
    expect_true("features_test" %in% names(metadata(result_features)))
    result_feat_data <- metadata(result_features)$features_test
    expect_true(all(unique(result_feat_data$rownames) %in% selected_features))
    
    ## Test 9: Default name parameter
    result_default <- addPairwiseTest(
        tse,
        assay.type = "relabundance",
        group = "SampleType"
    )
    expect_true("pairwiseTest" %in% names(metadata(result_default)))
    
    ## Test 8: Input validation - error conditions
    # No variables specified
    expect_error(
        addPairwiseTest(tse, group = "SampleType"),
        "Please specify either 'assay.type', 'row.var', or 'col.var'"
    )
    
    # Multiple variables specified
    expect_error(
        addPairwiseTest(tse, assay.type = "counts", row.var = "Kingdom", 
                        group = "SampleType"),
        "Please specify either 'assay.type', 'row.var', or 'col.var'"
    )
    
    # Invalid grouping variable
    expect_error(
        addPairwiseTest(tse, assay.type = "counts", group = "nonexistent"),
        "must be.*character value from the following options"
    )
    
    # Invalid statistical method
    expect_error(
        addPairwiseTest(tse, assay.type = "counts", group = "SampleType",
                        significance.method = "invalid"),
        "'significance.method' must be one of"
    )
    
    # Invalid p.adjust.method
    expect_error(
        addPairwiseTest(tse, assay.type = "counts", group = "SampleType",
                        p.adjust.method = 123),
        "'p.adjust.method' must be a character string"
    )
})
