context("getDA")

test_that("getDA and add*DA methods work correctly", {
    library(dplyr)
    data(GlobalPatterns, package = "mia")
    tse <- GlobalPatterns
    tse <- transformAssay(tse, method = "relabundance")
    tse <- tse[690:700, ]
    
    ## Test 1: addPairwiseDA basic functionality
    result_tse <- addPairwiseDA(
        tse,
        assay.type = "relabundance",
        group = "SampleType",
        name = "basic_test"
    )
    expect_true("basic_test" %in% names(metadata(result_tse)))
    result <- metadata(result_tse)$basic_test
    expect_s3_class(result, "data.frame")
    
    # Check required columns (statistic and magnitude should be removed)
    expected_cols <- c(
        "group1",       "group2",       "p",            "p.adj",       
        "effsize",      "log2FC",       "mean_group1",  "mean_group2",  
        "n1",           "n2",           ".y.",          "p.adj.signif"
    )
    expect_true(all(expected_cols %in% colnames(result)))
    
    # p-values valid
    expect_true(is.numeric(result$p))
    expect_true(all(result$p >= 0 & result$p <= 1))
    expect_true(is.numeric(result$p.adj))
    expect_true(all(result$p.adj >= 0 & result$p.adj <= 1))
    
    ## Test 2: da.method variation
    for (method in c("wilcoxon", "ttest")) {
        name <- paste0(method, "_test")
        res_obj <- addPairwiseDA(
            tse,
            assay.type = "relabundance",
            group = "SampleType",
            da.method = method,
            name = name
        )
        expect_true(name %in% names(metadata(res_obj)))
    }
    
    ## Test 3: Posthoc (dunns) and omnibus attributes
    # Dunn's posthoc after kruskal
    result_kruskal <- addPosthocDA(
        tse,
        assay.type = "relabundance",
        group = "SampleType",
        da.method = "dunns",
        name = "kruskal_posthoc"
    )
    expect_true("kruskal_posthoc" %in% names(metadata(result_kruskal)))
    kruskal_res <- metadata(result_kruskal)$kruskal_posthoc
    expect_true(is.data.frame(kruskal_res))
    expect_true("global_test" %in% names(attributes(kruskal_res)))
    
    ## Test 4: p.adjust.method options
    for (padj in c("fdr", "bonferroni", "holm", "none")) {
        nm <- paste0("padj_", padj)
        res_p <- addPairwiseDA(
            tse,
            assay.type = "relabundance",
            group = "SampleType",
            p.adjust.method = padj,
            name = nm
        )
        expect_true(nm %in% names(metadata(res_p)))
        df_p <- metadata(res_p)[[nm]]
        expect_true("p.adj" %in% colnames(df_p))
        expect_true(is.numeric(df_p$p.adj))
    }
    
    ## Test 5: col.var option
    tse2 <- addAlpha(GlobalPatterns, assay.type = "counts", index = "shannon")
    result_col <- addPairwiseDA(
        tse2,
        col.var = "shannon",
        group = "SampleType",
        name = "col_test"
    )
    expect_true("col_test" %in% names(metadata(result_col)))
    df_col <- metadata(result_col)$col_test
    expect_s3_class(df_col, "data.frame")
    
    ## Test 6: include.effect TRUE/FALSE
    res_eff <- addPairwiseDA(
        tse,
        assay.type = "relabundance",
        group = "SampleType",
        include.effect = TRUE,
        name = "with_effect"
    )
    df_eff <- metadata(res_eff)$with_effect
    expect_true("effsize" %in% colnames(df_eff))
    
    res_noeff <- addPairwiseDA(
        tse,
        assay.type = "relabundance",
        group = "SampleType",
        include.effect = FALSE,
        name = "no_effect"
    )
    df_noeff <- metadata(res_noeff)$no_effect
    expect_false("effsize" %in% colnames(df_noeff))
    
    ## Test 7: features filtering
    feats <- rownames(tse)[1:2]
    res_feat <- addPairwiseDA(
        tse,
        assay.type = "relabundance",
        group = "SampleType",
        features = feats,
        name = "feat_test"
    )
    df_feat <- metadata(res_feat)$feat_test
    expect_true(all(df_feat$rownames %in% feats))
    
    ## Test 8: default name
    res_def <- addPairwiseDA(
        tse,
        assay.type = "relabundance",
        group = "SampleType"
    )
    # default name is "pairwiseDA"
    expect_true("pairwiseDA" %in% names(metadata(res_def)))
    
    ## Test 9: input validation errors
    expect_error(
        addPairwiseDA(tse, group = "SampleType"),
        "Please specify either 'assay.type', 'row.var', or 'col.var'"
    )
    expect_error(
        addPairwiseDA(tse, assay.type = "counts", row.var = "Kingdom",
                      group = "SampleType"),
        "Please specify either 'assay.type', 'row.var', or 'col.var'"
    )
    expect_error(
        addPairwiseDA(tse, assay.type = "counts", group = "SampleType",
                      da.method = "invalid"),
        "Unsupported method"
    )
    expect_error(
        addPairwiseDA(tse, assay.type = "counts", group = "SampleType",
                      p.adjust.method = 123),
        "'p.adjust.method' must be a character string"
    )
})
