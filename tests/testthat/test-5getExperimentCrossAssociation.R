
context("getExperimentCrossAssociation")

test_that("getExperimentCrossAssociation", {
    
    # Get data
    mae <- microbiomeDataSets::peerj32()
    ############################### Test input ###############################
    expect_error(getExperimentCrossAssociation(mae,
                                                experiment1 = 3,
                                                experiment2 = 2,
                                                abund_values1 = "counts",
                                                abund_values2 = "counts",
                                                method = "spearman",
                                                mode = "table",
                                                p_adj_method = "fdr",
                                                p_adj_threshold = 0.05,
                                                cor_threshold = NULL,
                                                sort = FALSE,
                                                filter_self_correlations = FALSE,
                                                verbose = TRUE))
     expect_error(getExperimentCrossAssociation(mae,
                                                experiment1 = 1,
                                                experiment2 = TRUE,
                                                abund_values1 = "counts",
                                                abund_values2 = "counts",
                                                method = "spearman",
                                                mode = "table",
                                                p_adj_method = "fdr",
                                                p_adj_threshold = 0.05,
                                                cor_threshold = NULL,
                                                sort = FALSE,
                                                filter_self_correlations = FALSE,
                                                verbose = TRUE))
     expect_error(getExperimentCrossAssociation(mae,
                                                experiment1 = 1,
                                                experiment2 = 2,
                                                abund_values1 = "test",
                                                abund_values2 = "counts",
                                                method = "spearman",
                                                mode = "table",
                                                p_adj_method = "fdr",
                                                p_adj_threshold = 0.05,
                                                cor_threshold = NULL,
                                                sort = FALSE,
                                                filter_self_correlations = FALSE,
                                                verbose = TRUE))
     expect_error(getExperimentCrossAssociation(mae,
                                                experiment1 = 1,
                                                experiment2 = 2,
                                                abund_values1 = "counts",
                                                abund_values2 = "counts",
                                                method = 1,
                                                mode = "table",
                                                p_adj_method = "fdr",
                                                p_adj_threshold = 0.05,
                                                cor_threshold = NULL,
                                                sort = FALSE,
                                                filter_self_correlations = FALSE,
                                                verbose = TRUE))
     expect_error(getExperimentCrossAssociation(mae,
                                                experiment1 = 1,
                                                experiment2 = 2,
                                                abund_values1 = "counts",
                                                abund_values2 = "counts",
                                                method = FALSE,
                                                mode = "table",
                                                p_adj_method = "fdr",
                                                p_adj_threshold = 0.05,
                                                cor_threshold = NULL,
                                                sort = FALSE,
                                                filter_self_correlations = FALSE,
                                                verbose = TRUE))
     expect_error(getExperimentCrossAssociation(mae,
                                                experiment1 = 1,
                                                experiment2 = 2,
                                                abund_values1 = "counts",
                                                abund_values2 = "counts",
                                                method = "pearson",
                                                mode = TRUE,
                                                p_adj_method = "fdr",
                                                p_adj_threshold = 0.05,
                                                cor_threshold = NULL,
                                                sort = FALSE,
                                                filter_self_correlations = FALSE,
                                                verbose = TRUE))
     expect_error(getExperimentCrossAssociation(mae,
                                                experiment1 = 1,
                                                experiment2 = 2,
                                                abund_values1 = "counts",
                                                abund_values2 = "counts",
                                                method = "pearson",
                                                mode = "matrix",
                                                p_adj_method = 1,
                                                p_adj_threshold = 0.05,
                                                cor_threshold = NULL,
                                                sort = FALSE,
                                                filter_self_correlations = FALSE,
                                                verbose = TRUE))
     expect_error(getExperimentCrossAssociation(mae,
                                                experiment1 = 1,
                                                experiment2 = 2,
                                                abund_values1 = "counts",
                                                abund_values2 = "counts",
                                                method = "pearson",
                                                mode = "matrix",
                                                p_adj_method = "fdr",
                                                p_adj_threshold = 2,
                                                cor_threshold = NULL,
                                                sort = FALSE,
                                                filter_self_correlations = FALSE,
                                                verbose = TRUE))
     expect_error(getExperimentCrossAssociation(mae,
                                                experiment1 = 1,
                                                experiment2 = 2,
                                                abund_values1 = "counts",
                                                abund_values2 = "counts",
                                                method = "pearson",
                                                mode = "matrix",
                                                p_adj_method = "fdr",
                                                p_adj_threshold = TRUE,
                                                cor_threshold = NULL,
                                                sort = FALSE,
                                                filter_self_correlations = FALSE,
                                                verbose = TRUE))
     expect_error(getExperimentCrossAssociation(mae,
                                                experiment1 = 1,
                                                experiment2 = 2,
                                                abund_values1 = "counts",
                                                abund_values2 = "counts",
                                                method = "pearson",
                                                mode = "matrix",
                                                p_adj_method = "fdr",
                                                p_adj_threshold = 0.1,
                                                cor_threshold = 2,
                                                sort = FALSE,
                                                filter_self_correlations = FALSE,
                                                verbose = TRUE))
     expect_error(getExperimentCrossAssociation(mae,
                                                experiment1 = 1,
                                                experiment2 = 2,
                                                abund_values1 = "counts",
                                                abund_values2 = "counts",
                                                method = "pearson",
                                                mode = "matrix",
                                                p_adj_method = "fdr",
                                                p_adj_threshold = 0.1,
                                                cor_threshold = TRUE,
                                                sort = FALSE,
                                                filter_self_correlations = FALSE,
                                                verbose = TRUE))
     expect_error(getExperimentCrossAssociation(mae,
                                                experiment1 = 1,
                                                experiment2 = 2,
                                                abund_values1 = "counts",
                                                abund_values2 = "counts",
                                                method = "pearson",
                                                mode = "matrix",
                                                p_adj_method = "fdr",
                                                p_adj_threshold = 0.1,
                                                cor_threshold = NULL,
                                                sort = 1,
                                                filter_self_correlations = FALSE,
                                                verbose = TRUE))
     expect_error(getExperimentCrossAssociation(mae,
                                                experiment1 = 1,
                                                experiment2 = 2,
                                                abund_values1 = "counts",
                                                abund_values2 = "counts",
                                                method = "pearson",
                                                mode = "matrix",
                                                p_adj_method = "fdr",
                                                p_adj_threshold = 0.1,
                                                cor_threshold = NULL,
                                                sort = TRUE,
                                                filter_self_correlations = 1,
                                                verbose = TRUE))
     expect_error(getExperimentCrossAssociation(mae,
                                                experiment1 = 1,
                                                experiment2 = 2,
                                                abund_values1 = "counts",
                                                abund_values2 = "counts",
                                                method = "pearson",
                                                mode = "matrix",
                                                p_adj_method = "fdr",
                                                p_adj_threshold = 0.1,
                                                cor_threshold = NULL,
                                                sort = TRUE,
                                                filter_self_correlations = TRUE,
                                                verbose = 1))
     expect_error(getExperimentCrossAssociation(mae[[1]],
                                                assay(mae1[[2]]),
                                                experiment1 = 1,
                                                experiment2 = 2,
                                                abund_values1 = "counts",
                                                abund_values2 = "counts",
                                                method = "pearson",
                                                mode = "matrix",
                                                p_adj_method = "fdr",
                                                p_adj_threshold = 0.1,
                                                cor_threshold = NULL,
                                                sort = TRUE,
                                                filter_self_correlations = TRUE,
                                                verbose = 1))
     expect_error(getExperimentCrossAssociation(mae[[1]],
                                                NULL,
                                                experiment1 = 1,
                                                experiment2 = 2,
                                                abund_values1 = "counts",
                                                abund_values2 = "counts",
                                                method = "pearson",
                                                mode = "matrix",
                                                p_adj_method = "fdr",
                                                p_adj_threshold = 0.1,
                                                cor_threshold = NULL,
                                                sort = TRUE,
                                                filter_self_correlations = TRUE,
                                                verbose = 1))
     ############################# Test input end #############################
     # Test that association is calculated correctly with numeric data
     # Result from
     # d1 <- t(assay(mae[[1]]))
     # d2 <- t(assay(mae[[2]]))
     # cc <- microbiome::associate(d1, d2, method='pearson', mode= "table")
    cor_compare <- c(-0.67473156, 0.19980484, -0.19942102, -0.17346641, -0.16954081, 
                     -0.15477367, -0.15279005, 0.08792788, 0.08443186)
    p_adj_compare <- c(0.001247967, 0.784472862, 0.785332288, 0.830362548, 
                       0.836148425, 0.856762552, 0.859203260, 0.938444366, 0.942610008)
    # Calculate correlation
    cor <- testExperimentCrossAssociation(mae, method = "pearson", 
                                          p_adj_threshold = NULL, show_warnings = FALSE)
    # Take only specific taxa and lipids
    df <- cor[cor$Var1 %in% c("Fusobacteria", "Campylobacter", "Actinomycetaceae") & 
                 cor$Var2 %in% c("PE(48:7)", "TG(50:0)", "SM(d18:1/18:0)"), ]
    # Sort the data, so that lowest p-values are first
    df <- df[order(df$p_adj), ]
    # Correlation values and p-values should be the same
    expect_equal(round(df$cor, 7), round(cor_compare, 7))
    expect_equal(round(df$p_adj, 7), round(p_adj_compare, 7))
    
    # Test that association is calculated correctly with factor data
    # Create a dummy data
    assay1 <- matrix(rep(c("A", "B", "B"), 20*30), 
                     nrow = 20, ncol = 30)
    assay2 <- matrix(rep(c("A", "B", "A", "A", "B", "B"), 20*30), 
                     nrow = 20, ncol = 30)
    # Reference
    # ref <- c()
    # for(i in 1:20){
    #     ref <- c(ref, GoodmanKruskal::GKtau(assay1[i, ], assay2[i, ])$tauxy)
    # }
    ref <- c(0.25, 1.00, 0.25, 1.00, 0.25, 1.00, 0.25, 1.00, 0.25, 1.00, 0.25, 
             1.00, 0.25, 1.00, 0.25, 1.00, 0.25, 1.00, 0.25, 1.00)
    # Calculate values for 20 feature-pairs
    result <- c()
    for(i in 1:20){
        result <- c(result, .calculate_gktau(assay1[i, ], assay2[i, ])$estimate)
    }
    # Values should be the same
    expect_equal(round(result, 4), round(ref, 4))
    
    mae_sub <- mae[1:10, 1:10]
    # Test that output is in correct type
    expect_true( is.data.frame(
        testExperimentCrossAssociation(mae_sub, p_adj_threshold = NULL, show_warnings = FALSE)) )
    expect_true( is.data.frame(
        getExperimentCrossAssociation(mae_sub, test_significance = TRUE, 
                                      p_adj_threshold = NULL, show_warnings = FALSE)) )
    expect_true( is.data.frame(getExperimentCrossAssociation(mae_sub, show_warnings = FALSE)) )
    # There should not be any p-values that are under 0
    expect_true( is.null(
        testExperimentCrossAssociation(mae_sub, p_adj_threshold = 0, show_warnings = FALSE)) )
    # Test that output is in correct type
    expect_true( is.list(
        testExperimentCrossAssociation(mae_sub, mode = "matrix", 
                                          p_adj_threshold = NULL, show_warnings = FALSE)) )
    expect_true( is.list(
        getExperimentCrossAssociation(mae_sub, test_significance = TRUE, 
                                      mode = "matrix", 
                                      p_adj_threshold = NULL, show_warnings = FALSE)) )
    expect_true( is.matrix(getExperimentCrossAssociation(mae_sub, mode = "matrix", 
                                                         show_warnings = FALSE)) )
    
    # There should not be any p-values that are under 0
    expect_true( is.null(
        testExperimentCrossAssociation(mae_sub, p_adj_threshold = 0, mode = "matrix", show_warnings = FALSE)) )
    
    # When correlation between same assay is calculated, calculation is made faster
    # by not calculating duplicates
    expect_error(testExperimentCrossAssociation(mae, experiment1 = 1, experiment2 = 1, 
                                                show_warnings = FALSE, symmetric = "TRUE"))
    expect_error(testExperimentCrossAssociation(mae, experiment1 = 1, experiment2 = 1, 
                                                show_warnings = FALSE, symmetric = 1))
    expect_error(testExperimentCrossAssociation(mae, experiment1 = 1, experiment2 = 1, 
                                                show_warnings = FALSE, symmetric = NULL))
    expect_error(testExperimentCrossAssociation(mae, experiment1 = 1, experiment2 = 1, 
                                                show_warnings = FALSE, symmetric = c(TRUE, TRUE)))
    
    time <- system.time(
        cor <-  testExperimentCrossAssociation(mae, experiment1 = 1, experiment2 = 1, 
                                           show_warnings = FALSE, symmetric = TRUE)
    )
    time2 <- system.time(
        cor2 <-  testExperimentCrossAssociation(mae, experiment1 = 1, experiment2 = 1, 
                                                show_warnings = FALSE)
    )
    # Get random variables and test that their duplicates are equal
    for(i in 1:10 ){
        random_var1 <- sample(cor$Var1, 1)
        random_var2 <- sample(cor$Var1, 1)
        expect_equal(as.numeric(cor[cor$Var1 == random_var1 & cor$Var2 == random_var2, c("cor", "pval", "p_adj")]),
                     as.numeric(cor[cor$Var1 == random_var2 & cor$Var2 == random_var1, c("cor", "pval", "p_adj")]))
    }
    expect_equal(cor, cor2)
    # Test that symmetric = TRUE was faster
    expect_true(time[3] < time2[3])
    # Test that paired samples work correctly
    tse1 <- mae[[1]]
    tse2 <- mae[[1]]
    # Convert assay to have random values
    mat <- matrix(sample(0:100, nrow(tse2)*ncol(tse2), replace = TRUE), 
                  nrow = nrow(tse2), ncol = ncol(tse2))
    colnames(mat) <- colnames(tse2)
    rownames(mat) <- rownames(tse2)
    assay(tse2) <- mat
    # Calculate with paired samples
    cor_paired <- testExperimentCrossAssociation(tse1,
                                                 experiment2 = tse2,  
                                                 paired = TRUE, 
                                                 direction = "col", 
                                                 show_warnings = FALSE)
    # Calculate all pairs
    cor <- testExperimentCrossAssociation(tse1,
                                          experiment2 = tse2,  
                                          direction = "col", 
                                          show_warnings = FALSE)
    # Take only pairs that are paired
    cor <- cor[cor$Var1 == cor$Var2, ]
    rownames(cor) <- NULL
    
    # Should be equal
    expect_equal(cor[, c("cor", "pval")], 
                 cor_paired[, c("cor", "pval")])
    
    # Test that result does not depend on names (if there are equal names)
    tse <- mae[[1]]
    rownames(tse)[1:10] <- rep("Unknown", 10)
    cor_table <- testExperimentCrossCorrelation(tse, show_warnings = TRUE)
    cor_table_ref <- testExperimentCrossCorrelation(mae[[1]], show_warnings = TRUE)
    expect_equal(cor_table[ , 3:5], cor_table_ref[ , 3:5])
})
