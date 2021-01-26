context("transformCounts")

# historic left overs
test_that("relabundance", {
    mat <- matrix(1:60, nrow = 6)
    df <- DataFrame(n = c(1:6))

    se <- SummarizedExperiment(assays = list(counts = mat),
                               rowData = df)
    actual <- relAbundanceCounts(se)

    expect_equal(assay(actual,"relabundance"),
                 relabundance(actual))
    rel_mat <- relabundance(actual)
    f <- rev(seq_len(ncol(relabundance(actual))))
    relabundance(actual) <- relabundance(actual)[,f]
    expect_equal(rel_mat,
                 relabundance(actual)[,f])
})


test_that("transformCounts", {

    # TSE object
    data(GlobalPatterns)
    tse <- GlobalPatterns

    # No method specified. Should be an error.
    expect_error(mia::transformCounts(tse))

    # Method is not provided. Should be an error.
    expect_error(mia::transformCounts(tse, method="test"))

    # Name is a vector of 2. Should be an error.
    expect_error(mia::transformCounts(tse, method="relabundance", name = c("123", "456")))

    # Pseudocount is a string. Should be an error.
    expect_error(mia::transformCounts(tse, method="relabundance", pseudocount = "pseudocount"))

    # Counts table should not be changed
    expect_equal(assays(mia::transformCounts(tse, method = "pa"))$counts, assays(tse)$counts)

    ################################################################

    # Calculates relative abundances. Should be equal.
    expect_equal(assays(mia::transformCounts(tse, method = "relabundance"))$relabundance,
                 apply(assays(tse)$counts, 2, FUN=function(x){
                     x/sum(x)
                 }))

    expect_equal(assays(mia::relAbundanceCounts(tse, pseudocount = 12))$relabundance,
                 apply((assays(tse)$counts + 12), 2, FUN=function(x){
                     x/sum(x)
                 }))

    mat <- matrix(1:60, nrow = 6)
    df <- DataFrame(n = c(1:6))
    expect_error(relAbundanceCounts(SummarizedExperiment(assays = list(mat = mat),
                                                         rowData = df)),
                 "'abund_values' must be a valid name of assays")

    se <- SummarizedExperiment(assays = list(counts = mat),
                               rowData = df)
    expect_error(relAbundanceCounts(se, name = FALSE),
                 "'name' must be a non-empty single character value")

    actual <- relAbundanceCounts(se)
    expect_named(assays(actual), c("counts", "relabundance"))

    expect_equal(assay(actual,"relabundance")[,1],
                 seq.int(1,6)/21)



    ##############################################################
    # Calculates log10 transformation with pseudocount = 1. Should be equal.
    expect_equal(assays(mia::transformCounts(tse, method = "log10", pseudocount = 1))$log10,
                 apply(assays(tse)$counts, 2, FUN=function(x){
                     log10(x+1)
                 }))

    # Calculates pa transformation. Should be equal.
    expect_equal(as.vector(assays(mia::transformCounts(tse, method = "pa"))$pa),
                 as.integer(assays(tse)$counts > 0))

    # TSE object
    data(esophagus)
    tse <- esophagus

    ###############################################################
    # Calculates Hellinger transformation. Should be equal.
    # Calculates relative abundance table
    relative <- assays(mia::transformCounts(tse, method = "relabundance", name = "relative"))$relative

    expect_equal(assays(mia::transformCounts(tse, method = "hellinger", name = "test_123"))$test_123,
                 apply(relative, 2, FUN=function(x){
                     sqrt(x)
                 }))

    ##############################################################
    # Calculates clr-transformation. Should be equal.

    # Random pseudocount
    pseudonumber <- runif(1, 1, 100)

    # Calculates relative abundance table
    relative <- assays(mia::transformCounts(tse, method = "relabundance", pseudocount = pseudonumber))$relabundance

    expect_equal(
        assays(mia::transformCounts(tse, method = "clr", pseudocount = pseudonumber))$clr,
        apply(relative, 2, FUN=function(x){
            log(x) - mean(log(x))
        }))

    #############################################################
    # Tests that samples have correct names
    expect_equal(colnames(assays(mia::transformCounts(tse, method = "clr", pseudocount = 1))$clr),
                 colnames(assays(tse)$counts))

    # Tests that otus have correct names
    expect_equal(rownames(assays(mia::transformCounts(tse, method = "hellinger", pseudocount = 1000))$hellinger),
                 rownames(assays(tse)$counts))

    #############################################################
    # Calculates Z-transformation for features
    # Information collected with microbiome package
    B <- c(0.9557828, -0.5773503, -0.1110960, -1.0750696, -0.5773503, -0.5773503, -0.3332626,
           -0.5773503,  0.9626491, -0.5773503, -0.5773503, -0.5773503, -0.5773503, -0.5773503,
           -0.5773503,  1.1547005, -0.5773503, -0.5773503,  1.1547005, -1.0750696)

    C <- c(0.08323189,  1.15470054,  1.05090886,  0.90245812, -0.57735027, -0.57735027, -0.79081431,
           1.15470054, -1.03357456, -0.57735027, -0.57735027, -0.57735027, -0.57735027, -0.57735027,
           1.15470054, -0.57735027,  1.15470054, -0.57735027, -0.57735027,  0.90245812)

    D <- c(-1.0390147, -0.5773503, -0.9398129,  0.1726115,  1.1547005,  1.1547005,  1.1240769, -0.5773503,
           0.0709255,  1.1547005,  1.1547005,  1.1547005,  1.1547005,  1.1547005, -0.5773503, -0.5773503,
           -0.5773503,  1.1547005, -0.5773503,  0.1726115)

    df <- data.frame(B, C, D)
    df <- round(df, 7)
    rownames(df) <- rownames(assays(tse)$counts[1:20,])

    # Rounded, because hard-coded values have only 7 decimals
    expect_equal(as.data.frame(round(assays(mia::ZTransform(tse, pseudocount = 1))$ZTransform,7))[1:20,],
                 df)

})

