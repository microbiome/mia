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

    testTransformations <- function(tse){
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
        expect_equal(as.matrix(assays(mia::transformCounts(tse, method = "relabundance"))$relabundance),
                     apply(as.matrix(assay(tse,"counts")), 2, FUN=function(x){
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
        expect_equal(as.matrix(assays(mia::transformCounts(tse, method = "log10", pseudocount = 1))$log10),
                     apply(as.matrix(assay(tse, "counts")), 2, FUN=function(x){
                         log10(x+1)
                     }))

        ###############################################################

        # Calculates pa transformation. Should be equal.
        actual <- assay(mia::transformCounts(tse, method = "pa"),"pa")
        expect_equal(as.vector(actual),
                     as.integer(as.matrix(assay(tse, "counts")) > 0))
        expect_equal(type(actual),"integer")

        # Calculates pa transformation. Should be equal.
        actual <- assay(mia::transformCounts(tse, method = "pa", threshold = 12.5),"pa")
        expect_equal(as.vector(actual),
                     as.integer(as.matrix(assay(tse, "counts")) > 12.5))
        expect_equal(type(actual),"integer")

        ###############################################################
        # Calculates Hellinger transformation. Should be equal.
        # Calculates relative abundance table
        relative <- assays(mia::transformCounts(tse, method = "relabundance", name = "relative"))$relative

        expect_equal(as.matrix(assays(mia::transformCounts(tse, method = "hellinger", name = "test_123"))$test_123),
                     apply(as.matrix(relative), 2, FUN=function(x){
                         sqrt(x)
                     }))

        ##############################################################
        # Calculates clr-transformation. Should be equal.

        # Random pseudocount
        pseudonumber <- runif(1, 1, 100)

        # Calculates relative abundance table
        relative <- assays(mia::transformCounts(tse, method = "relabundance", pseudocount = pseudonumber))$relabundance

        expect_equal(
            as.matrix(assays(mia::transformCounts(tse, method = "clr", pseudocount = pseudonumber))$clr),
            apply(as.matrix(relative) + pseudonumber, 2, FUN=function(x){
                log(x) - mean(log(x))
            }))
        
        
        tse2 <- relAbundanceCounts(tse)
        expect_true(all(round(as.matrix(assays(mia::transformCounts(tse2, method = "clr", 
                                                                    abund_values = "relabundance", pseudocount = pseudonumber))$clr) -
                                  as.matrix(assays(mia::transformCounts(tse, method = "clr", pseudocount = pseudonumber))$clr),
                              8) == 0))

        #############################################################
        # Tests that samples have correct names
        expect_equal(colnames(assays(transformCounts(tse, method = "clr", pseudocount = 1))$clr),
                     colnames(assays(tse)$counts))

        # Tests that otus have correct names
        expect_equal(rownames(assays(transformCounts(tse, method = "hellinger", pseudocount = 1000))$hellinger),
                     rownames(assays(tse)$counts))
        
        ############################################################
        # Calculates rank
        tse_rank <- transformCounts(tse, method = "rank")
        # Expect that assay contains count and rank table
        expect_equal(names(assays(tse_rank)), c("counts", "rank") )
        
        for( i in c(1:10) ){
            # Gets columns from 'rank' table
            ranks <- assay(tse_rank, "rank")[,i]
            # Gets columns from 'counts' table, and calculates ranks
            counts_compare <- assay(tse_rank, "counts")[,i]
            ranks_compare <- rank(counts_compare, na.last = "keep", ties.method = "min")
            # Expect that they are equal
            expect_equal(ranks, ranks_compare)
        }
        
        # Calculates rank with pseudocount
        tse_rank_pseudo <- transformCounts(tse, method = "rank", pseudocount = runif(1, 0, 1000))
        # Pseudocount should not change the rank
        expect_equal(tse_rank, tse_rank_pseudo)
        #############################################################
        # SE object
        data("esophagus")
        se <- esophagus
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
        rownames(df) <- rownames(assays(se)$counts[1:20,])

        # Rounded, because hard-coded values have only 7 decimals
        expect_equal(as.data.frame(round(assays(mia::ZTransform(se, pseudocount = 1))$ZTransform,7))[1:20,],
                     df)
    }

    # TSE object
    data(GlobalPatterns)
    tse <- GlobalPatterns
    testTransformations(GlobalPatterns)
    assay(tse,"counts") <- DelayedArray(assay(tse,"counts"))
    testTransformations(tse)
})

