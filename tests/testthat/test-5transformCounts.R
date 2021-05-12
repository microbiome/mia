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
        # Calculates log10 transformation with pseudocount. Should be equal.
	tmp <- mia::transformCounts(tse, method = "log10", pseudocount = 1)	
        ass <- assays(tmp)$log10
        expect_equal(as.matrix(ass),
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
            apply(as.matrix(relative), 2, FUN=function(x){
                log(x) - mean(log(x))
            }))

        tse <- relAbundanceCounts(tse)
	
        #############################################################
        # Tests that samples have correct names
        expect_equal(colnames(assays(transformCounts(tse, method = "clr", pseudocount = 1))$clr),
                     colnames(assays(tse)$counts))

        # Tests that otus have correct names
        expect_equal(rownames(assays(transformCounts(tse, method = "hellinger", pseudocount = 1000))$hellinger),
                     rownames(assays(tse)$counts))
        
        ############################################################
        # Calculates rank
        se_rank <- transformCounts(se, method = "rank")
        # Expect that assay contains count and rank table
        expect_equal(names(assays(se_rank)), c("counts", "rank") )
        
        for( i in c(1:10) ){
            # Gets columns from 'rank' table
            ranks <- assay(se_rank, "rank")[,i]
            # Gets columns from 'counts' table, and calculates ranks
            counts_compare <- assay(se_rank, "counts")[,i]
            ranks_compare <- rank(counts_compare, na.last = "keep", ties.method = "first")
            # Expect that they are equal
            expect_equal(ranks, ranks_compare)
        }
        
        # Calculates rank with pseudocount
        se_rank_pseudo <- transformCounts(se, method = "rank", pseudocount = runif(1, 0, 1000))
        # Pseudocount should not change the rank
        expect_equal(se_rank, se_rank_pseudo)

        # For other options for rank calculations, call the colRanks directly:
	# data("esophagus"); x <- esophagus;  
        # assay(x, "rank_average") <- t(colRanks(assay(x, "counts"), ties.method="average"))

        #############################################################
        # SE object
        data("esophagus")
        se <- esophagus
	
        # Calculates Z-transformation for features	
        xx <- t(scale(t(assay(se, "counts") + 1)))
	expect_equal(max(abs(assays(mia::ZTransform(se, pseudocount = 1))$ZTransform - xx), na.rm=TRUE), 0, tolerance = 1e-16)

    }

    # TSE object
    data(GlobalPatterns)
    tse <- GlobalPatterns
    testTransformations(GlobalPatterns)
    assay(tse,"counts") <- DelayedArray(assay(tse,"counts"))
    testTransformations(tse)
})

