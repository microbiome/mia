context("transformCounts")

test_that("transformCounts", {

    testTransformations <- function(tse){
        # No method specified. Should be an error.
        expect_error(mia::transformCounts(tse))
        # Method is not provided. Should be an error.
        expect_error(mia::transformCounts(tse, method="test"))

        # Name is a vector of 2. Should be an error.
        expect_error(mia::transformCounts(tse, method="relabundance", name = c("123", "456")))

        # Pseudocount is a string. Should be an error.
        expect_error(mia::transformCounts(tse, method="log10", pseudocount = "pseudocount"))
        expect_error(mia::transformCounts(tse, method="log2", pseudocount = FALSE))

        # Counts table should not be changed
        expect_equal(assays(mia::transformCounts(tse, method = "pa"))$counts, assays(tse)$counts,
                     check.attributes = FALSE)

        ############################# RELATIVE ABUNDANCE #######################
        # Calculates relative abundances. Should be equal.
        expect_equal(as.matrix(assays(mia::transformCounts(tse, method = "relabundance"))$relabundance),
                     apply(as.matrix(assay(tse,"counts")), 2, FUN=function(x){
                         x/sum(x)
                     }), check.attributes = FALSE)
        
        mat <- matrix(1:60, nrow = 6)
        df <- DataFrame(n = c(1:6))
        expect_error(relAbundanceCounts(SummarizedExperiment(assays = list(mat = mat),
                                                             rowData = df)),
                     "'assay.type' must be a valid name of assays")

        se <- SummarizedExperiment(assays = list(counts = mat),
                                   rowData = df)
        expect_error(relAbundanceCounts(se, name = FALSE),
                     "'name' must be a non-empty single character value")

        actual <- relAbundanceCounts(se)
        expect_named(assays(actual), c("counts", "relabundance"))

        expect_equal(assay(actual,"relabundance")[,1],
                     seq.int(1,6)/21, check.attributes = FALSE)

        ########################### LOG10 ######################################
        # Calculates log10 transformation with pseudocount. Should be equal.
	    tmp <- mia::transformCounts(tse, method = "log10", pseudocount = 1)	
        ass <- assays(tmp)$log10
        expect_equal(as.matrix(ass),
                     apply(as.matrix(assay(tse, "counts")), 2, FUN=function(x){
                         log10(x+1)
                     }), check.attributes = FALSE)
        # Tests transformCounts(MARGIN = "features"), , calculates log10 transformation with pseudocount.
        # Should be equal.
        tmp <- mia::transformCounts(tse, MARGIN = "features", method = "log10", pseudocount = 1)
        ass <- assays(tmp)$log10
        expect_equal(as.matrix(ass),
                     t(apply(as.matrix(t(assay(tse, "counts"))), 2, FUN=function(x){
                         log10(x+1)
                     })), check.attributes = FALSE)
        
        ########################### LOG2 ######################################
        # Calculates log2 transformation with pseudocount. Should be equal.
        tmp <- mia::transformCounts(tse, method = "log2", pseudocount = 5)	
        ass <- assays(tmp)$log2
        expect_equal(as.matrix(ass),
                     apply(as.matrix(assay(tse, "counts")), 2, FUN=function(x){
                         log2(x+5)
                     }), check.attributes = FALSE)
        
        ########################## PA ##########################################
        # Calculates pa transformation. Should be equal.
        actual <- assay(mia::transformCounts(tse, method = "pa"),"pa")
        expect_equal(as.vector(actual),
                     as.integer(as.matrix(assay(tse, "counts")) > 0),
                     check.attributes = FALSE)
        expect_equal(type(actual),"double")
        expect_true(all(actual == 1 | actual == 0))
        
        # Tests transformCounts(MARGIN = "features"), , calculates pa transformation. Should be equal.
        actual <- assay(mia::transformCounts(tse, MARGIN = "features", method = "pa"),"pa")
        expect_equal(as.vector(actual),
                     as.integer(t(as.matrix(t(assay(tse, "counts"))) > 0)))
        expect_equal(type(actual),"double")
        expect_true(all(actual == 1 | actual == 0))
        
        ######################## HELLINGER #####################################
        # Calculates Hellinger transformation. Should be equal.
        # Calculates relative abundance table
        relative <- assays(mia::transformCounts(tse, method = "relabundance", name = "relative"))$relative

        expect_equal(as.matrix(assays(mia::transformCounts(tse, method = "hellinger", name = "test_123"))$test_123),
                     apply(as.matrix(relative), 2, FUN=function(x){
                         sqrt(x)
                     }), check.attributes = FALSE)
        
        ############################### CLR ####################################
        # Calculates clr-transformation. Should be equal.

        # Random pseudocount
        pseudonumber <- runif(1, 1, 100)
        
        tse <- transformCounts(tse, method = "relabundance")
        # Calculates relative abundance table
        relative <- assay(tse, "relabundance")
        relative <- relative + pseudonumber
        # Tests clr
        # Calc CLRs
        mat <- assays(mia::transformCounts(tse, assay.type = "relabundance",
                                           method = "clr", pseudocount = pseudonumber))$clr
        mat_comp <- apply(as.matrix(relative), 2, FUN=function(x){
            log(x) - mean(log(x))
        })
        # Remove atributes since vegan adds additional ones
        attributes(mat) <- NULL
        attributes(mat_comp) <- NULL
        # Compare
        expect_equal( mat, mat_comp )
        
        # Tests rclr
        # Calc RCLRs
        assay <- assay(tse, "counts")
        suppressWarnings(
        mat <- assays(mia::transformCounts(tse, assay.type = "counts",
                                           method = "rclr"))$rclr
        )
        suppressWarnings(
        mat_comp <- apply(as.matrix(assay), 2, FUN=function(x){
            temp <- log(x)
            temp[is.infinite(temp)] <- NA
            temp <- log(x) - mean(temp, na.rm = TRUE)
            temp[is.infinite(temp)] <- 0
            return(temp)
        })
        )
        # Round
        mat <- round(mat, 4)
        mat_comp <- round(mat_comp, 4)
        # Compare
        expect_equal(mat, mat_comp, check.attributes = FALSE)
        
        # Expect that error occurs
            expect_error(mia::transformCounts(tse, method = "clr"))
        # Expect that error does not occur
        tse <- mia::transformCounts(tse, method = "rclr")
        
        # Tests transformCounts, tries to calculate clr. Should be an error, because of zeros.
        expect_error(mia::transformCounts(tse, method = "clr"))
        
        # Tests that clr robust gives values that are approximately same if only 
        # one value per sample are changed to zero
        tse <- transformCounts(tse, method = "relabundance")
        # Adds pseudocount
        assay(tse, "test") <- assay(tse, "relabundance")+1
        assay(tse, "test2") <- assay(tse, "test") 
        # First row is zeroes
        assay(tse, "test2")[1, ] <- 0
        
        # clr robust transformations
        test <- assay(transformCounts(tse, method = "rclr", assay_name = "test"), "rclr")
        test2 <- assay(transformCounts(tse, method = "rclr", assay_name = "test2"), "rclr")
        
        # Removes first rows
        test <- test[-1, ]
        test2 <- test2[-1, ]
        
        # Expect that under 10 values are unequal. Values have only one decimal.
        expect_true( sum(round(test, 1) != round(test2, 1), na.rm = TRUE) < 10 )
        
        tse <- transformCounts(tse, method = "relabundance")
        # Expect error when counts and zeroes
        expect_error(
            transformCounts(tse, assay_name = "counts", method = "clr"))

        # Expect error warning when zeroes
        tse <- transformCounts(tse, method = "relabundance")
        expect_error(transformCounts(tse, assay_name = "relabundance", 
                                     method = "clr") ) 

        # Expect no warning when pseudocount is added, colSums are over 1
        tse <- transformCounts(
            tse, assay_name = "relabundance", method = "clr", pseudocount = 1)
        # Expect no warning when pseudocount is added, colSums are 1
        tse <- transformCounts(tse, method = "relabundance", pseudocount = 1, 
                                name = "relabund2")
        expect_error(transformCounts(tse, assay_name = "relabund2", method = "clr"))
        # Expect warning when colSums are not equal
        tse <- transformCounts(
            tse, assay_name = "counts", method = "clr", pseudocount = 1)
        tse <- transformCounts(tse, assay_name = "counts", method = "rclr")
        
        # Test that CLR with counts equal to CLR with relabundance
        assay(tse, "pseudo") <- assay(tse, "counts") + 1
        tse <- transformCounts(
            tse, assay_name = "counts", method = "relabundance")
        tse <- transformCounts(
            tse, assay_name = "pseudo", method = "clr", name = "clr1")
        tse <- transformCounts(
            tse, assay_name = "counts", method = "clr", name = "clr2",
            pseudocount =1)
        expect_equal(assay(tse, "clr1"), assay(tse, "clr2"),
                     check.attributes = FALSE)

        ############################# NAMES ####################################
        # Tests that samples have correct names
        expect_equal(colnames(assays(transformCounts(tse, assay.type = "relabundance",
                                                     method = "clr", pseudocount = 1))$clr),
                     colnames(assays(tse)$relabundance))

        # Tests that otus have correct names
        expect_equal(rownames(assays(transformCounts(tse, MARGIN = "features", method = "log10", pseudocount = 1000))$log10),
                     rownames(assays(tse)$counts))
        
        ################################## RANK ################################
        # Calculates rank
        tse_rank <- transformCounts(tse, method = "rank")
        # Expect that assay contains count and rank table
        expect_true( all(c("counts", "rank") %in% assayNames(tse_rank)) )
        
        for( i in c(1:10) ){
            # Gets columns from 'rank' table
            ranks <- assay(tse_rank, "rank")[,i]
            # Gets columns from 'counts' table, and calculates ranks
            counts_compare <- assay(tse_rank, "counts")[,i]
            counts_compare[counts_compare == 0] <- NA
            ranks_compare <- rank(counts_compare, na.last = "keep")
            ranks_compare[is.na(ranks_compare)] <- 0
            # Expect that they are equal
            expect_equal(ranks, ranks_compare, check.attributes = FALSE)
        }
        
        # Calculates rank with pseudocount
        tse_rank_pseudo <- transformCounts(tse, method = "rank", pseudocount = runif(1, 0, 1000))
        # Pseudocount should not change the rank
        expect_equal(tse_rank, tse_rank_pseudo, check.attributes = FALSE)

        # For other options for rank calculations, call the colRanks directly:
	    # data(esophagus); x <- esophagus;  
        # assay(x, "rank_average") <- t(colRanks(assay(x, "counts"), ties.method="average"))
        
        ############################## Z TRANSFORMATION ########################
        # Calculates Z-transformation for features	
        xx <- t(scale(t(as.matrix(assay(tse, "counts")))))
        expect_warning(z_assay <- assays(mia::ZTransform(tse, pseudocount = 1))$z)
        expect_equal(max(abs(z_assay - xx), na.rm=TRUE), 0,
                     tolerance = 1e-14, check.attributes = FALSE)
        
        # Tests that ZTransform and transformCounts(MARGIN = "features"),  gives the same result
        expect_warning(
        expect_equal(assays(mia::ZTransform(tse))$z,
                     assays(mia::transformCounts(tse, MARGIN = "features", method = "z"))$z)
        )
        
        # Test that transformations are equal to ones directly from vegan
        # clr
        tse <- transformCounts(tse, method = "relabundance")
        tse <- transformCounts(tse, assay_name = "relabundance", method = "clr",
                                pseudocount = 4)
        actual <- assay(tse, "clr")
        compare <- vegan::decostand(assay(tse, "relabundance"), method = "clr",
                                    pseudocount = 4, MARGIN = 2)
        expect_equal(actual, compare)
        # rclr
        tse <- transformCounts(tse, assay_name = "relabundance", method = "rclr")

        actual <- assay(tse, "rclr")
        compare <- vegan::decostand(assay(tse, "relabundance"), method = "rclr",
                                    MARGIN = 2)
        expect_equal(actual, compare)
        # alr
        tse <- transformCounts(tse, assay_name = "relabundance", method = "alr",
                               pseudocount = 4, reference = 2)
        actual <- assay(tse, "alr")
        compare <- vegan::decostand(assay(tse, "relabundance"), method = "alr",
                                    pseudocount = 4, reference = 2, MARGIN = 2)
        # Add reference sample
        sample <- rep(NA, nrow(tse))
        colnames <- c(colnames(tse)[2], colnames(compare))
        compare <- cbind(sample, compare)
        colnames(compare) <- colnames
        compare <- compare[, colnames(tse)]
        expect_equal(actual, compare, check.attributes = FALSE)
        # hellinger
        tse <- transformCounts(tse, assay_name = "counts", method = "hellinger",
                               pseudocount = 2, reference = 2)
        actual <- assay(tse, "hellinger")
        compare <- vegan::decostand(assay(tse, "counts"), method = "hellinger",
                                    pseudocount = 4, MARGIN = 2)
        expect_equal(actual, compare)
        # chi.squared
        tse <- transformCounts(tse, assay.type = "counts", method = "chi.square")
        actual <- assay(tse, "chi.square")
        compare <- vegan::decostand(assay(tse, "counts"), method = "chi.square",
                                    MARGIN = 2)
        compare <- t(compare)
        expect_equal(actual, compare)
    }

    # TSE object
    data(GlobalPatterns)
    tse <- GlobalPatterns
    testTransformations(GlobalPatterns)
    assay(tse,"counts") <- DelayedArray(assay(tse,"counts"))
    testTransformations(tse)
})

