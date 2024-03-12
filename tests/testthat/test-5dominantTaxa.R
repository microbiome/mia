context("perSampleDominantFeatures")

test_that("perSampleDominantFeatures", {

    test_perSampleDominantFeatures <- function(tse){

        # names
        exp.names.one <- c("CL3", "CC1", "SV1", "M31Fcsw",
                           "M11Fcsw","M31Plmr","M11Plmr",
                           "F21Plmr","M31Tong","M11Tong",
                           "LMEpi24M", "SLEpi20M", "AQC1cm",
                           "AQC4cm", "AQC7cm")
        # Generic test for output
        exp.vals.one <- c("36155", "256977", "71074",
                          "331820", "331820", "98605",
                          "484436", "64396", "360229",
                          "114821", "279599", "329744",
                          "549656","549656","549656")
        names(exp.vals.one) <- exp.names.one

        expect_equal(perSampleDominantFeatures(tse)[1:15], exp.vals.one)

        # Test at taxonomic level for values are passed to agglomerateRanks
        perSampleDominantFeatures(tse, rank = "Genus", na.rm = FALSE)

        exp.vals.two <- c("Genus:CandidatusSolibacter", "Genus:MC18",
                          "Class:Chloracidobacteria", "Genus:Bacteroides",
                          "Genus:Bacteroides", "Genus:Streptococcus",
                          "Family:Moraxellaceae", "Genus:Streptococcus",
                          "Genus:Neisseria", "Genus:Veillonella",
                          "Genus:Dolichospermum", "Family:ACK-M1",
                          "Order:Stramenopiles","Order:Stramenopiles","Order:Stramenopiles")
        names(exp.vals.two) <- exp.names.one
        expect_equal(perSampleDominantFeatures(tse,
                                           rank = "Genus",
                                           na.rm = FALSE)[1:15],
                     exp.vals.two)

        # Check if DominantTaxa is added to coldata
        expect_equal(colData(addPerSampleDominantFeatures(tse,
                                                      name="dominant"))$dominant[1:15],
                     exp.vals.one)
        expect_equal(colData(addPerSampleDominantFeatures(tse,
                                                      rank = "Genus",
                                                      na.rm = FALSE,
                                                      name="dominant"))$dominant[1:15],
                     exp.vals.two)
        
        tse1 <- tse
        # Now data contains 2 dominant taxa in one sample
        assay(tse1)[1, 1] <- max(assay(tse1)[, 1])
        
        # Get dominant taxa
        dominant_taxa <- perSampleDominantFeatures(tse, complete = TRUE)
        dominant_taxa1 <- perSampleDominantFeatures(tse1, complete = TRUE)
        
        # dominant_taxa1 should have one additional element
        expect_equal( length(unlist(dominant_taxa1)), length(dominant_taxa)+1)
        
        # Add dominant taxa to colData and check that it equals to dominant taxa
        # that is got by perSampleDominantFeatures
        add_dom <- unlist(colData(addPerSampleDominantFeatures(tse, complete = FALSE))$dominant_taxa)
        dominant_taxa <- perSampleDominantFeatures(tse, complete = FALSE)
        expect_equal(unname(add_dom), unname(dominant_taxa))
        
        # Test alias
        alias <- unlist(colData(addPerSampleDominantFeatures(tse))$dominant_taxa)
        expect_equal(unname(add_dom), unname(alias))
        alias <- perSampleDominantFeatures(tse)
        expect_equal(alias, dominant_taxa)
        
        # Test that chosen dominant taxon is included in the list of dominant taxa for those samples with multiple dominant taxa
        expect_warning(add_dom <- unlist(colData(addPerSampleDominantFeatures(tse1, complete = FALSE))$dominant_taxa))
        expect_true(add_dom[[1]] %in% dominant_taxa1)
    }

    # TSE object
    data(GlobalPatterns, package="mia")
    test_perSampleDominantFeatures(GlobalPatterns)


    test_that("countDominantFeatures", {

        test_countDominantFeatures <- function(tse){
            expect_equal(countDominantFeatures(tse, group = "SampleType")$dominant_taxa,
                         c("331820", "549656", "550960", "319044", "189047",
                           "279599", "329744", "12812",  "534609", "557211",
                           "87194", "484436", "64396", "98605", "256977",
                           "36155","71074",  "114821", "360229"))

            expect_equal(countDominantFeatures(tse,
                                           rank = "Kingdom")$dominant_taxa[1],
                         c("Bacteria"))

            expect_equal(countDominantFeatures(tse, rank = "Order", digits = 3)$rel_freq,
                         c(0.231, 0.115, 0.077, 0.077, 0.077, 0.077, 0.038, 0.038, 0.038, 0.038, 0.038, 0.038, 0.038, 0.038, 0.038))

            # check sample type
            sample.type <- countDominantFeatures(tse, rank = "Class",
                                             group = "SampleType")$SampleType

            expect_equal(as.character(sample.type),
                         c("Freshwater (creek)", "Mock", "Feces", "Feces", "Sediment (estuary)",
                           "Skin", "Freshwater", "Freshwater", "Ocean", "Ocean", "Ocean",
                           "Sediment (estuary)","Skin", "Soil", "Soil", "Soil", "Tongue", "Tongue"))
            
            tse1 <- tse
            # Now data contains 2 dominant taxa in one sample
            assay(tse1)[1, 1] <- max(assay(tse1)[, 1])
            
            # Calculate info about dominant taxa
            count_dominant <- countDominantFeatures(tse)
            count_dominant1 <- countDominantFeatures(tse1)
            
            # count_dominant1 should have one extra row, since there are more dominant taxa
            expect_equal( nrow(count_dominant1), nrow(count_dominant) + 1)
            
            
            # Test if dominant taxa list with only one dominant taxon is included in the list that have multiple for one sample
            expect_true(count_dominant$dominant_taxa[1] %in% count_dominant1$dominant_taxa)
            
            # Now the row lengths should be equal
            count_dominant <- countDominantFeatures(tse, complete = F)
            expect_warning(count_dominant1 <- countDominantFeatures(tse1, complete = F))
            
            expect_equal(nrow(count_dominant1), nrow(count_dominant))

        }

        # TSE object
        data(GlobalPatterns, package="mia")
        test_countDominantFeatures(GlobalPatterns)

    })
})

