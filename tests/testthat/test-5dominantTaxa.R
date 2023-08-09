context("perSampleDominantTaxa")

test_that("perSampleDominantTaxa", {

    test_perSampleDominantTaxa <- function(tse){

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

        expect_equal(perSampleDominantTaxa(tse)[1:15], exp.vals.one)

        # Test at taxonomic level for values are passed to agglomerateRanks
        perSampleDominantTaxa(tse, rank = "Genus", na.rm = FALSE)

        exp.vals.two <- c("Genus:CandidatusSolibacter", "Genus:MC18",
                          "Class:Chloracidobacteria", "Genus:Bacteroides",
                          "Genus:Bacteroides", "Genus:Streptococcus",
                          "Family:Moraxellaceae", "Genus:Streptococcus",
                          "Genus:Neisseria", "Genus:Veillonella",
                          "Genus:Dolichospermum", "Family:ACK-M1",
                          "Order:Stramenopiles","Order:Stramenopiles","Order:Stramenopiles")
        names(exp.vals.two) <- exp.names.one
        expect_equal(perSampleDominantTaxa(tse,
                                           rank = "Genus",
                                           na.rm = FALSE)[1:15],
                     exp.vals.two)

        # Check if DominantTaxa is added to coldata
        expect_equal(colData(addPerSampleDominantTaxa(tse,
                                                      name="dominant"))$dominant[1:15],
                     exp.vals.one)
        expect_equal(colData(addPerSampleDominantTaxa(tse,
                                                      rank = "Genus",
                                                      na.rm = FALSE,
                                                      name="dominant"))$dominant[1:15],
                     exp.vals.two)
        
        tse1 <- tse
        # Now data contains 2 dominant taxa in one sample
        assay(tse1)[1, 1] <- max(assay(tse1)[, 1])
        
        # Get dominant taxa
        dominant_taxa <- perSampleDominantTaxa(tse)
        dominant_taxa1 <- perSampleDominantTaxa(tse1)
        
        # dominant_taxa1 should have one additioal element
        expect_equal( length(dominant_taxa1), length(dominant_taxa)+1 )
        
        # Remove additional dominant taxa
        dominant_taxa1_removed <- dominant_taxa1[ !(names(dominant_taxa1) == colnames(tse1)[1] & 
                                                dominant_taxa1 == rownames(tse1)[1]) ] 
        # Now they should be equal
        expect_equal(dominant_taxa1_removed, dominant_taxa)
        
        # Add dominant taxa to colData and check that it equals to dominant taxa
        # that is got by perSampleDominantTaxa
        add_dom1 <- unlist(colData(addPerSampleDominantTaxa(tse1))$dominant_taxa)
        expect_equal(unname(add_dom1), unname(dominant_taxa1))
        
        # Test alias
        alias <- unlist(colData(addPerSampleDominantFeatures(tse1))$dominant_taxa)
        expect_equal(unname(add_dom1), unname(alias))
        alias <- perSampleDominantFeatures(tse1)
        expect_equal(alias, dominant_taxa1)
    }

    # TSE object
    data(GlobalPatterns, package="mia")
    test_perSampleDominantTaxa(GlobalPatterns)


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

            expect_equal(countDominantFeatures(tse, rank = "Order", digits = 3)$rel.freq,
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
            
            # count_dominant should have one additional row
            expect_equal( nrow(count_dominant1), nrow(count_dominant)+1 )
            
            # Remove additional dominant taxa
            count_dominant1 <- count_dominant1[ !count_dominant1$dominant_taxa == rownames(tse1)[1], ]
            
            # Now the order of taxa should be equal
            expect_equal(count_dominant1$dominant_taxa, count_dominant$dominant_taxa)

        }

        # TSE object
        data(GlobalPatterns, package="mia")
        test_countDominantFeatures(GlobalPatterns)

    })
})

