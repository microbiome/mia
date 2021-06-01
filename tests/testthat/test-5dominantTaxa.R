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

    }

    # TSE object
    data(GlobalPatterns)
    test_perSampleDominantTaxa(GlobalPatterns)


    test_that("countDominantTaxa", {

        test_countDominantTaxa <- function(tse){
            expect_equal(countDominantTaxa(tse, group = "SampleType")$dominant_taxa,
                         c("331820", "549656", "550960", "319044", "189047",
                           "279599", "329744", "12812",  "534609", "557211",
                           "87194", "484436", "64396", "98605", "256977",
                           "36155","71074",  "114821", "360229"))

            expect_equal(countDominantTaxa(tse,
                                           rank = "Kingdom",
                                           name = "dominant_kingdom")$dominant_kingdom[1],
                         c("Bacteria"))

            expect_equal(countDominantTaxa(tse, rank = "Order",
                                           name = "dominant_order")$rel.freq,
                         c(23.1, 11.5, 7.7, 7.7, 7.7, 7.7, 3.8, 3.8, 3.8, 3.8, 3.8, 3.8, 3.8, 3.8, 3.8))

            expect_equal(countDominantTaxa(tse, rank = "Class",
                                           name = "dominant_class")$rel.freq.pct,
                         c("23%", "15%", "12%", "8%", "8%", "8%", "4%", "4%", "4%", "4%", "4%", "4%", "4%"))

            # check sample type
            sample.type <- countDominantTaxa(tse, rank = "Class",
                                             group = "SampleType",
                                             name = "dominant_class")$SampleType

            expect_equal(as.character(sample.type),
                         c("Freshwater (creek)", "Mock", "Feces", "Feces", "Sediment (estuary)",
                           "Skin", "Freshwater", "Freshwater", "Ocean", "Ocean", "Ocean",
                           "Sediment (estuary)","Skin", "Soil", "Soil", "Soil", "Tongue", "Tongue"))

        }

        # TSE object
        data(GlobalPatterns)
        test_countDominantTaxa(GlobalPatterns)

    })
})

