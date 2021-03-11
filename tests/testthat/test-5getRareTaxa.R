context("getRareTaxa")

test_that("getRareTaxa", {

    # TSE object
    data("GlobalPatterns")

    ############# Test that output type is correct #############
    expect_type(taxa <- getRareTaxa(GlobalPatterns,
                                    detection = 130,
                                    prevalence = 90/100), "character")

    expect_type(taxa <- getRare(GlobalPatterns,
                                detection = 130,
                                prevalence = 90/100), typeof(GlobalPatterns))
    ############################################################

    ### Test that getRare has same rownames than getRareTaxa ###

    expect_equal(getRareTaxa(GlobalPatterns, rank = "Phylum"),
                 rownames(getRare(GlobalPatterns, rank = "Phylum")))

    ############################################################

    ##### Test that getPrevalentTaxa and getRareTaxa has all the taxa ####

    # Gets rownames for all the taxa
    agglomerated <- agglomerateByRank(GlobalPatterns, rank = "Species")
    all_taxa <- rownames(agglomerated)

    # Gets prevalent taxa
    prevalent_taxa <- getPrevalentTaxa(GlobalPatterns, rank = "Species",
                                       detection = 0.125,
                                       prevalence = 90/100,
                                       as_relative = TRUE,
                                       include_lowest = FALSE)
    # Gets rare taxa
    rare_taxa <- getRareTaxa(GlobalPatterns, rank = "Species",
                                       detection = 0.125,
                                       prevalence = 90/100,
                                       as_relative = TRUE,
                                       include_lowest = TRUE)
    # Concatenates prevalent and rare taxa
    prevalent_and_rare_taxa <- c(prevalent_taxa, rare_taxa)

    # Not working correctly, this is not equal to all the taxa###########################
    expect_mapequal(all_taxa, prevalent_and_rare_taxa)





})

