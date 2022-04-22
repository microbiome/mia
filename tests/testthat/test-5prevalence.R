context("prevalence")

test_that("getPrevalence", {

    data(GlobalPatterns)
    expect_error(getPrevalence(GlobalPatterns, detection="test"),
                 "'detection' must be a single numeric value or coercible to one")
    expect_error(getPrevalence(GlobalPatterns, include_lowest="test"),
                 "'include_lowest' must be TRUE or FALSE")
    expect_error(getPrevalence(GlobalPatterns, sort="test"),
                 "'sort' must be TRUE or FALSE")
    expect_error(getPrevalence(GlobalPatterns, as_relative="test"),
                 "'as_relative' must be TRUE or FALSE")
    expect_error(getPrevalence(GlobalPatterns, abund_values="test"),
                 "'abund_values' must be a valid name")
    # Output should be always a frequency between 0 to 1
    pr <- getPrevalence(GlobalPatterns, detection=0.1/100, as_relative=TRUE)
    expect_true(min(pr) >= 0 && max(pr) <= 1)
    pr <- getPrevalence(GlobalPatterns, detection=0.1/100, as_relative=FALSE)
    expect_true(min(pr) >= 0 && max(pr) <= 1)

    # Same prevalences should be returned for as_relative T/F in certain cases.
    pr1 <- getPrevalence(GlobalPatterns, detection=1, include_lowest=TRUE, as_relative=FALSE)
    pr2 <- getPrevalence(GlobalPatterns, detection=0/100, include_lowest=FALSE, as_relative=TRUE)
    expect_true(all(pr1 == pr2))

    # Same prevalences should be returned for as_relative T/F in certain cases.
    pr1 <- getPrevalence(GlobalPatterns, detection=1, include_lowest=TRUE, as_relative=FALSE)
    pr2 <- getPrevalence(GlobalPatterns, detection=0, include_lowest=FALSE, as_relative=FALSE)
    expect_true(all(pr1 == pr2))

    # Different ways to use relative abundance should yield the same output
    pr2 <- getPrevalence(GlobalPatterns, as_relative=TRUE, abund_values = "counts")
    GlobalPatterns <- relAbundanceCounts(GlobalPatterns)
    pr1 <- getPrevalence(GlobalPatterns, as_relative=FALSE, abund_values = "relabundance")
    expect_true(all(pr1 == pr2))

    # Sorting should put the top values first
    pr <- getPrevalence(GlobalPatterns, sort=TRUE, detection = 0.1/100)
    expect_equal(as.vector(which.max(pr)), 1)
    pr <- names(head(getPrevalence(GlobalPatterns, sort=TRUE,  include_lowest = TRUE), 5L))
    actual <- getTopTaxa(GlobalPatterns,
                         method="prevalence",
                         top=5,
                         abund_values="counts")
    expect_equal(pr, actual)
    # Test alias
    alias<- getTopTaxa(GlobalPatterns,
                       method="prevalence",
                       top=5,
                       abund_values="counts")
    expect_equal(alias, actual)
    # Check that works also when rownames is NULL
    gp_null <- GlobalPatterns
    rownames(gp_null) <- NULL
    
    pr1 <- unname(getPrevalence(GlobalPatterns, detection=0.004, as_relative=TRUE))
    pr2 <- getPrevalence(gp_null, detection=0.004, as_relative=TRUE) 
    expect_equal(pr1, pr2)
    
    pr1 <- getPrevalence(GlobalPatterns, detection=0.004, as_relative=TRUE, rank = "Family")
    pr2 <- getPrevalence(gp_null, detection=0.004, as_relative=TRUE, rank = "Family") 
    expect_equal(pr1, pr2)
    
})


test_that("getPrevalentTaxa", {

    data(GlobalPatterns)
    expect_error(getPrevalentTaxa(GlobalPatterns, prevalence="test"),
                 "'prevalence' must be a single numeric value or coercible to one")
    # Results compatible with getPrevalence
    pr1 <- getPrevalentTaxa(GlobalPatterns, detection=0.1/100, as_relative=TRUE, sort=TRUE)
    pr2 <- names(getPrevalence(GlobalPatterns, rank = "Kingdom", detection=0.1/100, as_relative=TRUE, sort=TRUE))
    expect_true(all(pr1 == pr2))

    # Same sorting for toptaxa obtained in different ways
    pr1 <- getPrevalentTaxa(GlobalPatterns, detection=0.1/100, as_relative=TRUE, sort=TRUE)
    pr2 <- names(getPrevalence(GlobalPatterns, rank = "Kingdom", detection=0.1/100, as_relative=TRUE, sort=TRUE))
    expect_true(all(pr1 == pr2))

    # Retrieved taxa are the same for counts and relative abundances
    pr1 <- getPrevalentTaxa(GlobalPatterns, prevalence=0.1/100, as_relative=TRUE)
    pr2 <- getPrevalentTaxa(GlobalPatterns, prevalence=0.1/100, as_relative=FALSE)
    expect_true(all(pr1 == pr2))

    # Prevalence and detection threshold at 0 has the same impact on counts and relative abundances
    pr1 <- getPrevalentTaxa(GlobalPatterns, detection=0, prevalence=0, as_relative=TRUE)
    pr2 <- getPrevalentTaxa(GlobalPatterns, detection=0, prevalence=0, as_relative=FALSE)
    expect_true(all(pr1 == pr2))
    
    # Check that works also when rownames is NULL
    gp_null <- GlobalPatterns
    rownames(gp_null) <- NULL
    
    pr1 <- getPrevalentTaxa(GlobalPatterns, detection=0.0045, prevalence = 0.25, as_relative=TRUE)
    pr1 <- which(rownames(GlobalPatterns) %in% pr1)
    pr2 <- getPrevalentTaxa(gp_null, detection=0.0045, prevalence = 0.25, as_relative=TRUE) 
    expect_equal(pr1, pr2)
    
    # Test alias
    alias <- getPrevalentFeatures(gp_null, detection=0.0045, prevalence = 0.25, as_relative=TRUE) 
    expect_equal(pr1, alias)
    
    pr1 <- getPrevalentTaxa(GlobalPatterns, detection=0.004, prevalence = 0.1, 
                            as_relative=TRUE, rank = "Family")
    pr2 <- getPrevalentTaxa(gp_null, detection=0.004, prevalence = 0.1, 
                            as_relative=TRUE, rank = "Family") 
    expect_equal(pr1, pr2)

})

test_that("getRareTaxa", {

    data(GlobalPatterns)
    expect_error(getRareTaxa(GlobalPatterns, prevalence="test"),
                 "'prevalence' must be a single numeric value or coercible to one")

    ############# Test that output type is correct #############
    expect_type(taxa <- getRareTaxa(GlobalPatterns,
                                    detection = 130,
                                    prevalence = 90/100), "character")

    ##### Test that getPrevalentTaxa and getRareTaxa has all the taxa ####

    # Gets rownames for all the taxa
    all_taxa <- rownames(GlobalPatterns)

    # Gets prevalent taxa
    prevalent_taxa <- getPrevalentTaxa(GlobalPatterns,
                                       detection = 0,
                                       prevalence = 90/100,
                                       include_lowest = FALSE)
    # Gets rare taxa
    rare_taxa <- getRareTaxa(GlobalPatterns,
                             detection = 0,
                             prevalence = 90/100,
                             include_lowest = FALSE)

    # Concatenates prevalent and rare taxa
    prevalent_and_rare_taxa <- c(prevalent_taxa, rare_taxa)

    # If all the elements are in another vector, vector of TRUEs
    # Opposite --> negative of FALSEs
    # If every element is FALSE, any is FALSE --> opposite --> TRUE
    expect_true( !any( !(all_taxa %in% prevalent_and_rare_taxa) ) &&
                     !any( !( prevalent_and_rare_taxa %in% all_taxa) ) )

    ##### Test that getPrevalentTaxa and getRareTaxa has all the taxa, ####
    ##### but now with detection limit and with rank #####

    # Check that it works with all the ranks
    ranks <- taxonomyRanks(GlobalPatterns)

    for( rank in ranks ){

        # Agglomerates data by rank
        se <- agglomerateByRank(GlobalPatterns, rank = rank)

        # Gets rownames for all the taxa
        all_taxa <- rownames(se)

        # All but "Kingdom" can includes different taxa levels and e.g. "Species"
        # before their name
        if( rank != "Kingdom" ){
            # Takes e.g. only species and removes e.g. "Species:" from the names
            all_taxa <- stringr::str_remove(all_taxa[grepl(paste0(rank, ":"), all_taxa)], paste0(rank, ":"))
        }

        # Gets prevalent taxa
        prevalent_taxa <- getPrevalentTaxa(GlobalPatterns,
                                           prevalence = 0.05,
                                           detection = 0.1,
                                           rank = rank,
                                           include_lowest = TRUE, as_relative = TRUE)
        # Gets rare taxa
        rare_taxa <- getRareTaxa(GlobalPatterns,
                                 prevalence = 0.05,
                                 detection = 0.1,
                                 rank = rank,
                                 include_lowest = TRUE, as_relative = TRUE)

        # Concatenates prevalent and rare taxa
        prevalent_and_rare_taxa <- c(prevalent_taxa, rare_taxa)

        # If all the elements are in another vector, vector of TRUEs
        # Opposite --> negative of FALSEs
        # If every element is FALSE, any is FALSE --> opposite --> TRUE
        expect_true( !any( !(all_taxa %in% prevalent_and_rare_taxa) ) &&
                         !any( !( prevalent_and_rare_taxa %in% all_taxa) ) )

        # Expect that there are no duplicates
        expect_true(!anyDuplicated(prevalent_taxa))
        expect_true(!anyDuplicated(rare_taxa))

    }
    
    # Check that works also when rownames is NULL
    gp_null <- GlobalPatterns
    rownames(gp_null) <- NULL
    
    pr1 <- getRareTaxa(GlobalPatterns, detection=0.0045, prevalence = 0.25, as_relative=TRUE)
    pr1 <- which(rownames(GlobalPatterns) %in% pr1)
    pr2 <- getRareTaxa(gp_null, detection=0.0045, prevalence = 0.25, as_relative=TRUE) 
    expect_equal(pr1, pr2)
    
    # Test lias
    alias <- getRareFeatures(gp_null, detection=0.0045, prevalence = 0.25, as_relative=TRUE)
    expect_equal(pr1, alias)
    
    pr1 <- getRareTaxa(GlobalPatterns, detection=0.004, prevalence = 0.1, 
                       as_relative=TRUE, rank = "Family")
    pr2 <- getRareTaxa(gp_null, detection=0.004, prevalence = 0.1, 
                       as_relative=TRUE, rank = "Family") 
    expect_equal(pr1, pr2)

})

test_that("subsetByPrevalentTaxa", {
    data(GlobalPatterns)
    expect_error(subsetByPrevalentTaxa(GlobalPatterns, prevalence="test"),
                 "'prevalence' must be a single numeric value or coercible to one")
    # Expect TSE object
    expect_equal(class(subsetByPrevalentTaxa(GlobalPatterns)), class(GlobalPatterns))
    
    # Results compatible with getPrevalentTaxa
    pr1 <- rownames(subsetByPrevalentTaxa(GlobalPatterns, rank = "Class", detection=0.1/100, 
                                          as_relative=TRUE, sort=TRUE))
    pr2 <- getPrevalentTaxa(GlobalPatterns, rank = "Class", detection=0.1/100, 
                            as_relative=TRUE, sort=TRUE)
    expect_true(all(pr1 == pr2))
    
    # Retrieved taxa are the same for counts and relative abundances
    pr1 <- assay(subsetByPrevalentTaxa(GlobalPatterns, prevalence=0.1/100, as_relative=TRUE), "counts")
    pr2 <- assay(subsetByPrevalentTaxa(GlobalPatterns, prevalence=0.1/100, as_relative=FALSE), "counts")
    expect_true(all(pr1 == pr2))
    
    # Prevalence and detection threshold at 0 has the same impact on counts and relative abundances
    pr1 <- rownames(subsetByPrevalentTaxa(GlobalPatterns, detection=0, prevalence=0, as_relative=TRUE))
    pr2 <- rownames(subsetByPrevalentTaxa(GlobalPatterns, detection=0, prevalence=0, as_relative=FALSE))
    expect_true(all(pr1 == pr2))
    
    # Check that works also when rownames is NULL
    gp_null <- GlobalPatterns
    rownames(gp_null) <- NULL
    
    pr1 <- subsetByPrevalentTaxa(GlobalPatterns, detection=12, prevalence = 0.33)
    pr1 <- unname(assay(pr1, "counts"))
    pr2 <- subsetByPrevalentTaxa(gp_null, detection=12, prevalence = 0.33) 
    pr2 <- unname(assay(pr2, "counts"))
    expect_equal(pr1, pr2)
    
    pr1 <- subsetByPrevalentTaxa(GlobalPatterns, detection=5, prevalence = 0.33, rank = "Phylum")
    pr1 <- unname(assay(pr1, "counts"))
    pr2 <- subsetByPrevalentTaxa(gp_null, detection=5, prevalence = 0.33, rank = "Phylum") 
    pr2 <- unname(assay(pr2, "counts"))
    expect_equal(pr1, pr2)
    
    # Test alias
    alias <- subsetByPrevalentFeatures(gp_null, detection=5, prevalence = 0.33, rank = "Phylum") 
    alias <- unname(assay(alias, "counts"))
    expect_equal(alias, pr2)
    
})

test_that("subsetByRareTaxa", {
    data(GlobalPatterns)
    expect_error(subsetByRareTaxa(GlobalPatterns, prevalence="test"),
                 "'prevalence' must be a single numeric value or coercible to one")
    # Expect TSE object
    expect_equal(class(subsetByRareTaxa(GlobalPatterns)), class(GlobalPatterns))
    
    # Results compatible with getRareTaxa
    pr1 <- rownames(subsetByRareTaxa(GlobalPatterns, rank = "Phylum", detection=0.1/100, 
                                     as_relative=TRUE, sort=TRUE))
    pr2 <- getRareTaxa(GlobalPatterns, rank = "Phylum", detection=0.1/100, 
                       as_relative=TRUE, sort=TRUE)
    expect_true(all(pr1 == pr2))
    
    # Retrieved taxa are the same for counts and relative abundances
    pr1 <- assay(subsetByRareTaxa(GlobalPatterns, prevalence=0.1/100, as_relative=TRUE), "counts")
    pr2 <- assay(subsetByRareTaxa(GlobalPatterns, prevalence=0.1/100, as_relative=FALSE), "counts")
    expect_true(all(pr1 == pr2))
    
    # Prevalence and detection threshold at 0 has the same impact on counts and relative abundances
    pr1 <- rownames(subsetByRareTaxa(GlobalPatterns, detection=0, prevalence=0, as_relative=TRUE))
    pr2 <- rownames(subsetByRareTaxa(GlobalPatterns, detection=0, prevalence=0, as_relative=FALSE))
    expect_true(all(pr1 == pr2))
    
    # subsetByRareTaxa + subsetByPrevalentTaxa should include all the taxa in OTU level
    d <- runif(1, 0.0001, 0.1)
    p <- runif(1, 0.0001, 0.5)
    rare <- rownames(subsetByRareTaxa(GlobalPatterns, detection=d, prevalence=p, 
                                      as_relative=TRUE))
    
    prevalent <- rownames(subsetByPrevalentTaxa(GlobalPatterns, detection=d, prevalence=p, 
                                                as_relative=TRUE))
    
    all_taxa <- c(rare, prevalent)
    
    expect_true( all(all_taxa %in% rownames(GlobalPatterns)) )
    
    # Check that works also when rownames is NULL
    gp_null <- GlobalPatterns
    rownames(gp_null) <- NULL
    
    pr1 <- subsetByRareTaxa(GlobalPatterns, detection=12, prevalence = 0.33)
    pr1 <- unname(assay(pr1, "counts"))
    pr2 <- subsetByRareTaxa(gp_null, detection=12, prevalence = 0.33) 
    pr2 <- unname(assay(pr2, "counts"))
    expect_equal(pr1, pr2)
    
    pr1 <- subsetByRareTaxa(GlobalPatterns, detection=5, prevalence = 0.33, rank = "Phylum")
    pr1 <- unname(assay(pr1, "counts"))
    pr2 <- subsetByRareTaxa(gp_null, detection=5, prevalence = 0.33, rank = "Phylum") 
    pr2 <- unname(assay(pr2, "counts"))
    expect_equal(pr1, pr2)
    
    # Test alias
    alias <- subsetByRareFeatures(gp_null, detection=5, prevalence = 0.33, rank = "Phylum") 
    alias <- unname(assay(alias, "counts"))
    expect_equal(alias, pr2)
    
})

test_that("agglomerateByPrevalence", {

    data(GlobalPatterns)
    expect_error(agglomerateByPrevalence(GlobalPatterns, other_label=TRUE),
                 "'other_label' must be a single character value")
    actual <- agglomerateByPrevalence(GlobalPatterns)
    expect_s4_class(actual,class(GlobalPatterns))
    expect_equal(dim(actual),c(2,26))

    actual <- agglomerateByPrevalence(GlobalPatterns,
                                      rank = "Phylum",
                                      detection = 1/100,
                                      prevalence = 50/100,
                                      as_relative = TRUE,
                                      other_label = "test")
    expect_s4_class(actual,class(GlobalPatterns))
    expect_equal(dim(actual),c(6,26))
    expect_equal(rowData(actual)[6,"Phylum"],"test")

    actual <- agglomerateByPrevalence(GlobalPatterns,
                                      rank = NULL,
                                      detection = 0.0001,
                                      prevalence = 50/100,
                                      as_relative = TRUE,
                                      other_label = "test")
    expect_s4_class(actual,class(GlobalPatterns))
    expect_equal(dim(actual),c(6,26))
    expect_true(all(is.na(rowData(actual)[6,])))
})
