context("makeTreeSEFromBiom")
test_that("makeTreeSEFromBiom", {
    biom_object <- biomformat::read_biom(
        system.file("extdata/testdata/Aggregated_humanization2.biom",
                    package="mia")
        )
    library(mia)
    library(dplyr)
    
    tse <- makeTreeSEFromBiom(biom_object,
                              removeTaxaPrefixes=FALSE,
                              rankFromPrefix=FALSE,
                              cleanTaxaPattern = "\"")
    # Testing no prefixes removed
    expect_true(rowData(tse) %>%
                    apply(2,grepl,pattern="sk__|([dkpcofgs]+)__") %>%
                    all())
    # Testing no taxonomy ranks parsed
    expect_false(
        sapply(tolower(colnames(rowData(tse))),
               function(x) x %in% TAXONOMY_RANKS) %>% all())
    # Testing the cleanTaxaPattern, since the original artifact in the biom file 
    # is '\"'
    expect_false(apply(rowData(tse), 2, grepl, pattern="^\"") %>% all())
    
    # Testing prefixes removed
    tse <- makeTreeSEFromBiom(biom_object,
                              removeTaxaPrefixes=TRUE,
                              rankFromPrefix=FALSE,
                              cleanTaxaPattern = "\"")
    expect_false(rowData(tse) %>%
                    apply(2,grepl,pattern="sk__|([dkpcofgs]+)__") %>%
                    all())
    
    # Testing parsing taxonomy ranks from prefixes
    tse <- makeTreeSEFromBiom(biom_object,
                              removeTaxaPrefixes=FALSE,
                              rankFromPrefix=TRUE,
                              cleanTaxaPattern = "\"")
    expect_true(
        sapply(tolower(colnames(rowData(tse))),
               function(x) x %in% TAXONOMY_RANKS) %>% all())
    
    # Testing the cleanTaxaPattern, the original artifact in the biom file 
    # is '\"', as a test we rather try remove a non existing pattern.
    tse <- makeTreeSEFromBiom(biom_object,
                              removeTaxaPrefixes=FALSE,
                              rankFromPrefix=FALSE,
                              cleanTaxaPattern = "\\*|\\?")
    # with wrong pattern artifact not cleaned
    expect_true(apply(rowData(tse), 2, grepl, pattern="\"") %>% any())
    # Testing the cleanTaxaPattern, with the value 'auto' to automatically 
    # detect the artifact and remove it (in our case the artifact is '\"').
    tse <- makeTreeSEFromBiom(biom_object,
                              removeTaxaPrefixes=FALSE,
                              rankFromPrefix=FALSE,
                              cleanTaxaPattern = "auto")
    # Checking if 'auto' has detected and cleaned the artifact
    expect_false(apply(rowData(tse), 2, grepl, pattern="\"") %>% any())
    # Testing the cleanTaxaPattern, with the value NULL to not detect or clean 
    # anything.
    tse <- makeTreeSEFromBiom(biom_object,
                              removeTaxaPrefixes=FALSE,
                              rankFromPrefix=FALSE,
                              cleanTaxaPattern = NULL)
    # Checking if the '\"' artifact still exists.
    expect_true(apply(rowData(tse), 2, grepl, pattern="\"") %>% any())
    
    # General final test
    tse <- makeTreeSEFromBiom(biom_object,
                              removeTaxaPrefixes=TRUE,
                              rankFromPrefix=TRUE,
                              cleanTaxaPattern = 'auto')
    # check if '\"' cleaned
    expect_false(apply(rowData(tse), 2, grepl, pattern="\"") %>% any())
    # check if taxa prefixes removed
    expect_false(rowData(tse) %>%
                     apply(2,grepl,pattern="sk__|([dkpcofgs]+)__") %>%
                     all())
    # Check if rank names were parsed correctly
    expect_true(
        sapply(tolower(colnames(rowData(tse))),
               function(x) x %in% TAXONOMY_RANKS) %>% all())
})