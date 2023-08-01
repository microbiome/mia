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
                              clean.taxa.names = "\"")
    # Testing no prefixes removed
    expect_true(rowData(tse) %>%
                    apply(2,grepl,pattern="sk__|([dkpcofgs]+)__") %>%
                    all())
    # Testing no taxonomy ranks parsed
    expect_false(
        sapply(tolower(colnames(rowData(tse))),
               function(x) x %in% TAXONOMY_RANKS) %>% all())
    # Testing the clean.taxa.names, since the original artifact in the biom file 
    # is '\"'
    expect_false(apply(rowData(tse), 2, grepl, pattern="^\"") %>% all())
    
    # Testing prefixes removed
    tse <- makeTreeSEFromBiom(biom_object,
                              removeTaxaPrefixes=TRUE,
                              rankFromPrefix=FALSE,
                              clean.taxa.names = "\"")
    expect_false(rowData(tse) %>%
                    apply(2,grepl,pattern="sk__|([dkpcofgs]+)__") %>%
                    all())
    
    # Testing parsing taxonomy ranks from prefixes
    tse <- makeTreeSEFromBiom(biom_object,
                              removeTaxaPrefixes=FALSE,
                              rankFromPrefix=TRUE,
                              clean.taxa.names = "\"")
    expect_true(
        sapply(tolower(colnames(rowData(tse))),
               function(x) x %in% TAXONOMY_RANKS) %>% all())
    
    # Testing the clean.taxa.names, the original artifact in the biom file 
    # is '\"', as a test we rather try remove a non existing pattern.
    tse <- makeTreeSEFromBiom(biom_object,
                              removeTaxaPrefixes=FALSE,
                              rankFromPrefix=FALSE,
                              clean.taxa.names = "\\*|\\?")
    # with wrong pattern artifact not cleaned
    expect_true(apply(rowData(tse), 2, grepl, pattern="\"") %>% any())
    # Testing the clean.taxa.names, with the value 'auto' to automatically 
    # detect the artifact and remove it (in our case the artifact is '\"').
    tse <- makeTreeSEFromBiom(biom_object,
                              removeTaxaPrefixes=FALSE,
                              rankFromPrefix=FALSE,
                              clean.taxa.names = "auto")
    # Checking if 'auto' has detected and cleaned the artifact
    expect_false(apply(rowData(tse), 2, grepl, pattern="\"") %>% any())
    # Testing the clean.taxa.names, with the value NULL to not detect or clean 
    # anything.
    tse <- makeTreeSEFromBiom(biom_object,
                              removeTaxaPrefixes=FALSE,
                              rankFromPrefix=FALSE,
                              clean.taxa.names = NULL)
    # Checking if the '\"' artifact still exists.
    expect_true(apply(rowData(tse), 2, grepl, pattern="\"") %>% any())
    
    # General final test
    tse <- makeTreeSEFromBiom(biom_object,
                              removeTaxaPrefixes=TRUE,
                              rankFromPrefix=TRUE,
                              clean.taxa.names = 'auto')
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
    
    # General final test with another biom file
    biom_object <- biomformat::read_biom(
        system.file("extdata", "rich_dense_otu_table.biom",
                    package = "biomformat")
    )
    tse <- makeTreeSEFromBiom(biom_object,
                              removeTaxaPrefixes=TRUE,
                              rankFromPrefix=TRUE,
                              clean.taxa.names = 'auto')
    # check if taxa prefixes removed
    expect_false(rowData(tse) %>%
                     apply(2,grepl,pattern="sk__|([dkpcofgs]+)__") %>%
                     all())
    # Check if rank names were parsed correctly
    expect_true(
        sapply(tolower(colnames(rowData(tse))),
               function(x) x %in% TAXONOMY_RANKS) %>% all())
    
})