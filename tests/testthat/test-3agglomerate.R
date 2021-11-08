context("agglomerate")
test_that("agglomerate", {
    gr <- GRanges("chr1",rep("1-6",11))
    df <- DataFrame(n = c(1:11))
    mcols(gr) <- df
    grl <- splitAsList(gr,1:11)
    mat <- matrix(1:110, nrow = 11)
    xtse <- TreeSummarizedExperiment(assays = list(mat = mat),
                                     rowRanges = unname(grl))
    tax_data <- DataFrame(Phylum = c(rep("a",3),rep("b",3),rep("c",3),rep("b",2)),
                          score = 1:11,
                          Family = c("c",NA,"d","e","f","g","h",NA,"h","e","f"),
                          n = 7:17)
    rowData(xtse) <- tax_data
    # mergeRows for agglomerateByRank
    tax_factors <- mia:::.get_tax_groups(xtse, col = 2)
    actual_family <- actual <- mergeRows(xtse, f = tax_factors)
    expect_s4_class(actual,class(xtse))
    expect_equal(dim(actual),c(8,10))
    expect_equal(assays(actual)$mat[8,1],c(8))
    expect_equal(assays(actual)$mat[7,1],c(16))
    tax_factors <- mia:::.get_tax_groups(xtse, col = 1)
    actual_phylum <- actual <- mergeRows(xtse, f = tax_factors)
    expect_s4_class(actual,class(xtse))
    expect_equal(dim(actual),c(3,10))
    expect_equal(assays(actual)$mat[1,1],c(6))
    expect_equal(assays(actual)$mat[2,1],c(36))
    expect_equal(assays(actual)$mat[3,1],c(24))
    #
    expect_error(agglomerateByRank(xtse,"",na.rm=FALSE),
                 "'rank' must be an non empty single character value")
    expect_error(agglomerateByRank(xtse,"Family",na.rm=""),
                 "'na.rm' must be TRUE or FALSE")
    expect_error(agglomerateByRank(xtse,"Family",na.rm=FALSE,agglomerateTree=""),
                 "'agglomerateTree' must be TRUE or FALSE")
    xtse2 <- xtse
    rowData(xtse2) <- NULL
    expect_error(agglomerateByRank(xtse2,"Family",na.rm=FALSE),
                 "taxonomyData needs to be populated")
    #
    actual <- agglomerateByRank(xtse,"Family",na.rm=FALSE)
    expect_equivalent(rowData(actual),rowData(actual_family))
    actual <- agglomerateByRank(xtse,"Phylum",na.rm=FALSE)
    expect_equivalent(rowData(actual),rowData(actual_phylum))
    #
    actual <- agglomerateByRank(xtse,"Family", na.rm = TRUE)
    expect_equal(dim(actual),c(6,10))
    expect_equal(rowData(actual)$Family,c("c","d","e","f","g","h"))
    actual <- agglomerateByRank(xtse,"Family", na.rm = FALSE) # the default
    expect_equal(dim(actual),c(8,10))
    expect_equal(rowData(actual)$Family,c("c",NA,"d","e","f","g","h",NA))
    actual <- agglomerateByRank(xtse,"Phylum")
    expect_equivalent(rowData(actual),rowData(actual_phylum))

    # Only one rank available in the object -
    # the same dimensionality is retained
    data(enterotype)
    expect_equal(length(unique(rowData(enterotype)[,"Genus"])),
                 nrow(agglomerateByRank(enterotype,"Genus")))

    # agglomeration in all its forms
    data(GlobalPatterns)
    se <- GlobalPatterns
    actual <- agglomerateByRank(se, rank = "Family")
    expect_equal(dim(actual),c(603,26))
    expect_equal(length(rowTree(actual)$tip.label),
                 length(rowTree(se)$tip.label))
    actual <- agglomerateByRank(se, rank = "Family", mergeTree = TRUE)
    expect_equal(dim(actual),c(603,26))
    expect_equal(length(rowTree(actual)$tip.label),
                 603)
    actual <- expect_warning(agglomerateByRank(se, rank = "Family",
                                               agglomerateTree = TRUE))
    expect_equal(dim(actual),c(603,26))
    expect_equal(length(rowTree(actual)$tip.label),
                 496)

    # checking reference consensus sequence generation
    se <- microbiomeDataSets::SilvermanAGutData()
    actual <- agglomerateByRank(se,"Genus", mergeRefSeq = FALSE)
    expect_equal(as.character(referenceSeq(actual)[[1]]),
                 paste0("TCAAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGTAGGCGGTTTGATAA",
                        "GTTAGAGGTGAAATCCCGGGGCTTAACTCCGGAACTGCCTCTAATACTGTTAG",
                        "ACTAGAGAGTAGTTGCGGTAGGCGGAATGTATGGTGTAGCGGTGAAATGCTTA",
                        "GAGATCATACAGAACACCGATTGCGAAGGCAGCTTACCAAACTATATCTGACG",
                        "TTGAGGCACGAAAGCGTGGGG"))
    actual <- agglomerateByRank(se,"Genus", mergeRefSeq = TRUE)
    expect_equal(as.character(referenceSeq(actual)[[1]]),
                 paste0("BCNMKCKTTVWYCKKMHTTMYTKKKYKTMMMKNKHDYKYMKDYKKNHNNNYMM",
                        "KHHNDNNKTKMMMDNBHNBKKCTYMMCHNBNDDDNKSSHBNNRWDMYKKBNND",
                        "NYTDRRKDVHNKNDRVGRNDRSBRRAWTBYNHRKKKWRSSRKKRAAWKSSKWR",
                        "RWDWTNDBRVRRAMHHCMRDKKSSRARGSSVSYYHNYBRRVHNDNNHYKRMVV",
                        "YKVRDNNNSRAARSBDKGGKK"))
})




