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
    # mergeRows for mergeFeaturesByRank
    tax_factors <- mia:::.get_tax_groups(xtse, col = 2)
    actual_family <- actual <- mergeRows(xtse, f = tax_factors)
    expect_s4_class(actual,class(xtse))
    expect_equal(dim(actual),c(8,10))
    expect_equal(assays(actual)$mat[8,1],c(c_NA = 8))
    expect_equal(assays(actual)$mat[7,1],c(c_h = 16))
    tax_factors <- mia:::.get_tax_groups(xtse, col = 1)
    actual_phylum <- actual <- mergeRows(xtse, f = tax_factors)
    expect_s4_class(actual,class(xtse))
    expect_equal(dim(actual),c(3,10))
    expect_equal(assays(actual)$mat[1,1],c(a = 6))
    expect_equal(assays(actual)$mat[2,1],c(b = 36))
    expect_equal(assays(actual)$mat[3,1],c(c = 24))
    #
    expect_error(mergeFeaturesByRank(xtse,"",na.rm=FALSE),
                 "'rank' must be an non empty single character value")
    expect_error(mergeFeaturesByRank(xtse,"Family",na.rm=""),
                 "'na.rm' must be TRUE or FALSE")
    expect_error(
        mergeFeaturesByRank(xtse,"Family",na.rm=FALSE,agglomerate.tree=""),
        "'agglomerate.tree' must be TRUE or FALSE")
    xtse2 <- xtse
    rowData(xtse2) <- NULL
    expect_error(mergeFeaturesByRank(xtse2,"Family",na.rm=FALSE),
                 "taxonomyData needs to be populated")
    #
    actual <- mergeFeaturesByRank(xtse,"Family",na.rm=FALSE)
    expect_equivalent(rowData(actual),rowData(actual_family))
    actual <- mergeFeaturesByRank(xtse,"Phylum",na.rm=FALSE)
    expect_equivalent(rowData(actual),rowData(actual_phylum))
    #
    actual <- mergeFeaturesByRank(xtse,"Family", na.rm = TRUE)
    expect_equal(dim(actual),c(6,10))
    expect_equal(rowData(actual)$Family,c("c","d","e","f","g","h"))
    actual <- mergeFeaturesByRank(xtse,"Family", na.rm = FALSE) # the default
    expect_equal(dim(actual),c(8,10))
    expect_equal(rowData(actual)$Family,c("c",NA,"d","e","f","g","h",NA))
    actual <- mergeFeaturesByRank(xtse,"Phylum")
    expect_equivalent(rowData(actual),rowData(actual_phylum))
    #
    actual1 <- mergeFeaturesByRank(xtse,"Family")
    actual2 <- .merge_features(xtse, merge.by = "Family")
    expect_equal(actual1, actual2)
    expect_equal(mergeFeaturesByRank(xtse,"Family"), mergeFeaturesByRank(xtse,"Family"))

    # Only one rank available in the object -
    # the same dimensionality is retained
    data(enterotype, package="mia")
    expect_equal(length(unique(rowData(enterotype)[,"Genus"])),
                 nrow(mergeFeaturesByRank(enterotype,"Genus")))

    # agglomeration in all its forms
    data(GlobalPatterns, package="mia")
    se <- GlobalPatterns
    actual <- mergeFeaturesByRank(se, rank = "Family")
    expect_equal(dim(actual),c(603,26))
    expect_equal(length(rowTree(actual)$tip.label),
                 length(rowTree(se)$tip.label))
    actual <- mergeFeaturesByRank(se, rank = "Family", agglomerate.tree = TRUE)
    expect_equal(dim(actual),c(603,26))
    expect_equal(length(rowTree(actual)$tip.label), 603)
    actual <- mergeFeaturesByRank(se, rank = "Family", agglomerate.tree = TRUE)
    expect_equal(dim(actual),c(603,26))
    expect_equal(length(rowTree(actual)$tip.label), nrow(actual))
    # Test that warning occurs when assay contian binary or negative values
    se1 <- transformAssay(se, method = "pa")
    se2 <- se1
    assay(se2, "pa")[1, 1] <- -1
    expect_warning(mergeFeaturesByRank(se1, rank = "Phylum"))
    expect_warning(mergeFeaturesByRank(se1, rank = "Order"))

    # Load data from miaTime package
    skip_if_not(require("miaTime", quietly = TRUE))
    data(SilvermanAGutData)
    se <- SilvermanAGutData
    # checking reference consensus sequence generation
    actual <- mergeFeaturesByRank(se,"Genus", mergeRefSeq = FALSE)
    expect_equal(as.character(referenceSeq(actual)[[1]]),
                 paste0("TCAAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGTAGGCGGTTTGATAA",
                        "GTTAGAGGTGAAATCCCGGGGCTTAACTCCGGAACTGCCTCTAATACTGTTAG",
                        "ACTAGAGAGTAGTTGCGGTAGGCGGAATGTATGGTGTAGCGGTGAAATGCTTA",
                        "GAGATCATACAGAACACCGATTGCGAAGGCAGCTTACCAAACTATATCTGACG",
                        "TTGAGGCACGAAAGCGTGGGG"))
    actual <- mergeFeaturesByRank(se,"Genus", mergeRefSeq = TRUE)
    expect_equal(as.character(referenceSeq(actual)[[1]]),
                 paste0("BCNMKCKTTVWYCKKMHTTMYTKKKYKTMMMKNKHDYKYMKDYKKNHNNNYMM",
                        "KHHNDNNKTKMMMDNBHNBKKCTYMMCHNBNDDDNKSSHBNNRWDMYKKBNND",
                        "NYTDRRKDVHNKNDRVGRNDRSBRRAWTBYNHRKKKWRSSRKKRAAWKSSKWR",
                        "RWDWTNDBRVRRAMHHCMRDKKSSRARGSSVSYYHNYBRRVHNDNNHYKRMVV",
                        "YKVRDNNNSRAARSBDKGGKK"))
    # Test that remove_empty_ranks work
    expect_error(mergeFeaturesByRank(se, rank = "Class", remove_empty_ranks = NULL))
    expect_error(mergeFeaturesByRank(se, rank = "Class", remove_empty_ranks = "NULL"))
    expect_error(mergeFeaturesByRank(se, rank = "Class", remove_empty_ranks = 1))
    expect_error(mergeFeaturesByRank(se, rank = "Class", remove_empty_ranks = c(TRUE, TRUE)))
    x <- mergeFeaturesByRank(se, rank = "Class")
    rd1 <- rowData(x)[, 1:3]
    x <- mergeFeaturesByRank(se, rank = "Class", remove_empty_ranks = TRUE)
    rd2 <- rowData(x)
    expect_equal(rd1, rd2)
    # Test that make_unique work
    uniq <- mergeFeaturesByRank(se, rank = "Species")
    not_uniq <- mergeFeaturesByRank(se, rank = "Species", make_unique = FALSE)
    expect_true( !any( duplicated(rownames(uniq)) ) )
    expect_true( any( duplicated(rownames(not_uniq)) ) )
})




