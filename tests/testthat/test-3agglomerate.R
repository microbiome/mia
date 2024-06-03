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
    actual_family <- actual <- agglomerateByVariable(xtse, MARGIN = "rows",
                                                    f = tax_factors)
    expect_s4_class(actual,class(xtse))
    expect_equal(dim(actual),c(8,10))
    expect_equal(assays(actual)$mat[8,1],c(c_NA = 8))
    expect_equal(assays(actual)$mat[7,1],c(c_h = 16))
    tax_factors <- mia:::.get_tax_groups(xtse, col = 1)
    actual_phylum <- actual <- agglomerateByVariable(xtse, MARGIN = "rows", 
                                                    f = tax_factors)
    expect_s4_class(actual,class(xtse))
    expect_equal(dim(actual),c(3,10))
    expect_equal(assays(actual)$mat[1,1],c(a = 6))
    expect_equal(assays(actual)$mat[2,1],c(b = 36))
    expect_equal(assays(actual)$mat[3,1],c(c = 24))
    #
    expect_error(agglomerateByRank(xtse,"",na.rm=FALSE),
                 "'rank' must be a non-empty single character value")
    expect_error(agglomerateByRank(xtse,"Family",na.rm=""),
                 "'na.rm' must be TRUE or FALSE")
    expect_error(
        agglomerateByRank(xtse,"Family",na.rm=FALSE,agglomerate.tree=""),
        "'agglomerate.tree' must be TRUE or FALSE")
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
    actual <- agglomerateByRank(xtse,"Family", onRankOnly = FALSE, na.rm = TRUE)
    expect_equal(dim(actual),c(6,10))
    expect_equal(rowData(actual)$Family,c("c","d","e","f","g","h"))
    actual <- agglomerateByRank(xtse,"Family", onRankOnly = FALSE, na.rm = FALSE) # the default
    expect_equal(dim(actual),c(8,10))
    expect_equal(rowData(actual)$Family,c("c","d","e","f","g","h",NA,NA))
    actual <- agglomerateByRank(xtse,"Phylum")
    expect_equivalent(rowData(actual),rowData(actual_phylum))
    #
    actual1 <- agglomerateByRank(xtse,"Family")
    actual2 <- .merge_features(xtse, merge.by = "Family")
    expect_equal(actual1, actual2)
    expect_equal(agglomerateByRank(xtse,"Family"), agglomerateByRank(xtse,"Family"))

    # Only one rank available in the object -
    # the same dimensionality is retained
    data(enterotype, package="mia")
    expect_equal(length(unique(rowData(enterotype)[,"Genus"])),
                 nrow(agglomerateByRank(enterotype,"Genus", onRankOnly = FALSE, 
                 na.rm = FALSE)))

    # agglomeration in all its forms
    data(GlobalPatterns, package="mia")
    se <- GlobalPatterns
    actual <- agglomerateByRank(se, rank = "Family", 
        onRankOnly = FALSE, na.rm = FALSE)
    expect_equal(dim(actual),c(603,26))
    expect_equal(length(rowTree(actual)$tip.label),
                 length(rowTree(se)$tip.label))
    actual <- agglomerateByRank(se, rank = "Family", 
        onRankOnly = FALSE, na.rm = FALSE, agglomerate.tree = TRUE)
    expect_equal(dim(actual),c(603,26))
    expect_equal(length(rowTree(actual)$tip.label), 603)
    actual <- agglomerateByRank(se, rank = "Family", 
        onRankOnly = FALSE, na.rm = FALSE, agglomerate.tree = TRUE)
    expect_equal(dim(actual),c(603,26))
    expect_equal(length(rowTree(actual)$tip.label), nrow(actual))
    # Test that warning occurs when assay contian binary or negative values
    se1 <- transformAssay(se, method = "pa")
    se2 <- se1
    assay(se2, "pa")[1, 1] <- -1
    expect_warning(agglomerateByRank(se1, rank = "Phylum"))
    expect_warning(agglomerateByRank(se1, rank = "Order"))

    # Load data from miaTime package
    skip_if_not(require("miaTime", quietly = TRUE))
    data(SilvermanAGutData)
    se <- SilvermanAGutData
    # checking reference consensus sequence generation
    actual <- agglomerateByRank(se,"Genus", mergeRefSeq = FALSE)
    # There should be only one exact match for each sequence
    seqs_test <- as.character( referenceSeq(actual) )
    seqs_ref <- as.character( referenceSeq(se) )
    expect_true(all(vapply(
        seqs_test, function(seq) sum(seqs_ref %in% seq) == 1,
        FUN.VALUE = logical(1) )) )
    # Merging creates concensus sequences.
    th <- runif(1, 0, 1)
    actual <- agglomerateByRank(
        se, "Genus", mergeRefSeq = TRUE, threshold = th)
    seqs_test <- referenceSeq(actual)
    # Get single taxon as reference. Merge those sequences and test that it
    # equals to one that is output of agglomerateByRank
    seqs_ref <- referenceSeq(se)
    feature <- sample(na.omit(rowData(se)[["Genus"]]), 1)
    seqs_ref <- seqs_ref[ rowData(se)[["Genus"]] %in% feature ]
    seqs_ref <- .merge_refseq(
        seqs_ref, factor(rep(feature, length(seqs_ref))), rownames(seqs_ref),
        threshold = th)
    seqs_test <- seqs_test[ names(seqs_test) %in% feature ]
    expect_equal(seqs_test, seqs_ref)
    # Test that remove_empty_ranks work
    expect_error(agglomerateByRank(se, rank = "Class", remove_empty_ranks = NULL))
    expect_error(agglomerateByRank(se, rank = "Class", remove_empty_ranks = "NULL"))
    expect_error(agglomerateByRank(se, rank = "Class", remove_empty_ranks = 1))
    expect_error(agglomerateByRank(se, rank = "Class", remove_empty_ranks = c(TRUE, TRUE)))
    x <- agglomerateByRank(se, rank = "Class")
    rd1 <- rowData(x)[, 1:3]
    x <- agglomerateByRank(se, rank = "Class", remove_empty_ranks = TRUE)
    rd2 <- rowData(x)
    expect_equal(rd1, rd2)
    # Test that make_unique work
    uniq <- agglomerateByRank(se, rank = "Species")
    not_uniq <- agglomerateByRank(se, rank = "Species", make_unique = FALSE)
    expect_true( !any( duplicated(rownames(uniq)) ) )
    expect_true( any( duplicated(rownames(not_uniq)) ) )
    
    data(GlobalPatterns, package="mia")
    tse <- GlobalPatterns
    # Test that dimentionality is the same for merging object by agglomerateByRank
    # and agglomerateByVariable.
    expect_equal(length(unique(rowData(tse)[,"Family"])),
                 nrow(agglomerateByRank(tse, rank="Family", na.rm = FALSE)))
    expect_equal(length(unique(rowData(tse)[,"Family"])),
                 nrow(agglomerateByVariable(tse, f="Family", MARGIN = 1, agg.na.rm = FALSE)))
    
    # Test that dimentionality is the same when NA values are removed.
    expect_equal(length((unique(rowData(tse)[,"Family"]))[ !is.na(unique(rowData(tse)[,"Family"])) ]),
                 nrow(agglomerateByRank(tse, rank="Family", na.rm = TRUE)))
      expect_equal(length((unique(rowData(tse)[,"Family"]))[ !is.na(unique(rowData(tse)[,"Family"])) ]),
                 nrow(agglomerateByVariable(tse, f="Family", MARGIN = 1, agg.na.rm = TRUE)))
    
})
