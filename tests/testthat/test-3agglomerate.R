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
    actual_family <- actual <- agglomerateByVariable(xtse, by = "rows", group = tax_factors)
    expect_s4_class(actual,class(xtse))
    expect_equal(dim(actual),c(8,10))
    expect_equal(assays(actual)$mat[8,1],c(c_NA = 8))
    expect_equal(assays(actual)$mat[7,1],c(c_h = 16))
    tax_factors <- mia:::.get_tax_groups(xtse, col = 1)
    actual_phylum <- actual <- agglomerateByVariable(xtse, by = "rows", group = tax_factors)
    expect_s4_class(actual,class(xtse))
    expect_equal(dim(actual),c(3,10))
    expect_equal(assays(actual)$mat[1,1],c(a = 6))
    expect_equal(assays(actual)$mat[2,1],c(b = 36))
    expect_equal(assays(actual)$mat[3,1],c(c = 24))
    #
    expect_error(agglomerateByRank(xtse,"",empty.rm=FALSE),
                 "'rank' must be a non-empty single character value")
    expect_error(agglomerateByRank(xtse,"Family",empty.rm=""),
                 "'empty.rm' must be TRUE or FALSE")
    expect_error(
        agglomerateByRank(xtse,"Family",empty.rm=FALSE,update.tree=""),
        "'update.tree' must be TRUE or FALSE")
    xtse2 <- xtse
    rowData(xtse2) <- NULL
    expect_error(agglomerateByRank(xtse2,"Family",empty.rm=FALSE),
                 "taxonomyData needs to be populated")
    #
    actual <- agglomerateByRank(xtse,"Family",empty.rm=FALSE)
    expect_equivalent(rowData(actual),rowData(actual_family))
    actual <- agglomerateByRank(xtse,"Phylum",empty.rm=FALSE)
    expect_equivalent(rowData(actual),rowData(actual_phylum))
    #
    actual <- agglomerateByRank(xtse,"Family", ignore.taxonomy = FALSE, empty.rm = TRUE)
    expect_equal(dim(actual),c(6,10))
    expect_equal(rowData(actual)$Family,c("c","d","e","f","g","h"))
    actual <- agglomerateByRank(xtse,"Family", ignore.taxonomy = FALSE, empty.rm = FALSE)
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
                 nrow(agglomerateByRank(enterotype,"Genus", ignore.taxonomy = FALSE,
                 empty.rm = FALSE)))

    # agglomeration in all its forms
    data(GlobalPatterns, package="mia")
    se <- GlobalPatterns
    actual <- agglomerateByRank(se, rank = "Family",
        ignore.taxonomy = FALSE, empty.rm = FALSE, update.tree = FALSE)
    expect_equal(dim(actual),c(603,26))
    # Tree is dropped if it is not updated
    expect_null(rowTree(actual))
    # Families are not monophyletic, so the tree cannot be updated. It is
    # removed with a warning.
    expect_warning(
        actual <- agglomerateByRank(se, rank = "Family",
            ignore.taxonomy = FALSE, empty.rm = FALSE, update.tree = TRUE),
        "not monophyletic")
    expect_equal(dim(actual),c(603,26))
    expect_null(rowTree(actual))
    # Test that warning occurs when assay contian binary or negative values
    se1 <- transformAssay(se, method = "pa")
    se2 <- se1
    assay(se2, "pa")[1, 1] <- -1
    expect_warning(agglomerateByRank(se1, rank = "Phylum"))
    expect_warning(agglomerateByRank(se1, rank = "Order"))

    # Load data
    data(GlobalPatterns, package="mia")
    tse <- GlobalPatterns

    # Check that empty.rm works
    # Get all phyla
    all_phyla <- unique( rowData(tse)$Phylum )

    # When empty.rm = FALSE, then phyla should also include NA --> one extra row
    test0 <- agglomerateByVariable(tse, by = 1, group = "Phylum", empty.rm = FALSE)
    test1 <- agglomerateByRank(tse, rank = "Phylum", empty.rm = FALSE)

    # Test that dimentionality is the same for merging object by agglomerateByRank
    # and agglomerateByVariable.
    expect_equal(nrow(test0), length(all_phyla))
    expect_equal(nrow(test1), length(all_phyla))

    # When empty.rm = TRUE, there should be as many rows as there are non-NA phyla
    test0 <- agglomerateByVariable(tse, by = 1, group = "Phylum", empty.rm = TRUE)
    test1 <- agglomerateByRank(tse, rank = "Phylum", empty.rm = TRUE)

    # Test that dimensionality is the same when NA values are removed.
    expect_equal(nrow(test0), length( all_phyla[!is.na(all_phyla)] ))
    expect_equal(nrow(test1), length( all_phyla[!is.na(all_phyla)] ))

    # Check that the names of tree tips correspond to rownames
    expect_true( all(rowTree(test0)$tip %in% rownames(test0)) )
    expect_true( all(rowTree(test1)$tip %in% rownames(test1)) )

    # Check that there are more taxa when agglomeration is to "Species" level
    test0 <- agglomerateByVariable(tse, by = 1, group = "Species", empty.rm = FALSE)
    test1 <- agglomerateByRank(tse, rank = "Species", empty.rm = FALSE)
    expect_equal(nrow(test0), 945)
    expect_equal(nrow(test1), 2307)

    # Test that empty.ranks.rm work
    expect_error(
        agglomerateByRank(tse, rank = "Class", empty.ranks.rm = NULL))
    expect_error(
        agglomerateByRank(tse, rank = "Class", empty.ranks.rm = "NULL"))
    expect_error(
        agglomerateByRank(tse, rank = "Class", empty.ranks.rm = 1))
    expect_error(
        agglomerateByRank(
            tse, rank = "Class", empty.ranks.rm = c(TRUE, TRUE)))

    # Add a column to rowData(se) to test that only NA rank columns are removed
    # when empty.ranks.rm = TRUE
    rank <- "Class"
    rowData(tse)[["test"]] <- rep(NA, nrow(rowData(tse)))
    x <- agglomerateByRank(tse, rank = rank)
    rd1 <- rowData(x)
    x <- agglomerateByRank(tse, rank = rank, empty.ranks.rm = TRUE)
    rd2 <- rowData(x)
    cols <- taxonomyRanks(tse)[ seq_len(which(taxonomyRanks(tse) == "Class")) ]
    cols <- c(cols, "test")
    expect_equal(rd1[, cols], rd2[, cols])
    expect_true( ncol(rd1) > ncol(rd2) )
    # Test that make.unique work
    uniq <- agglomerateByRank(tse, rank = "Species", empty.rm = FALSE)
    not_uniq <- agglomerateByRank(
        tse, rank = "Species", make.unique = FALSE, empty.rm = FALSE)
    expect_true( !any( duplicated(rownames(uniq)) ) )
    expect_true( any( duplicated(rownames(not_uniq)) ) )

    ## START test agglomerateBySet ##

    # Generate random sets
    N_set <- 30L
    sets <- sample(
        c(TRUE, FALSE),
        size = nrow(tse) * N_set,
        prob = c(0.2, 0.8),
        replace = TRUE
    )
    # Add sets to rowData
    sets <- sets |> matrix(nrow = nrow(tse))
    # Add row and column names
    colnames(sets) <- paste0("set_", seq_len(ncol(sets)))
    # Add sets to rowData
    rowData(tse) <- cbind(rowData(tse), sets)
    # Extract set columns
    set_columns <- grep("set_", colnames(rowData(tse)), value = TRUE)
    # Add pseudocount assay
    tse <- transformAssay(tse, assay.type = "counts", method = "pseudocount")
    # Introduce NA to pseudocount assay
    assay(tse, "pseudocount")[1, 1] <- NA
    # Agglomerate based on sets
    tse_set <- agglomerateBySet(
        tse, by = 1, group = set_columns, na.rm = TRUE
    )
    # Compute reference for counts assay
    set_counts <- crossprod(sets, assay(tse, "counts"))
    # Simulate NA removal by replacing it with 0
    assay(tse, "pseudocount")[1, 1] <- 0
    # Compute reference for pseudocount assay
    set_pseudocount <- crossprod(sets, assay(tse, "pseudocount"))
    # Compare output of agglomerateBySet with reference assays
    expect_equal(assay(tse_set, "counts"), set_counts,
        check.attributes = FALSE)
    expect_equal(assay(tse_set, "pseudocount"), set_pseudocount,
                 check.attributes = FALSE)

    ## END test agglomerateBySet ##

    # Load data from miaTime package
    data("SilvermanAGutData", package = "miaTime")
    se <- SilvermanAGutData

    # checking reference consensus sequence generation
    actual <- agglomerateByRank(se, "Genus", update.refseq = FALSE)
    # There should be only one exact match for each sequence
    seqs_test <- as.character( referenceSeq(actual) )
    seqs_ref <- as.character( referenceSeq(se) )
    expect_true(all(vapply(
        seqs_test, function(seq) sum(seqs_ref %in% seq) == 1,
        FUN.VALUE = logical(1) )) )
    # Merging creates consensus sequences.
    th <- runif(1, 0, 1)
    actual <- agglomerateByRank(
        se, "Genus", update.refseq = TRUE, threshold = th)
    seqs_test <- referenceSeq(actual)
    # Get single taxon as reference. Merge those sequences and test that it
    # equals to one that is output of agglomerateByRank
    seqs_ref <- referenceSeq(se)
    feature <- sample(na.omit(rowData(se)[["Genus"]]), 1)
    seqs_ref <- seqs_ref[ rowData(se)[["Genus"]] %in% feature ]
    seqs_ref <- .merge_refseq(
        seqs_ref, factor(rep(feature, length(seqs_ref))),
        threshold = th)
    seqs_test <- seqs_test[ names(seqs_test) %in% feature ]
    expect_equal(seqs_test, seqs_ref)

    # checking reference consensus sequence generation using 'Genus:Alistipes'
    actual <- agglomerateByRank(se,"Genus", update.refseq = FALSE)
    expect_equal(as.character(referenceSeq(actual)[["Alistipes"]]),
                 paste0("TCAAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGTAGGCGGTTTGATAA",
                        "GTTAGAGGTGAAATCCCGGGGCTTAACTCCGGAACTGCCTCTAATACTGTTAG",
                        "ACTAGAGAGTAGTTGCGGTAGGCGGAATGTATGGTGTAGCGGTGAAATGCTTA",
                        "GAGATCATACAGAACACCGATTGCGAAGGCAGCTTACCAAACTATATCTGACG",
                        "TTGAGGCACGAAAGCGTGGGG"))
    actual <- agglomerateByRank(se,"Genus", update.refseq = TRUE)
    expect_equal(as.character(referenceSeq(actual)[["Alistipes"]]),
                 paste0("BCNMKCKTTVWYCKKMHTTMYTKKKYKTMMMKNKHDYKYMKDYKKNHNNNYMM",
                        "KHHNDNNKTKMMMDNBHNBKKCTYMMCHNBNDDDNKSSHBNNRWDMYKKBNND",
                        "NYTDRRKDVHNKNDRVGRNDRSBRRAWTBYNHRKKKWRSSRKKRAAWKSSKWR",
                        "RWDWTNDBRVRRAMHHCMRDKKSSRARGSSVSYYHNYBRRVHNDNNHYKRMVV",
                        "YKVRDNNNSRAARSBDKGGKK"))
})

test_that("Tree agglomeration requires monophyletic groups", {
    # Example from #867: (A, B) form one clade and (C, D, E) another
    tree <- ape::read.tree(text = "((A:1,B:1)n1:1,(C:1,D:1,E:1)n2:1)root;")
    tips <- match(c("A", "B", "C", "D", "E"), tree$tip.label)
    df <- .tree_to_table(tree)
    # Monophyletic groups
    expect_true(.check_monophyly(df, tips, c("x", "x", "y", "y", "y")))
    # Single-member groups and the whole tree as one group
    expect_true(.check_monophyly(df, tips, c("a", "b", "c", "d", "e")))
    expect_true(.check_monophyly(df, tips, rep("x", 5)))
    # D belongs to the group of A and B, but it is in the other clade. Neither
    # group is monophyletic: Firmicutes is in two clades, and the clade of
    # Bacteroidota (C, E) also includes D.
    groups <- c("Firmicutes", "Firmicutes", "Bacteroidota", "Firmicutes",
        "Bacteroidota")
    expect_false(.check_monophyly(df, tips, groups))
    expect_setequal(
        .get_non_monophyletic_groups(df, tips, groups),
        c("Firmicutes", "Bacteroidota"))
    # A single group that is split by another group. Only that group is
    # reported.
    expect_false(.check_monophyly(df, tips, c("x", "x", "y", "z", "y")))
    expect_identical(
        .get_non_monophyletic_groups(df, tips, c("x", "x", "y", "z", "y")),
        "y")
    expect_identical(
        .get_non_monophyletic_groups(df, tips, c("x", "x", "y", "y", "y")),
        character(0))
    # Tips that are not linked to any feature are ignored
    expect_true(.check_monophyly(df, tips[c(1, 3, 5)], c("x", "y", "y")))
    # Features can be linked to internal nodes
    clade_ab <- ape::getMRCA(tree, c("A", "B"))
    expect_true(.check_monophyly(df, c(clade_ab, tips[3]), c("x", "y")))
    expect_false(.check_monophyly(df, c(clade_ab, tips[1]), c("x", "y")))
    # Multiple features linked to the same node with different groups
    expect_false(.check_monophyly(df, tips[c(1, 1)], c("x", "y")))
    # Through the user-facing function: monophyletic groups are agglomerated
    tse <- TreeSummarizedExperiment(
        assays = list(counts = matrix(1:10, nrow = 5,
            dimnames = list(c("A", "B", "C", "D", "E"), c("s1", "s2")))),
        rowTree = tree)
    res <- agglomerateByVariable(
        tse, by = "rows", group = c("x", "x", "y", "y", "y"),
        update.tree = TRUE)
    expect_identical(rowLinks(res)$nodeLab, c("n1", "n2"))
    expect_identical(sort(rowTree(res)$tip.label), c("n1", "n2"))
    # Non-monophyletic groups: warning that names the groups, and the tree is
    # removed
    expect_warning(
        res <- agglomerateByVariable(
            tse, by = "rows", group = groups, update.tree = TRUE),
        "2 group\\(s\\) are not monophyletic: '(Firmicutes|Bacteroidota)'")
    expect_null(rowTree(res))
    expect_identical(nrow(res), 2L)
})

test_that(".get_levels_of_nodes", {
    # Branch lengths do not affect levels, also zero-length branches
    tree <- ape::read.tree(text = "((A:0,B:5)n1:0,(C:1,(D:1,E:1)n3:0)n2:1)root;")
    df <- .tree_to_table(tree)
    parent <- match(df[["parent"]], df[["node"]])
    level <- .get_levels_of_nodes(parent)
    labels <- df[["label"]]
    expect_identical(
        setNames(level, labels)[c("root", "n1", "n2", "n3", "A", "C", "D")],
        c(root = 0L, n1 = 1L, n2 = 1L, n3 = 2L, A = 2L, C = 2L, D = 3L))
    # Children always have a higher level than their parents
    not_root <- parent != seq_along(parent)
    expect_true(all(level[not_root] == level[ parent[not_root] ] + 1L))
    # Multiple roots or no root
    expect_error(.get_levels_of_nodes(c(1L, 2L, 1L, 2L)), "single root")
    expect_error(.get_levels_of_nodes(c(2L, 1L, 2L)), "single root")
    # Nodes that are not connected to the root (3 and 4 form a loop)
    expect_error(.get_levels_of_nodes(c(1L, 1L, 4L, 3L)), "connected")
})

test_that("Tree agglomeration with multiple trees, single group and tree name", {
    tree1 <- ape::read.tree(text = "((A:1,B:1)a1:1,(C:1,D:1)a2:1)a0;")
    tree2 <- ape::read.tree(text = "((E:1,F:1)b1:1,G:1)b0;")
    .create_tse <- function(tree){
        n <- length(tree$tip.label)
        mat <- matrix(seq_len(2*n), ncol = 2,
            dimnames = list(tree$tip.label, c("s1", "s2")))
        TreeSummarizedExperiment(assays = list(counts = mat), rowTree = tree)
    }
    # Tree agglomeration supports only a single tree. With multiple trees, the
    # trees are removed with a warning, and the data is merged as usual.
    tse <- rbind(.create_tse(tree1), .create_tse(tree2))
    group <- c("x", "x", "y", "y", "z", "z", "G")
    expect_warning(
        res <- agglomerateByVariable(
            tse, by = "rows", group = group, update.tree = TRUE),
        "multiple trees")
    ref <- agglomerateByVariable(
        tse, by = "rows", group = group, update.tree = FALSE)
    expect_null(rowTree(res))
    expect_identical(assay(res), assay(ref))
    # The name of the tree is kept
    tse <- .create_tse(tree1)
    rowTreeNames(tse) <- "my_tree"
    res <- agglomerateByVariable(
        tse, by = "rows", group = c("x", "x", "y", "y"), update.tree = TRUE)
    expect_identical(rowTreeNames(res), "my_tree")
    expect_identical(unname(rowLinks(res)$nodeLab), c("a1", "a2"))
    # Single group: its MRCA is the root, so the tree has a single tip
    res <- agglomerateByVariable(
        .create_tse(tree1), by = "rows", group = rep("x", 4),
        update.tree = TRUE)
    expect_identical(rowTree(res)$tip.label, "a0")
    expect_identical(unname(rowLinks(res)$nodeLab), "a0")
    # Each feature is its own group, so the data is not merged. Rows are still
    # linked correctly.
    res <- agglomerateByVariable(
        .create_tse(tree1), by = "rows", group = c("p", "q", "r", "s"),
        update.tree = TRUE)
    expect_identical(unname(rowLinks(res)$nodeLab), c("A", "B", "C", "D"))
})
