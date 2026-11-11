context("Unifrac beta diversity")
test_that("Unifrac beta diversity", {
    data(esophagus, package="mia")
    tse <- esophagus
    tse <- transformAssay(tse, assay.type="counts", method="relabundance")

    expect_error(
        getDissimilarity(tse, method = "unifrac",assay.type = "test",
                          tree.name = "phylo", weighted = FALSE)
    )
    expect_error(
        getDissimilarity(tse, method = "unifrac", assay.type = 2,
                          tree.name = "phylo", weighted = FALSE)
    )
    expect_error(
        getDissimilarity(tse, method = "unifrac", assay.type = TRUE,
                          tree.name = "phylo", weighted = FALSE)
    )
    expect_error(
      getDissimilarity(tse, method = "unifrac", assay.type = "counts",
                       tree.name = "test", weighted = FALSE)
    )
    expect_error(
      getDissimilarity(tse, method = "unifrac", assay.type = "counts",
                       tree.name = 1, weighted = FALSE)
    )
    expect_error(
      getDissimilarity(tse, method = "unifrac", assay.type = "counts",
                       tree.name = TRUE, weighted = FALSE)
    )
    expect_error(
      getDissimilarity(tse, method = "unifrac", assay.type = "counts",
                       tree.name = "phylo", weighted = 1)
    )

    data(GlobalPatterns, package="mia")
    tse <- GlobalPatterns
    # Calculate unweighted unifrac
    unifrac_mia <- as.matrix(getDissimilarity(tse, method = "unifrac",
                                              weighted = FALSE))
    unifrac_ecodive <- as.matrix(ecodive::unweighted_unifrac(t(assay(tse)),
                                              rowTree(tse)))
    expect_equal(unifrac_mia, unifrac_ecodive)
    # Calculate weighted unifrac. Allow tolerance since weighted unifrac
    # calculation may have some stochasticity.
    unifrac_mia <- as.matrix(getDissimilarity(tse, method = "unifrac",
                                              weighted = TRUE))
    unifrac_ecodive <- as.matrix(ecodive::weighted_unifrac(t(assay(tse)),
                                              rowTree(tse)))
    expect_equal(unifrac_mia, unifrac_ecodive, tolerance = 5e-3)

    # Test that the function works correctly when there are multiple trees.
    # The function should subset the data based on tree.
    tse <- GlobalPatterns
    trees <- list(phylo = rowTree(tse), tree2 = rowTree(tse))
    links <- rowLinks(tse)
    links[ 1:500 , "whichTree"] <- "tree2"
    tse@rowTree <- trees
    tse@rowLinks <- links
    tse_ref <- tse
    tse_ref <- tse_ref[ rowLinks(tse_ref)[["whichTree"]] == "tree2", ]
    # Calculate unweighted unifrac
    expect_warning(
    unifrac_mia <- getDissimilarity(tse, method = "unifrac",
                                    weighted = FALSE, tree.name = "tree2")
    )
    unifrac_mia <- as.matrix(unifrac_mia)

    skip_if_not(requireNamespace("ecodive", quietly = TRUE))
    unifrac_ecodive <- as.matrix(ecodive::unweighted_unifrac(t(assay(tse_ref)),
                                              rowTree(tse_ref)))
    expect_equal(unifrac_mia, unifrac_ecodive)
    # Calculate weighted unifrac
    expect_warning(
    unifrac_mia <- getDissimilarity(tse, method = "unifrac",
                                    weighted = TRUE, tree.name = "tree2")
    )
    unifrac_mia <- as.matrix(unifrac_mia)
    unifrac_ecodive <- as.matrix(ecodive::weighted_unifrac(t(assay(tse_ref)),
                                              rowTree(tse_ref)))
    expect_equal(unifrac_mia, unifrac_ecodive, tolerance = 1e-3)

    # Test that rows are linked to tips based on rowLinks, and that the tree is
    # pruned so that it includes only tips that are linked to rows.
    tree <- ape::read.tree(
        text = "((A:1,B:2)n1:1,(C:1,D:3)n2:2,(E:1,F:1)n3:1)root;")
    mat <- matrix(
        c(10, 5, 3, 1, 4, 8, 7, 6, 0), nrow = 3,
        dimnames = list(c("x", "y", "z"), c("s1", "s2", "s3")))
    tse <- TreeSummarizedExperiment(
        assays = list(counts = mat), rowTree = tree,
        rowNodeLab = c("A", "C", "D"))
    # Reference: rows are named by their tips, and B, E and F are dropped from
    # the tree
    tree_ref <- ape::read.tree(text = "(A:2,(C:1,D:3)n2:2)root;")
    mat_ref <- mat
    rownames(mat_ref) <- c("A", "C", "D")
    # Calculate unweighted unifrac
    expect_warning(
        unifrac_mia <- getDissimilarity(
            tse, method = "unifrac", weighted = FALSE),
        "Pruning tree")
    unifrac_ecodive <- ecodive::unweighted_unifrac(t(mat_ref), tree_ref)
    expect_equal(as.matrix(unifrac_mia), as.matrix(unifrac_ecodive))
    # Calculate weighted unifrac
    expect_warning(
        unifrac_mia <- getDissimilarity(
            tse, method = "unifrac", weighted = TRUE),
        "Pruning tree")
    unifrac_ecodive <- ecodive::weighted_unifrac(t(mat_ref), tree_ref)
    expect_equal(as.matrix(unifrac_mia), as.matrix(unifrac_ecodive))
})

test_that("Unifrac does not depend on the scale of abundances", {
    # S1 holds A and C in proportions 1/3 and 2/3, S2 holds B and C in
    # proportions 2/3 and 1/3, and every branch has length 1. Weighted Unifrac
    # sums |p1 - p2| over branches: A 1/3 + B 2/3 + (A,B) 1/3 + C 1/3 = 5/3.
    # Unweighted Unifrac is 2/4: A and B are unique, (A,B) and C are shared.
    tree <- ape::read.tree(text = "((A:1,B:1):1,C:1);")
    mat <- cbind(S1 = c(A = 1, B = 0, C = 2), S2 = c(A = 0, B = 2, C = 1))
    tse <- TreeSummarizedExperiment(assays = list(counts = mat),
                                    rowTree = tree)
    tse <- transformAssay(tse, method = "relabundance")
    # Percentages whose integer parts do not sum to 100
    assay(tse, "percent") <- assay(tse, "relabundance") * 100
    for( assay.type in assayNames(tse) ){
        weighted <- getDissimilarity(tse, method = "unifrac", weighted = TRUE,
                                     assay.type = assay.type)
        unweighted <- getDissimilarity(tse, method = "unifrac",
                                       weighted = FALSE,
                                       assay.type = assay.type)
        expect_equal(as.numeric(weighted), 5/3, info = assay.type)
        expect_equal(as.numeric(unweighted), 1/2, info = assay.type)
    }
})

test_that("Unifrac merges rows that are linked to the same node", {
    skip_if_not(requireNamespace("ecodive", quietly = TRUE))
    tree <- ape::read.tree(text = "((A:1,B:2)n1:1,(C:1,D:3)n2:2)root;")
    mat <- matrix(
        c(10, 0, 5, 3, 1, 0, 4, 8, 7, 2, 6, 0), nrow = 4,
        dimnames = list(c("x", "y", "z", "w"), c("s1", "s2", "s3")))
    # Rows x and y are linked to the same tip A
    node_lab <- c("A", "A", "C", "D")
    # Reference: rows linked to the same tip are summed
    mat_ref <- rbind(
        A = mat["x", ] + mat["y", ], C = mat["z", ], D = mat["w", ])
    tree_ref <- ape::read.tree(text = "((A:1)n1:1,(C:1,D:3)n2:2)root;")
    tree_ref <- ape::collapse.singles(tree_ref)
    for( weighted in c(FALSE, TRUE) ){
        res <- .get_unifrac(
            t(mat), tree, weighted = weighted, node.label = node_lab) |>
            suppressWarnings()
        ref <- if( weighted ){
            ecodive::weighted_unifrac(t(mat_ref), tree_ref)
        } else{
            ecodive::unweighted_unifrac(t(mat_ref), tree_ref)
        }
        expect_equal(as.matrix(res), as.matrix(ref))
    }
})
