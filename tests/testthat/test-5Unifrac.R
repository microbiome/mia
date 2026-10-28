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

    # Test the function with agglomerated data. .get_unifrac renames
    # rownames based on tips and links to them. Then it also prunes the tree so
    # that rows are in tips.
    tse <- GlobalPatterns
    tse <- agglomerateByRank(tse, rank = "Species")
    tse_ref <- tse
    rownames(tse_ref) <- rowLinks(tse_ref)[["nodeLab"]]
    rowTree(tse_ref) <- .prune_tree(rowTree(tse_ref), rowLinks(tse_ref)[["nodeLab"]])
    # Calculate unweighted unifrac
    unifrac_mia <- as.matrix(getDissimilarity(tse, method = "unifrac",
                                              weighted = FALSE))
    unifrac_ecodive <- as.matrix(ecodive::unweighted_unifrac(t(assay(tse_ref)),
                                              rowTree(tse_ref)))
    # Calculate weighted unifrac. No tolerance needed since the tree has
    # simpler structure after pruning.
    unifrac_mia <- as.matrix(getDissimilarity(tse, method = "unifrac",
                                              weighted = TRUE))
    unifrac_ecodive <- as.matrix(ecodive::weighted_unifrac(t(assay(tse_ref)),
                                              rowTree(tse_ref)))
    expect_equal(unifrac_mia, unifrac_ecodive)
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

test_that("Normalized weighted Unifrac", {
    data(GlobalPatterns, package="mia")
    tse <- GlobalPatterns
    expect_error(
        getDissimilarity(tse, method = "unifrac", weighted = TRUE,
                         normalized = 1)
    )
    expect_error(
        getDissimilarity(tse, method = "unifrac", weighted = TRUE,
                         normalized = c(TRUE, FALSE))
    )
    # normalized is ignored for unweighted Unifrac
    expect_equal(
        getDissimilarity(tse, method = "unifrac", weighted = FALSE,
                         normalized = TRUE),
        getDissimilarity(tse, method = "unifrac", weighted = FALSE)
    )
    # Compare against the definition on a small tree:
    # sum_b len_b * |pA_b - pB_b| / sum_tips d_tip * (pA_tip + pB_tip)
    tree <- ape::read.tree(text = "((a:1,b:2):0.5,(c:1.5,d:0.5):1);")
    mat <- matrix(
        c(10, 0, 5, 1,
          0, 3, 2, 8,
          4, 4, 0, 0),
        nrow = 4, dimnames = list(c("a", "b", "c", "d"), c("s1", "s2", "s3")))
    n_tip <- ape::Ntip(tree)
    .get_tips <- function(node){
        if( node <= n_tip ) tree$tip.label[node] else
            ape::extract.clade(tree, node)$tip.label
    }
    depth <- ape::node.depth.edgelength(tree)[seq_len(n_tip)]
    names(depth) <- tree$tip.label
    ref <- function(x, y){
        p <- x / sum(x)
        q <- y / sum(y)
        num <- sum(vapply(seq_len(nrow(tree$edge)), function(i){
            tips <- .get_tips(tree$edge[i, 2])
            tree$edge.length[i] * abs(sum(p[tips]) - sum(q[tips]))
        }, numeric(1)))
        den <- sum(depth[names(p)] * (p + q))
        num / den
    }
    res <- as.matrix(mia:::.get_unifrac(
        t(mat), tree, weighted = TRUE, normalized = TRUE))
    expect_equal(res["s1", "s2"], ref(mat[, "s1"], mat[, "s2"]))
    expect_equal(res["s1", "s3"], ref(mat[, "s1"], mat[, "s3"]))
    expect_equal(res["s2", "s3"], ref(mat[, "s2"], mat[, "s3"]))
    expect_true(all(res >= 0 & res <= 1))

    # Compare against ecodive
    skip_if_not(requireNamespace("ecodive", quietly = TRUE))
    unifrac_mia <- as.matrix(getDissimilarity(tse, method = "unifrac",
                                              weighted = TRUE,
                                              normalized = TRUE))
    unifrac_ecodive <- as.matrix(ecodive::normalized_unifrac(t(assay(tse)),
                                              rowTree(tse)))
    expect_equal(unifrac_mia, unifrac_ecodive)
    expect_true(all(unifrac_mia >= 0 & unifrac_mia <= 1))
})
