
test_that("applyByModule", {
    # Create sample experiment
    set.seed(123)
    tse <- TreeSummarizedExperiment::makeTSE()
    # Set assay name and alpha indices
    assayNames(tse) <- "counts"
    alpha.index <- c("observed", "shannon")
    # Check group depletion
    expect_error(
        applyByModule(
            tse, "rows", "var1", addAlpha, index = alpha.index, min.group.size = 2
        ),
        "No '.group' has size above 'min.group.size'.",
        fixed = TRUE
    )
    # Add numeric variable to rowData
    rowData(tse)$var1 <- rep(letters[1:2], 5)
    # Run getAlpha
    tse <- applyByModule(
        tse, "rows", "var1", getAlpha, index = alpha.index, min.group.size = 2
    )
    # Check alpha results
    out.cols <- interaction(c("observed", "shannon"), letters[1:2], sep = "_")
    expect_contains(names(colData(tse)), out.cols)
    # Run addAlpha
    tse <- applyByModule(
        tse, "rows", "var1", addAlpha, index = alpha.index, min.group.size = 2
    )
    # Check embedded results
    expect_in(altExpNames(tse), letters[1:2])
    # Define set names
    set.vars <- c("set1", "set2")
    # Add set vars to rowData
    rowData(tse)[set.vars] <- matrix(
        sample(c(0, 1), 2 * nrow(tse), replace = TRUE), ncol = 2
    )
    # Run getDominant set-wise
    tse <- applyByModule(
        tse, "rows", set.vars, getDominant,
        group = "var1"
    )
    # Check dominance results
    expect_contains(names(colData(tse)), paste("dominant", set.vars, sep = "_"))
    # Run plotRowTree set-wise
    tse <- applyByModule(
        tse, "rows", set.vars,  miaViz::plotRowTree, meta.name = "tree_plots"
    )
    # Check tree results
    expect_in(names(metadata(tse)), "tree_plots")
    # applyByModule works with workers in separate sessions
    old <- BiocParallel::bpparam()
    BiocParallel::register(BiocParallel::SnowParam(2))
    on.exit(BiocParallel::register(old))
    # Workers do not have bluster attached, so KmeansParam must be evaluated
    # before the arguments are sent to them
    tse <- applyByModule(
        tse, "rows", set.vars, getCluster,
        BLUSPARAM = bluster::KmeansParam(centers = 3),
        assay.type = "counts", by = "cols", clust.col = "kmeans"
    )
    # Check clustering results
    expect_contains(names(colData(tse)), paste("kmeans", set.vars, sep = "_"))
})
