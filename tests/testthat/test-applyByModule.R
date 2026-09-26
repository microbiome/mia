
test_that("applyByModule", {
    
    set.seed(123)
    tse <- TreeSummarizedExperiment::makeTSE()
    
    assayNames(tse) <- "counts"
    alpha.index <- c("observed", "shannon")
    
    expect_error(
        applyByModule(
            tse, "rows", "var1", addAlpha, index = alpha.index, min.group.size = 2
        ),
        "No '.group' has size above 'min.group.size'.",
        fixed = TRUE
    )
    
    rowData(tse)$var1 <- rep(letters[1:2], 5)
    
    tse <- applyByModule(
        tse, "rows", "var1", getAlpha, index = alpha.index, min.group.size = 2
    )
    
    out.cols <- interaction(c("observed", "shannon"), letters[1:2])
    expect_contains(names(colData(tse)), out.cols)
    
    tse <- applyByModule(
        tse, "rows", "var1", addAlpha, index = alpha.index, min.group.size = 2
    )
    
    expect_in(altExpNames(tse), letters[1:2])
    
    mod.vars <- c("mod1", "mod2")
    
    rowData(tse)[mod.vars] <- matrix(
        sample(c(0, 1), 2 * nrow(tse), replace = TRUE), ncol = 2
    )
    
    tse <- applyByModule(
        tse, "rows", mod.vars, getCluster,
        BLUSPARAM = bluster::KmeansParam(centers = 3),
        assay.type = "counts", by = "cols", clust.col = "kmeans"
    )
    
    expect_contains(names(colData(tse)), paste("kmeans", mod.vars, sep = "."))
    
    tse <- applyByModule(
        tse, "rows", mod.vars, getDominant,
        group = "var1"
    )
    
    expect_contains(names(colData(tse)), paste("dominant", mod.vars, sep = "."))
    
    tse <- applyByModule(
        tse, "rows", mod.vars,  miaViz::plotRowTree, meta.name = "tree_plots"
    )
    
    expect_in(names(metadata(tse)), "tree_plots")
})
