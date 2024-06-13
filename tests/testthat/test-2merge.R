context("merge")
test_that("merge", {
    # .check_f
    expect_error(mia:::.norm_f(),
                 'argument "f" is missing')
    expect_error(mia:::.norm_f(6),
                 'argument "f" is missing')
    expect_error(mia:::.norm_f(6,5),
                 "'f' must be a factor or character vector")
    f <- factor(c(rep("a",3),rep("b",3)))
    expect_true(is.factor(mia:::.norm_f(6,f)))
    # .check_archetype
    expect_error(mia:::.norm_archetype(),
                 'argument "archetype" is missing')
    expect_error(mia:::.norm_archetype(f),
                 'argument "archetype" is missing')
    expect_error(mia:::.norm_archetype(f),
                 'argument "archetype" is missing')
    expect_equal(mia:::.norm_archetype(f, 1),c(1,1))
    expect_equal(mia:::.norm_archetype(f, c(1,2)),c(1,2))
    expect_error(mia:::.norm_archetype(f, c(1,2,3)),
                 "length of 'archetype' must have the same length as levels")
    expect_error(mia:::.norm_archetype(f, c(5)),
                 "'archetype' out of bounds for some levels of 'f'")
    # .norm_archetype
    expect_error(mia:::.norm_archetype(),
                 'argument "archetype" is missing')
    expect_error(mia:::.norm_archetype(f),
                 'argument "archetype" is missing')
    actual <- mia:::.norm_archetype(f, c(1,2))
    expect_equal(actual, c(1,2))
    actual <- mia:::.norm_archetype(f, c(1))
    expect_equal(actual, c(1,1))

    # .get_element_pos
    expect_error(mia:::.get_element_pos(),
                 'argument "archetype" is missing')
    expect_error(mia:::.get_element_pos(f),
                 'argument "archetype" is missing')

    actual <- mia:::.get_element_pos(f, archetype = mia:::.norm_archetype(f, 1))
    expect_equal(actual,c(a = 1, b = 4))
    actual <- mia:::.get_element_pos(f, archetype = mia:::.norm_archetype(f, 2))
    expect_equal(actual,c(a = 2, b = 5))
    actual <- mia:::.get_element_pos(f, archetype = c(2,1))
    expect_equal(actual,c(a = 2, b = 4))

    # .merge_rows
    mat <- matrix(1:60, nrow = 6)
    gr <- GRanges("chr1",rep("1-6",6))
    df <- DataFrame(n = c(1:6))
    mcols(gr) <- df
    grl <- splitAsList(gr,1:6)
    expect_error(mia:::.merge_rows(),
                 'argument "f" is missing')
    x <- SummarizedExperiment(assays = list(mat = mat))
    xr <- SummarizedExperiment(assays = list(mat = mat),
                               rowRanges = gr)
    xrl <- SummarizedExperiment(assays = list(mat = mat),
                                rowRanges = unname(grl))
    expect_error(mia:::.merge_rows(x),
                 'argument "f" is missing')
    FUN_check_x <- function(x,archetype=1){
        actual <- agglomerateByVariable(x, MARGIN = "rows", f, archetype)
        expect_s4_class(actual,class(x))
        expect_equal(dim(actual),c(2,10))
    }
    lapply(list(x,xr,xrl),FUN_check_x)
    lapply(list(x,xr,xrl),FUN_check_x,archetype=2)
    #
    f <- factor(c(rep("a",3),rep("b",3)))
    mat <- matrix(1:60, nrow = 6)
    gr <- GRanges("chr1",rep("1-6",6))
    df <- DataFrame(n = c(1:6))
    mcols(gr) <- df
    grl <- splitAsList(gr,1:6)
    xtse <- TreeSummarizedExperiment(assays = list(mat = mat),
                                     rowRanges = unname(grl))
    FUN_check_x <- function(x,archetype=1){
        actual <- agglomerateByVariable(x, MARGIN = "rows", f, archetype, 
            agglomerate.tree = FALSE)
        expect_s4_class(actual,class(x))
        expect_equal(dim(actual),c(2,10))
    }
    lapply(list(xtse),FUN_check_x)
    lapply(list(xtse),FUN_check_x,archetype=2)
    # Check multiple rowTrees
    data(esophagus, package="mia")
    data(GlobalPatterns, package="mia")
    # Add arbitrary groups
    rowData(esophagus)$group <- c(rep(c("A", "B", "C"), each = nrow(esophagus)/3),
                                  rep("A", nrow(esophagus)-round(nrow(esophagus)/3)*3) )
    rowData(esophagus)$group2 <- c(rep(c("A", "B", "C"), each = nrow(esophagus)/3),
                                   rep("A", nrow(esophagus)-round(nrow(esophagus)/3)*3) )
    rowData(GlobalPatterns)$group <- c(rep(c("C", "D", "E"), each = nrow(GlobalPatterns)/3),
                                       rep("C", nrow(GlobalPatterns)-round(nrow(GlobalPatterns)/3)*3) )
    # Merge
    tse <- mergeSEs(esophagus, GlobalPatterns, assay.type="counts")
    # Reorder data since mergeSEs does not order the data based on original order
    # (trees are pruned differently --> first instance represent specific branch)
    tse <- tse[c(rownames(esophagus), rownames(GlobalPatterns)), ]
    # Only esophagus has these groups --> the merge should contain only esophagus
    merged  <- agglomerateByVariable(tse, MARGIN = "rows",
                                    f = rowData(tse)$group2, agglomerate.tree=TRUE)
    merged2 <- agglomerateByVariable(tse, MARGIN = "rows",
                                    f = rowData(tse)$group2, agglomerate.tree = FALSE)
    merged3 <- agglomerateByVariable(esophagus, MARGIN = "rows",
                                    f = rowData(esophagus)$group2,
                                    agglomerate.tree = TRUE)
    merged4 <- .merge_features(tse, merge.by = rowData(tse)$group2,
                                    agglomerate.tree = TRUE)
    merged5 <- agglomerateByVariable(tse, MARGIN = "rows",
                                    f = rowData(tse)$group2, agglomerate.tree = TRUE)
    expect_equal( rowLinks(merged)$whichTree,
                  rowLinks(merged2)$whichTree )
    expect_false( all(rowLinks(merged) == rowLinks(merged2)) )
    expect_equal(rowTree(tse), rowTree(merged2))
    expect_equal(rowTree(merged), rowTree(merged3))
    expect_equal(merged4, merged5)
    expect_equal(agglomerateByVariable(tse, MARGIN = "rows",
                                        f=rowData(tse)$group2),
                agglomerateByVariable(tse, MARGIN = "rows",
                                        f=rowData(tse)$group2))

    # Both datasets have group variable
    merged <- agglomerateByVariable(tse, MARGIN = "rows",
                                    f = rowData(tse)$group, agglomerate.tree = TRUE)
    merged2 <- agglomerateByVariable(tse, MARGIN = "rows",
                                    f = rowData(tse)$group, agglomerate.tree = FALSE)
    expect_equal( rowLinks(merged)$whichTree,
                  rowLinks(merged2)$whichTree )
    expect_false( all(rowLinks(merged) == rowLinks(merged2)) )
    expect_true( rowTree(merged, "phylo")$Nnode < rowTree(merged2, "phylo")$Nnode )
    expect_true( rowTree(merged, "phylo.1")$Nnode < rowTree(merged2, "phylo.1")$Nnode )
})
