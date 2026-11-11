context("splitOn")
test_that("splitOn", {
    data(GlobalPatterns, package="mia")
    x <- GlobalPatterns

    ################################## splitOn #################################
    # Test that throughs an error
    expect_error(splitOn(x, "test"))
    expect_error(splitOn(x[1:10, 1:10], x$SampleType))
    set.seed(103)
    rowData(x)$group <- sample(1:10, nrow(x), replace = TRUE)
    set.seed(374)
    colData(x)$group <- sample(1:10, ncol(x), replace = TRUE)
    expect_error(splitOn(x, "group"))
    expect_error(splitOn(x, "SampleType", by = 1))
    expect_error(splitOn(x, "Phylum", by = 2))
    expect_error(splitOn(x, x$SampleType, by = 1))
    expect_error(splitOn(x, rowData(x)$Phylum, by = 2))
    expect_error(splitOn(x))
    expect_error(splitOn(assay(x), x$SampleType))
    expect_error(splitOn(x, "SampleType", use.names = 1))
    expect_error(splitOn(x, "SampleType", use.names = "TRUE"))

    # Test that names of elemetns are correct
    list <- splitOn(x, "SampleType")
    expect_equal(names(list), as.character(unique(x$SampleType)) )
    list <- splitOn(x, "SampleType", use.names = FALSE)
    expect_equal(names(list), NULL )

    # Test that col-wie split is done correctly
    list <- splitOn(x, "group", by = 1)
    expect_equal( colnames(list[[1]]), colnames(x) )
    expect_true( length(list) == 10 )

    # Test that row-wise split is done correctly
    list <- splitOn(x, "group", by = "cols")
    expect_equal( rownames(list[[1]]), rownames(x) )
    expect_true( length(list) == 10 )

    ################################# unsplitOn ################################
    # Test that error occurs
    expect_error( unsplitOn(x) )
    mod_list <- list
    mod_list[[1]] <- mod_list[[1]][1:2, 1:2]
    expect_error( unsplitOn(mod_list) )

    # Test that works
    x_sub <- x[1:100, 1:10]
    # Split
    list <- splitOn(x_sub, "group", by = "features")
    # Unsplit
    unsplitted <- unsplitOn(list)
    # Order the data
    unsplitted <- unsplitted[ rownames(x_sub), colnames(x_sub) ]
    # Convert delayed matrix to normal
    assay(unsplitted) <- as.matrix( assay(unsplitted) )
    expect_equal(assay(x_sub), assay(unsplitted) )

    x_sub <- x[1:100, 1:10]
    # Split
    list <- splitOn(x_sub, "group", by = 2)
    # Unsplit
    unsplitted <- unsplitOn(list)
    # Order the data
    unsplitted <- unsplitted[ rownames(x_sub), colnames(x_sub) ]
    # Convert delayed matrix to normal
    assay(unsplitted) <- as.matrix( assay(unsplitted) )
    expect_equal(assay(x_sub), assay(unsplitted) )

    list <- splitOn(x, "SampleType")
    unsplitted <- unsplitOn(list)
    # Order the data
    unsplitted <- unsplitted[ rownames(x), colnames(x) ]
    # Convert delayed matrix to normal
    assay(unsplitted) <- as.matrix( assay(unsplitted) )
    expect_equal(assay(x), assay(unsplitted) )

    # Split
    list <- splitOn(x, "Phylum")
    # Unsplit
    unsplitted <- unsplitOn(list)
    # Order the data
    unsplitted <- unsplitted[ rownames(x), colnames(x) ]
    # Convert delayed matrix to normal
    assay(unsplitted) <- as.matrix( assay(unsplitted) )
    expect_equal(assay(x), assay(unsplitted) )
})

test_that("splitOn and unsplitOn do not alter trees", {
    data(GlobalPatterns, package = "mia")
    x <- GlobalPatterns[1:100, 1:10]
    set.seed(103)
    rowData(x)$group <- sample(1:3, nrow(x), replace = TRUE)
    # Split by rows: each element keeps the tree and links unchanged
    list <- splitOn(x, "group", by = 1)
    for( se in as.list(list) ){
        expect_identical(rowTree(se), rowTree(x))
        expect_identical(rowLinks(se), rowLinks(x)[match(rownames(se), rownames(x)), ])
    }
    # Unsplit: the tree and links are the same as in the original data
    res <- unsplitOn(list)
    expect_identical(rowTree(res), rowTree(x))
    expect_identical(rowLinks(res), rowLinks(x)[match(rownames(res), rownames(x)), ])
    # Split by columns: rows and their links are unchanged
    list <- splitOn(x, "SampleType")
    for( se in as.list(list) ){
        expect_identical(rowTree(se), rowTree(x))
        expect_identical(rowLinks(se), rowLinks(x))
    }
    res <- unsplitOn(list)
    expect_identical(rowTree(res), rowTree(x))
    expect_identical(rowLinks(res), rowLinks(x))
})

test_that("unsplitOn keeps multiple trees", {
    t1 <- ape::read.tree(text = "((A:1,B:1)a1:1,(C:1,D:1)a2:1)a0;")
    t2 <- ape::read.tree(text = "((E:1,F:1)b1:1,G:1)b0;")
    .create_tse <- function(tree){
        n <- length(tree$tip.label)
        mat <- matrix(seq_len(2*n), ncol = 2,
            dimnames = list(tree$tip.label, c("s1", "s2")))
        TreeSummarizedExperiment(assays = list(counts = mat), rowTree = tree)
    }
    x <- rbind(.create_tse(t1), .create_tse(t2))
    # The first element has rows only from the first tree, so subsetting drops
    # the second tree from it
    list <- splitOn(x, c("p", "p", "p", "q", "q", "q", "q"), by = 1)
    expect_identical(rowTreeNames(list[[1]]), "phylo")
    res <- unsplitOn(list)
    expect_identical(rowTreeNames(res), rowTreeNames(x))
    expect_identical(rowTree(res, "phylo"), t1)
    expect_identical(rowTree(res, "phylo.1"), t2)
    expect_identical(
        rowLinks(res), rowLinks(x)[match(rownames(res), rownames(x)), ])
    # The result is valid, e.g., it can be subsetted
    expect_identical(rowTreeNames(res[5:7, ]), "phylo.1")
})
