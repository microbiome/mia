
test_that("exporters", {
    
    tse <- TreeSummarizedExperiment::makeTSE()
    assayNames(tse) <- "counts"
    
    # Test exportRaw
    dpath <- tempfile()
    exportRaw(tse, dpath)
    
    expect_true(dir.exists(dpath))
    
    row_data <- read.table(file.path(dpath, "rowdata.tsv"))
    col_data <- read.table(file.path(dpath, "coldata.tsv"))
    assay_tab <- read.table(file.path(dpath, "assays", "counts.tsv"))
    
    row_tree <- ape::read.tree(file.path(dpath, "row_trees", "phylo.nwk"))
    col_tree <- ape::read.tree(file.path(dpath, "col_trees", "phylo.nwk"))
    
    tse2 <- TreeSummarizedExperiment::TreeSummarizedExperiment(
        assays = S4Vectors::SimpleList(counts = as.matrix(assay_tab)),
        rowData = row_data,
        colData = col_data,
        rowTree = row_tree,
        colTree = col_tree
    )
    
    expect_equal(tse, tse2)
    
    # Test exportQIIME2
    dpath <- tempfile()
    exportQIIME2(tse, dpath)
    
    expect_true(dir.exists(dpath))
    
    assay_tab <- read.table(
        file.path(dpath, "counts.tsv"),
        sep = "\t",
        header = TRUE,
        comment.char = "",
        check.names = FALSE
    )
    
    expect_equal(colnames(assay_tab)[[1L]], "#OTU ID")
    
    row_data <- read.table(
        file.path(dpath, "taxonomy.tsv"),
        sep = "\t",
        header = TRUE,
        comment.char = "",
        check.names = FALSE
    )
    
    expect_named(row_data, c("Feature ID", "Taxon", "Confidence"))
    
    col_data <- read.table(
        file.path(dpath, "metadata.tsv"),
        sep = "\t",
        header = TRUE,
        comment.char = "",
        check.names = FALSE
    )
    
    expect_named(col_data, c("sample-id", "ID", "group"))
    expect_equal(col_data$`sample-id`[1L], "#q2:types")
    
    # Test exportMothur
    dpath <- tempfile()
    exportMothur(tse, dpath)
    
    expect_true(dir.exists(dpath))
    
    row_data <- read.table(
        file.path(dpath, "taxonomy.tsv"),
        sep = "\t",
        header = TRUE,
        comment.char = "",
        check.names = FALSE
    )
    
    expect_named(row_data, c("OTU", "Size", "Taxonomy"))
    
    expect_type(row_data$OTU, "character")
    expect_type(row_data$Size, "integer")
    # expect_type(row_data$Taxonomy, "character")
    
    col_data <- read.table(
        file.path(dpath, "metadata.tsv"),
        sep = "\t",
        header = TRUE,
        comment.char = "",
        check.names = FALSE
    )
    
    # exportMothur should give warning that group varname is taken
    expect_named(col_data, c("group", "ID", "group.1"))
})
