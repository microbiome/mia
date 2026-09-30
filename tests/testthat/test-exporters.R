
test_that("exporters", {
    
    set.seed(123)
    tse <- TreeSummarizedExperiment::makeTSE()
    assayNames(tse) <- "counts"
    
    taxranks <- c("Family", "Genus", "Species")
    rowData(tse)[taxranks] <- lapply(
        taxranks, function(x) sample(letters, nrow(tse), replace = TRUE)
    )
    
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
    
    expect_equal(colnames(assay_tab), c("#OTU ID", colnames(tse)))
    
    row_data <- read.table(
        file.path(dpath, "taxonomy.tsv"),
        sep = "\t",
        header = TRUE,
        comment.char = "",
        check.names = FALSE
    )
    
    expect_named(row_data, c("Feature ID", "Taxon", "Confidence"))
    expect_match(row_data$Taxon, "[a-z];_[a-z];_[a-z]")
    
    col_data <- read.table(
        file.path(dpath, "metadata.tsv"),
        sep = "\t",
        header = TRUE,
        comment.char = "",
        check.names = FALSE
    )
    
    expect_named(col_data, c("sample-id", colnames(colData(tse))))
    expect_identical(col_data[["sample-id"]], c("#q2:types", colnames(tse)))
    expect_in(col_data[1, -1], c("categorical", "numeric"))
    
    # Test exportMothur
    dpath <- tempfile()
    
    expect_warning(
        exportMothur(tse, dpath),
        "'group' is a reserved name in Mothur. colData variables with that name were made unique.",
        fixed = TRUE
    )
    
    expect_true(dir.exists(dpath))
    
    assay_tab <- read.table(
        file.path(dpath, "counts.tsv"),
        sep = "\t",
        header = TRUE,
        comment.char = "",
        check.names = FALSE
    )
    
    expect_equal(
        colnames(assay_tab),
        c("Representative_Sequence", "total", colnames(tse))
    )
    
    row_data <- read.table(
        file.path(dpath, "taxonomy.tsv"),
        sep = "\t",
        header = TRUE,
        comment.char = "",
        check.names = FALSE
    )
    
    expect_named(row_data, c("OTU", "Size", "Taxonomy"))
    expect_match(row_data$Taxonomy, "[a-z];[a-z];[a-z]")
    
    expect_identical(row_data$OTU, rownames(tse))
    expect_type(row_data$Size, "integer")
    
    col_data <- read.table(
        file.path(dpath, "metadata.tsv"),
        sep = "\t",
        header = TRUE,
        comment.char = "",
        check.names = FALSE
    )
    
    expect_named(col_data, c("group", "ID", "group.1"))
    expect_equal(colData(tse)$ID, col_data$ID)
})
