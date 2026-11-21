context("mergeSEs")

# Small objects with known content. Features a-e are tips of the same tree.
.tree <- ape::read.tree(
    text = "((a:1,b:1)n1:1,(c:1,(d:1,e:1)n3:1)n2:1)root;")
.seqs <- Biostrings::DNAStringSet(
    c(a = "AAAA", b = "CCCC", c = "GGGG", d = "TTTT", e = "ACGT"))
.genus <- c(a = "G1", b = "G1", c = "G2", d = "G3", e = "G3")

.make_se <- function(rows, cols, start = 1L, rowdata = TRUE){
    mat <- matrix(
        seq(start, length.out = length(rows) * length(cols)),
        nrow = length(rows), dimnames = list(rows, cols))
    rd <- S4Vectors::DataFrame(row.names = rows)
    if( rowdata ){
        rd$Genus <- unname(.genus[rows])
    }
    cd <- S4Vectors::DataFrame(group = cols, row.names = cols)
    SummarizedExperiment(
        assays = list(counts = mat), rowData = rd, colData = cd)
}

.make_tse <- function(rows, cols, start = 1L, tree = TRUE, seqs = TRUE){
    tse <- as(.make_se(rows, cols, start), "TreeSummarizedExperiment")
    if( tree ){
        rowTree(tse) <- .tree
    }
    if( seqs ){
        referenceSeq(tse) <- .seqs[rows]
    }
    tse
}

test_that("mergeSEs checks input", {
    se1 <- .make_se(c("a", "b"), c("s1", "s2"))
    se2 <- .make_se(c("b", "c"), c("s3", "s4"))
    # Wrong input objects
    expect_error(mergeSEs(se1))
    expect_error(mergeSEs(se1, se2, se1))
    expect_error(mergeSEs(list(se1, "test")), "not 'SummarizedExperiment'")
    se_no_names <- se1
    rownames(se_no_names) <- NULL
    expect_error(mergeSEs(se1, se_no_names), "rownames and colnames")
    # Wrong arguments
    for( join in list(1, TRUE, NA, "test", c("full", "inner")) ){
        expect_error(mergeSEs(se1, se2, join = join), "'join'")
    }
    for( join in c("left", "right") ){
        expect_error(mergeSEs(list(se1, se2, se1), join = join), "not available")
    }
    for( missing.values in list(TRUE, 1, c(0, 0), NULL) ){
        expect_error(
            mergeSEs(se1, se2, missing.values = missing.values),
            "'missing.values'")
    }
    for( value in list(NA, 1, "test", NULL) ){
        expect_error(mergeSEs(se1, se2, collapse.cols = value), "collapse.cols")
        expect_error(mergeSEs(se1, se2, collapse.rows = value), "collapse.rows")
    }
    # Assays must be found from all the objects
    expect_error(mergeSEs(se1, se2, assay.type = "test"), "'assay.type'")
    assay(se1, "relabundance") <- assay(se1)
    expect_warning(
        tse <- mergeSEs(se1, se2, assay.type = c("counts", "relabundance")),
        "'relabundance'")
    expect_equal(assayNames(tse), "counts")
})

test_that("mergeSEs gives same result for all input formats", {
    se1 <- .make_se(c("a", "b"), c("s1", "s2"))
    se2 <- .make_se(c("b", "c"), c("s3", "s4"))
    ref <- mergeSEs(se1, se2)
    expect_equal(mergeSEs(list(se1, se2)), ref)
    expect_equal(mergeSEs(SimpleList(se1, se2)), ref)
    # Deprecated argument names
    expect_equal(mergeSEs(se1, se2, assay_name = "counts"), ref)
})

test_that("mergeSEs merges assays", {
    se1 <- .make_se(c("a", "b"), c("s1", "s2"), start = 1L)
    se2 <- .make_se(c("b", "c"), c("s3"), start = 11L)
    tse <- mergeSEs(se1, se2)
    expect_equal(dimnames(tse), list(c("a", "b", "c"), c("s1", "s2", "s3")))
    ref <- matrix(
        c(1L, 2L, NA, 3L, 4L, NA, NA, 11L, 12L), nrow = 3L,
        dimnames = dimnames(tse))
    expect_identical(assay(tse), ref)
    # Missing values
    ref[is.na(ref)] <- 0L
    expect_equal(assay(mergeSEs(se1, se2, missing.values = 0)), ref)
    ref <- assay(mergeSEs(se1, se2, missing.values = "MISSING"))
    expect_equal(ref[is.na(assay(tse))], rep("MISSING", 3L))
    # Multiple assays
    assay(se1, "relabundance") <- assay(se1) / 10
    assay(se2, "relabundance") <- assay(se2) / 10
    tse <- mergeSEs(se1, se2, assay.type = c("counts", "relabundance"))
    expect_equal(assayNames(tse), c("counts", "relabundance"))
    expect_equal(assay(tse, "relabundance"), assay(tse, "counts") / 10)
})

test_that("mergeSEs selects rows based on join", {
    se1 <- .make_se(c("a", "b", "c"), "s1")
    se2 <- .make_se(c("b", "c", "d"), "s2")
    se3 <- .make_se(c("c", "d", "e"), "s3")
    expect_equal(rownames(mergeSEs(se1, se2, join = "full")),
        c("a", "b", "c", "d"))
    expect_equal(rownames(mergeSEs(se1, se2, join = "inner")), c("b", "c"))
    expect_equal(rownames(mergeSEs(se1, se2, join = "left")), rownames(se1))
    expect_equal(rownames(mergeSEs(se1, se2, join = "right")), rownames(se2))
    tse <- mergeSEs(list(se1, se2, se3), join = "inner")
    expect_equal(dimnames(tse), list("c", c("s1", "s2", "s3")))
    # No shared rows
    tse <- mergeSEs(se1, .make_se("e", "s2"), join = "inner")
    expect_equal(dim(tse), c(0L, 2L))
})

test_that("mergeSEs collapses columns and rows", {
    se1 <- .make_se(c("a", "b"), c("s1", "s2"), start = 1L)
    se2 <- .make_se(c("b", "c"), c("s2", "s3"), start = 11L)
    # Equally named samples are unique by default
    tse <- mergeSEs(se1, se2)
    expect_equal(colnames(tse), c("s1", "s2", "s2.1", "s3"))
    # Collapsed samples. The value of the first object is kept.
    tse <- mergeSEs(se1, se2, collapse.cols = TRUE)
    expect_equal(colnames(tse), c("s1", "s2", "s3"))
    expect_equal(assay(tse)["b", "s2"], assay(se1)["b", "s2"])
    expect_equal(assay(tse)["c", "s2"], assay(se2)["c", "s2"])
    # Equally named features are collapsed by default. If not, they are
    # made unique.
    tse <- mergeSEs(se1, se2, collapse.rows = FALSE)
    expect_equal(rownames(tse), c("a", "b", "b.1", "c"))
    expect_equal(rowData(tse)$Genus, unname(.genus[c("a", "b", "b", "c")]))
    # Duplicated names within an object are made unique
    se_dupl <- se1
    dimnames(se_dupl) <- list(c("a", "a"), c("s", "s"))
    tse <- mergeSEs(se_dupl, se2)
    expect_equal(dimnames(tse), list(
        c("a", "a.1", "b", "c"), c("s", "s.1", "s2", "s3")))
})

test_that("mergeSEs matches rows based on rownames and taxonomy", {
    se1 <- .make_se(c("a", "b"), "s1")
    se2 <- .make_se(c("a", "b"), "s2")
    rowData(se2)["b", "Genus"] <- "other"
    tse <- mergeSEs(se1, se2)
    expect_equal(rownames(tse), c("a", "b", "b.1"))
    expect_equal(rowData(tse)$Genus, c("G1", "G1", "other"))
    expect_equal(unname(is.na(assay(tse))), cbind(
        c(FALSE, FALSE, TRUE), c(FALSE, TRUE, FALSE)))
    # Missing rank is equal to missing value. Ranks are matched
    # case-insensitively and named by taxonomy ranks.
    se3 <- .make_se(c("a", "b"), "s3", rowdata = FALSE)
    rowData(se3)$genus <- c("G1", NA)
    rowData(se1)["b", "Genus"] <- NA
    se4 <- .make_se(c("a", "b"), "s4", rowdata = FALSE)
    tse <- mergeSEs(list(se1, se3, se4))
    expect_equal(rownames(tse), c("a", "b", "a.1"))
    expect_equal(colnames(rowData(tse)), "Genus")
    expect_equal(rowData(tse)$Genus, c("G1", NA, NA))
    expect_false(anyNA(assay(tse)[c("a", "b"), c("s1", "s3")]))
})

test_that("mergeSEs merges rowData and colData", {
    se1 <- .make_se(c("a", "b"), c("s1", "s2"))
    se2 <- .make_se(c("b", "c"), c("s2", "s3"))
    colData(se1)$only1 <- c(1, NA)
    colData(se2)$only1 <- c(10, 20)
    colData(se2)$only2 <- c("x", "y")
    tse <- mergeSEs(se1, se2, collapse.cols = TRUE)
    # Columns of all the objects are included. The first non-missing value is
    # used.
    expect_equal(colnames(colData(tse)), c("group", "only1", "only2"))
    expect_equal(colData(tse)$only1, c(1, 10, 20))
    expect_equal(colData(tse)$only2, c(NA, "x", "y"))
    # Objects without rowData
    se3 <- .make_se("a", "s1", rowdata = FALSE)
    tse <- mergeSEs(se3, se3)
    expect_equal(dim(rowData(tse)), c(1L, 0L))
    expect_equal(rownames(rowData(tse)), "a")
})

test_that("mergeSEs combines variables with different classes", {
    se1 <- .make_se("a", "s1")
    se2 <- .make_se("a", "s2")
    cd1 <- DataFrame(
        fct = factor("x", levels = c("x", "z")), num = 1L, chr = 1,
        row.names = "s1")
    cd2 <- DataFrame(
        fct = factor("y"), num = 2.5, chr = "b", row.names = "s2")
    colData(se1) <- cd1
    colData(se2) <- cd2
    cd <- colData(mergeSEs(se1, se2))
    # Factors get all the levels
    expect_equal(cd$fct, factor(c("x", "y"), levels = c("x", "z", "y")))
    # Numeric values stay numeric
    expect_equal(cd$num, c(1, 2.5))
    # Other combinations are converted to character
    expect_equal(cd$chr, c("1", "b"))
})

test_that("mergeSEs combines metadata", {
    se1 <- .make_se("a", "s1")
    se2 <- .make_se("a", "s2")
    metadata(se1) <- list(m1 = 1)
    metadata(se2) <- list(m2 = "x")
    expect_equal(metadata(mergeSEs(se1, se2)), list(m1 = 1, m2 = "x"))
})

test_that("mergeSEs returns the shared class of input", {
    se <- .make_se("a", "s1")
    sce <- as(se, "SingleCellExperiment")
    tse <- .make_tse("a", "s1", tree = FALSE, seqs = FALSE)
    expect_s4_class(mergeSEs(tse, tse), "TreeSummarizedExperiment")
    expect_identical(class(mergeSEs(sce, sce))[[1L]], "SingleCellExperiment")
    expect_identical(class(mergeSEs(se, se))[[1L]], "SummarizedExperiment")
    expect_warning(res <- mergeSEs(tse, sce), "multiple classes")
    expect_identical(class(res)[[1L]], "SingleCellExperiment")
    expect_warning(res <- mergeSEs(list(tse, sce, se)), "multiple classes")
    expect_identical(class(res)[[1L]], "SummarizedExperiment")
})

test_that("mergeSEs uses sparse matrix only if output is large", {
    se1 <- .make_se(c("a", "b"), c("s1", "s2"), start = 1L)
    se2 <- .make_se(c("b", "c"), c("s2", "s3"), start = 11L)
    dense <- mergeSEs(se1, se2, missing.values = 0, collapse.cols = TRUE)
    expect_true(is.matrix(assay(dense)))
    # Large output, i.e., output that has more cells than 'dense'
    sparse <- mergeSEs(
        se1, se2, missing.values = 0, collapse.cols = TRUE, dense = 0)
    expect_s4_class(assay(sparse), "dgCMatrix")
    expect_equal(as.matrix(assay(sparse)), assay(dense) * 1)
    # Sparse matrix cannot hold NAs
    expect_true(is.matrix(assay(mergeSEs(se1, se2, dense = 0))))
})

test_that("mergeSEs adds tree only if single tree includes all the rows", {
    tse1 <- .make_tse(c("a", "b"), "s1")
    tse2 <- .make_tse(c("c", "e"), "s2")
    # The tree of both objects includes all the rows
    tse <- mergeSEs(tse1, tse2)
    expect_equal(rowTree(tse), .tree)
    expect_equal(rowLinks(tse)$nodeLab, c("a", "b", "c", "e"))
    # Only the second tree includes all the rows
    tse1 <- .make_tse(c("a", "b"), "s1")
    rowTree(tse1) <- ape::keep.tip(.tree, c("a", "b"))
    tse <- mergeSEs(tse1, tse2)
    expect_equal(rowTree(tse), .tree)
    # Rows of an object without tree are not linked to the tree
    expect_warning(
        tse <- mergeSEs(tse2, .make_tse("d", "s3", tree = FALSE)),
        "rowTree is dropped")
    expect_null(rowTree(tse))
    # None of the trees includes all the rows
    tse2 <- .make_tse(c("c", "e"), "s2", tree = FALSE)
    rowTree(tse2) <- ape::keep.tip(.tree, c("c", "e"))
    expect_warning(tse <- mergeSEs(tse1, tse2), "rowTree is dropped")
    expect_null(rowTree(tse))
    # No rows
    tse <- mergeSEs(tse1, tse2, join = "inner")
    expect_null(rowTree(tse))
    # colTree
    tse1 <- .make_tse("a", c("a", "b"), tree = FALSE, seqs = FALSE)
    tse2 <- .make_tse("a", c("c", "d"), tree = FALSE, seqs = FALSE)
    colTree(tse1) <- .tree
    colTree(tse2) <- ape::keep.tip(.tree, c("c", "d"))
    tse <- mergeSEs(tse1, tse2)
    expect_equal(colTree(tse), .tree)
    expect_null(rowTree(tse))
})

test_that("mergeSEs merges reference sequences", {
    tse1 <- .make_tse(c("a", "b"), "s1", tree = FALSE)
    tse2 <- .make_tse(c("b", "e"), "s2", tree = FALSE)
    tse <- mergeSEs(tse1, tse2)
    expect_equal(referenceSeq(tse), .seqs[c("a", "b", "e")])
    # Sequences that are not available for all the rows are dropped
    tse3 <- .make_tse("d", "s3", tree = FALSE, seqs = FALSE)
    expect_warning(tse <- mergeSEs(tse1, tse3), "referenceSeq is dropped")
    expect_null(referenceSeq(tse))
    # No sequences
    expect_null(referenceSeq(mergeSEs(tse3, tse3)))
})
