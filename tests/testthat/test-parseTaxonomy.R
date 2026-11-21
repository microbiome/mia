test_that("parseTaxonomy splits taxonomy strings into ranks", {
    taxa <- c(
        ASV1 = "k__Bacteria; p__Firmicutes; c__Bacilli; g__Lactobacillus",
        ASV2 = "k__Bacteria; p__Bacteroidetes; c__Bacteroidia")
    res <- parseTaxonomy(taxa, remove.prefix = TRUE)
    expect_s4_class(res, "DataFrame")
    expect_equal(colnames(res), c("kingdom", "phylum", "class", "genus"))
    expect_equal(rownames(res), c("ASV1", "ASV2"))
    expect_equal(res$phylum, c("Firmicutes", "Bacteroidetes"))
    expect_equal(res$genus, c("Lactobacillus", NA))
    # Prefixes are kept by default
    expect_equal(parseTaxonomy(taxa)$kingdom, rep("k__Bacteria", 2L))
    # Unnamed input gives no row names
    expect_null(rownames(parseTaxonomy(unname(taxa))))
    # Other separators
    res <- parseTaxonomy("k__Bacteria|p__Firmicutes", sep = "\\|",
        remove.prefix = TRUE)
    expect_equal(unlist(as.list(res[1L, ])),
        c(kingdom = "Bacteria", phylum = "Firmicutes"))
    # Tables
    df <- data.frame(id = 1:2, tax = unname(taxa), row.names = names(taxa))
    expect_equal(parseTaxonomy(df, col.name = "tax"), parseTaxonomy(taxa))
    expect_equal(parseTaxonomy(S4Vectors::DataFrame(df), col.name = "tax"),
        parseTaxonomy(taxa))
    # Same result as the importers' internal parser
    expect_equal(as.data.frame(parseTaxonomy(unname(taxa))),
        as.data.frame(mia:::.parse_taxonomy(data.frame(Taxon = unname(taxa)))),
        ignore_attr = TRUE)
    # Errors
    expect_error(parseTaxonomy(1:2), "'x' must be a character vector")
    expect_error(parseTaxonomy(df, col.name = "wrong"),
        "'col.name' must be a column of 'x'.")
    expect_error(parseTaxonomy(c("k__Bacteria", NA)), "must not contain NA")
})
