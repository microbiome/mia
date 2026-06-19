################################################################################
# Helper: create small fixture files for the data-processing helpers
################################################################################

# Minimal long-format MetaPhlAn4 profile (3 cols: sample, clade, abundance)
.make_assay_fixture <- function(dir) {
    path <- file.path(dir, "assay.tsv")
    lines <- c(
        "sample_alias\tclade_name\trel_abund",
        "S1\tk__Bacteria\t100.0",
        "S1\tt__SGB1234\t60.5",
        "S1\tt__SGB5678\t39.5",
        "S2\tt__SGB1234\t80.0",
        "S2\tt__SGB5678\t20.0",
        "S3\tt__SGB1234\t50.0",
        "S3\tt__SGB9999\t50.0"
    )
    writeLines(lines, path)
    path
}

# Minimal long-format metadata
.make_metadata_fixture <- function(dir) {
    path <- file.path(dir, "metadata.tsv")
    lines <- c(
        "sample_alias\tmetadata_item\tvalue",
        "S1\tage\t30",
        "S1\tcountry\tFI",
        "S2\tage\t45",
        "S2\tcountry\tDE",
        "S3\tage\t25",
        "S3\tcountry\tSE"
    )
    writeLines(lines, path)
    path
}

# Minimal taxonomy mapping database
.make_taxdb_fixture <- function(dir) {
    path <- file.path(dir, "taxdb.tsv")
    lines <- c(
        paste("clade_name", "NCBI_taxids", "lineage", sep = "\t"),
        paste("t__SGB1234", "12345",
            "k__Bacteria|p__Firmicutes|c__Bacilli|o__Lactobacillales|f__Lactobacillaceae|g__Lactobacillus|s__L_acidophilus|t__SGB1234",
            sep = "\t"),
        paste("t__SGB5678", "56789",
            "k__Bacteria|p__Firmicutes|c__Bacilli|o__Lactobacillales|f__Streptococcaceae|g__Streptococcus|s__S_thermophilus|t__SGB5678",
            sep = "\t"),
        paste("t__SGB9999", "99999",
            "k__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|o__Enterobacterales|f__Enterobacteriaceae|g__Escherichia|s__E_coli|t__SGB9999",
            sep = "\t"),
        paste("k__Bacteria", "2", "k__Bacteria", sep = "\t")
    )
    writeLines(lines, path)
    path
}

# Minimal sample list (csv)
.make_samplelist_fixture <- function(dir, samples = c("S1", "S2")) {
    path <- file.path(dir, "samplelist.csv")
    df <- data.frame(sample_alias = samples)
    write.csv(df, path, row.names = FALSE)
    path
}

################################################################################
# Input validation tests
################################################################################

test_that("fetchMetalogTSE rejects invalid collection", {
    expect_error(fetchMetalogTSE("invalid_collection"),
        "'collection' must be one of")
    expect_error(fetchMetalogTSE(123),
        "'collection' must be one of")
    expect_error(fetchMetalogTSE(c("human", "animal")),
        "'collection' must be one of")
})

test_that("fetchMetalogTSE rejects invalid meta.type", {
    expect_error(fetchMetalogTSE("human", meta.type = "bad"),
        "'meta.type' must be one of")
    expect_error(fetchMetalogTSE("human", meta.type = 42),
        "'meta.type' must be one of")
})

test_that("fetchMetalogTSE rejects invalid use.cache", {
    expect_error(fetchMetalogTSE("human", use.cache = "yes"),
        "'use.cache' must be TRUE or FALSE")
    expect_error(fetchMetalogTSE("human", use.cache = NA),
        "'use.cache' must be TRUE or FALSE")
})

test_that("fetchMetalogTSE rejects non-string samplelist", {
    expect_error(fetchMetalogTSE("human", samplelist = 123),
        "'samplelist' must be a single character value or NULL")
})

test_that("fetchMetalogTSE rejects non-existent samplelist", {
    expect_error(
        fetchMetalogTSE("human", samplelist = "no_such_file.csv"),
        "'samplelist' file does not exist")
})

test_that("fetchMetalogTSE rejects unsupported samplelist extension", {
    tmp <- tempfile(fileext = ".xlsx")
    writeLines("placeholder", tmp)
    on.exit(unlink(tmp))
    expect_error(fetchMetalogTSE("human", samplelist = tmp),
        "'samplelist' file type must be one of")
})

################################################################################
# .load_metalog_assay
################################################################################

test_that(".load_metalog_assay returns correct structure", {
    dir <- tempdir()
    path <- .make_assay_fixture(dir)
    on.exit(unlink(path))
    result <- mia:::.load_metalog_assay(path)
    expect_type(result, "list")
    expect_named(result, c("assay", "taxa", "samples"))
    expect_s4_class(result$assay, "dgCMatrix")
})

test_that(".load_metalog_assay filters to SGB rows only", {
    dir <- tempdir()
    path <- .make_assay_fixture(dir)
    on.exit(unlink(path))
    result <- mia:::.load_metalog_assay(path)
    # k__Bacteria row should be excluded
    expect_true(all(startsWith(result$taxa, "t__SGB")))
    expect_equal(length(result$taxa), 3L)
})

test_that(".load_metalog_assay has correct dimensions", {
    dir <- tempdir()
    path <- .make_assay_fixture(dir)
    on.exit(unlink(path))
    result <- mia:::.load_metalog_assay(path)
    # 3 taxa (SGB1234, SGB5678, SGB9999), 3 samples (S1, S2, S3)
    expect_equal(nrow(result$assay), 3L)
    expect_equal(ncol(result$assay), 3L)
    expect_equal(length(result$samples), 3L)
})

test_that(".load_metalog_assay aggregates duplicate entries", {
    dir <- tempdir()
    path <- file.path(dir, "assay_dup.tsv")
    lines <- c(
        "sample_alias\tclade_name\trel_abund",
        "S1\tt__SGB1234\t30.0",
        "S1\tt__SGB1234\t20.0"
    )
    writeLines(lines, path)
    on.exit(unlink(path))
    result <- mia:::.load_metalog_assay(path)
    # Should sum to 50.0
    expect_equal(as.numeric(result$assay["t__SGB1234", "S1"]), 50.0)
})

################################################################################
# .load_metalog_metadata
################################################################################

test_that(".load_metalog_metadata returns wide data.frame", {
    dir <- tempdir()
    path <- .make_metadata_fixture(dir)
    on.exit(unlink(path))
    result <- mia:::.load_metalog_metadata(path, c("S1", "S2"))
    expect_s3_class(result, "data.frame")
    expect_true("age" %in% colnames(result))
    expect_true("country" %in% colnames(result))
    expect_equal(nrow(result), 2L)
})

test_that(".load_metalog_metadata subsets to requested samples", {
    dir <- tempdir()
    path <- .make_metadata_fixture(dir)
    on.exit(unlink(path))
    result <- mia:::.load_metalog_metadata(path, c("S2"))
    expect_equal(nrow(result), 1L)
    expect_equal(rownames(result), "S2")
})

test_that(".load_metalog_metadata warns on missing samples", {
    dir <- tempdir()
    path <- .make_metadata_fixture(dir)
    on.exit(unlink(path))
    expect_warning(
        mia:::.load_metalog_metadata(path, c("S1", "MISSING")),
        "sample\\(s\\) present in assay data but missing"
    )
})

################################################################################
# .construct_metalog_taxmap
################################################################################

test_that(".construct_metalog_taxmap returns taxonomy data.frame", {
    dir <- tempdir()
    db_path <- .make_taxdb_fixture(dir)
    on.exit(unlink(db_path))
    taxa <- c("t__SGB1234", "t__SGB5678")
    result <- mia:::.construct_metalog_taxmap(db_path, taxa)
    expect_s3_class(result, "data.frame")
    expect_equal(nrow(result), 2L)
    expect_equal(rownames(result), taxa)
    expect_equal(
        colnames(result),
        c("Kingdom", "Phylum", "Class", "Order",
            "Family", "Genus", "Species", "SGB")
    )
})

test_that(".construct_metalog_taxmap parses lineage correctly", {
    dir <- tempdir()
    db_path <- .make_taxdb_fixture(dir)
    on.exit(unlink(db_path))
    result <- mia:::.construct_metalog_taxmap(db_path, "t__SGB1234")
    expect_equal(result["t__SGB1234", "Kingdom"], "k__Bacteria")
    expect_equal(result["t__SGB1234", "Phylum"], "p__Firmicutes")
    expect_equal(result["t__SGB1234", "Genus"], "g__Lactobacillus")
})

test_that(".construct_metalog_taxmap preserves taxa order", {
    dir <- tempdir()
    db_path <- .make_taxdb_fixture(dir)
    on.exit(unlink(db_path))
    taxa <- c("t__SGB5678", "t__SGB1234")
    result <- mia:::.construct_metalog_taxmap(db_path, taxa)
    expect_equal(rownames(result), taxa)
})

test_that(".construct_metalog_taxmap warns on unmatched taxa", {
    dir <- tempdir()
    db_path <- .make_taxdb_fixture(dir)
    on.exit(unlink(db_path))
    taxa <- c("t__SGB1234", "t__SGB0000")
    expect_warning(
        mia:::.construct_metalog_taxmap(db_path, taxa),
        "1 of 2 taxa could not be mapped"
    )
})

################################################################################
# .parse_metalog_date
################################################################################

test_that(".parse_metalog_date extracts date from filename", {
    expect_equal(
        mia:::.parse_metalog_date(
            "/cache/human_metaphlan4_2025-03-15.tsv.gz"),
        "2025-03-15"
    )
    expect_equal(
        mia:::.parse_metalog_date(
            "/cache/human_core_long_2024-12-01.tsv.gz"),
        "2024-12-01"
    )
})

test_that(".parse_metalog_date returns NA when no date present", {
    expect_true(is.na(
        mia:::.parse_metalog_date("some_file_without_date.tsv.gz")
    ))
})

################################################################################
# .filter_metalog_samples
################################################################################

test_that(".filter_metalog_samples subsets to requested samples", {
    dir <- tempdir()
    assay_path <- .make_assay_fixture(dir)
    sl_path <- .make_samplelist_fixture(dir, samples = c("S1", "S2"))
    on.exit(unlink(c(assay_path, sl_path)))
    assay_list <- mia:::.load_metalog_assay(assay_path)
    result <- mia:::.filter_metalog_samples(assay_list, sl_path)
    expect_equal(sort(result$samples), c("S1", "S2"))
    expect_equal(ncol(result$assay), 2L)
})

test_that(".filter_metalog_samples drops zero-abundance taxa", {
    dir <- tempdir()
    assay_path <- .make_assay_fixture(dir)
    # S1 and S2 have SGB1234 and SGB5678 but NOT SGB9999
    sl_path <- .make_samplelist_fixture(dir, samples = c("S1", "S2"))
    on.exit(unlink(c(assay_path, sl_path)))
    assay_list <- mia:::.load_metalog_assay(assay_path)
    result <- mia:::.filter_metalog_samples(assay_list, sl_path)
    # SGB9999 only in S3, so it should be dropped
    expect_false("t__SGB9999" %in% result$taxa)
    expect_equal(length(result$taxa), 2L)
})

test_that(".filter_metalog_samples errors when no samples match", {
    dir <- tempdir()
    assay_path <- .make_assay_fixture(dir)
    sl_path <- .make_samplelist_fixture(dir, samples = c("NONEXISTENT"))
    on.exit(unlink(c(assay_path, sl_path)))
    assay_list <- mia:::.load_metalog_assay(assay_path)
    expect_error(
        mia:::.filter_metalog_samples(assay_list, sl_path),
        "None of the samples"
    )
})
