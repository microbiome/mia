context("deprecate")

# Mutes the deprecation warning, so that other warnings can still be tested.
.quiet_deprecated <- function(expr){
    withCallingHandlers(expr, warning = function(w){
        if( grepl("is deprecated", conditionMessage(w)) ){
            invokeRestart("muffleWarning")
        }
    })
}

test_that("Deprecated prevalence functions give a deprecation warning", {
    data(GlobalPatterns, package = "mia")
    tse <- GlobalPatterns
    expect_warning(
        subsetByPrevalent(tse), "'subsetByPrevalent' is deprecated")
    expect_warning(subsetByRare(tse), "'subsetByRare' is deprecated")
    expect_warning(
        agglomerateByPrevalence(tse, rank = "Phylum"),
        "'agglomerateByPrevalence' is deprecated")
})

test_that("subsetByPrevalent and subsetByRare match subsetBy", {
    data(GlobalPatterns, package = "mia")
    tse <- GlobalPatterns
    tse <- transformAssay(tse, method = "relabundance")
    # Default thresholds: detection 0 and prevalence 50%
    expect_identical(
        .quiet_deprecated(subsetByPrevalent(tse)),
        subsetBy(tse, rows = rowMeans(counts > 0) > 0.5))
    expect_identical(
        .quiet_deprecated(subsetByRare(tse)),
        subsetBy(tse, rows = rowMeans(counts > 0) <= 0.5))
    # Detection, prevalence and assay
    res <- .quiet_deprecated(subsetByPrevalent(
        tse, assay.type = "relabundance", detection = 0.001,
        prevalence = 0.2))
    ref <- subsetBy(tse, rows = rowMeans(relabundance > 0.001) > 0.2)
    expect_identical(rownames(res), rownames(ref))
    # include.lowest includes the prevalence threshold
    res <- .quiet_deprecated(subsetByPrevalent(
        tse, prevalence = 0.5, include.lowest = TRUE))
    ref <- subsetBy(tse, rows = rowMeans(counts > 0) >= 0.5)
    expect_identical(rownames(res), rownames(ref))
    # Agglomeration by rank
    res <- .quiet_deprecated(subsetByPrevalent(
        tse, rank = "Family", detection = 1, prevalence = 0.2))
    ref <- agglomerateByRank(tse, rank = "Family")
    ref <- subsetBy(ref, rows = rowMeans(counts > 1) > 0.2)
    expect_identical(rownames(res), rownames(ref))
    # Prevalent and rare features together include all features
    prevalent <- .quiet_deprecated(subsetByPrevalent(tse, prevalence = 0.3))
    rare <- .quiet_deprecated(subsetByRare(tse, prevalence = 0.3))
    expect_setequal(c(rownames(prevalent), rownames(rare)), rownames(tse))
    expect_length(intersect(rownames(prevalent), rownames(rare)), 0L)
    # rowData variables do not mask the condition
    tse_keep <- tse
    rowData(tse_keep)[["keep"]] <- TRUE
    res <- .quiet_deprecated(subsetByPrevalent(tse_keep))
    ref <- subsetBy(tse, rows = rowMeans(counts > 0) > 0.5)
    expect_identical(rownames(res), rownames(ref))
})

test_that("agglomerateByPrevalence matches subsetBy with collapse", {
    data(GlobalPatterns, package = "mia")
    tse <- GlobalPatterns
    res <- .quiet_deprecated(agglomerateByPrevalence(
        tse, rank = "Phylum", detection = 1, prevalence = 0.5))
    ref <- agglomerateByRank(tse, rank = "Phylum")
    ref <- subsetBy(
        ref, rows = rowMeans(counts > 1) > 0.5, collapse = TRUE)
    expect_identical(rownames(res), rownames(ref))
    expect_equal(assay(res, "counts"), assay(ref, "counts"))
    # Merging into "Other" keeps the sample totals of the agglomerated data
    tse_phylum <- agglomerateByRank(tse, rank = "Phylum")
    expect_equal(
        colSums(assay(res, "counts")), colSums(assay(tse_phylum, "counts")))
    # The "Other" row has taxonomy only for the agglomeration rank
    other <- rowData(res)["Other", ]
    expect_identical(other[["Phylum"]], "Other")
    expect_true(all(is.na(unlist(
        other[, setdiff(colnames(other), "Phylum")]))))
    # Custom name for the merged feature
    res <- .quiet_deprecated(agglomerateByPrevalence(
        tse, rank = "Phylum", other.name = "Rest"))
    expect_identical(tail(rownames(res), 1L), "Rest")
    expect_identical(rowData(res)["Rest", "Phylum"], "Rest")
    # No "Other" row if all features are prevalent
    res <- .quiet_deprecated(agglomerateByPrevalence(
        tse, rank = "Kingdom", prevalence = 0))
    expect_false("Other" %in% rownames(res))
    # Works with SummarizedExperiment without rowData
    se <- SummarizedExperiment(
        assays = list(counts = assay(tse, "counts")[1:200, ]))
    res <- .quiet_deprecated(agglomerateByPrevalence(se, prevalence = 0.3))
    expect_identical(tail(rownames(res), 1L), "Other")
    expect_equal(
        colSums(assay(res, "counts")), colSums(assay(se, "counts")))
})

test_that("Deprecated prevalence functions check input", {
    data(GlobalPatterns, package = "mia")
    tse <- GlobalPatterns
    expect_error(
        .quiet_deprecated(subsetByPrevalent(tse, prevalence = "test")),
        "'prevalence' must be a single numeric value")
    expect_error(
        .quiet_deprecated(subsetByRare(tse, include.lowest = "test")),
        "'include.lowest' must be TRUE or FALSE")
    expect_error(
        .quiet_deprecated(agglomerateByPrevalence(tse, other.name = TRUE)),
        "'other.name' must be a single character value")
})

test_that("Defunct functions give an error naming the replacement", {
    data(GlobalPatterns, package = "mia")
    tse <- GlobalPatterns
    defunct <- c(
        cluster = "addCluster", getTopTaxa = "getTop",
        loadFromBiom = "importBIOM", mergeRows = "agglomerateByVariable",
        subsetByRareTaxa = "subsetBy", transformCounts = "transformAssay")
    for( fun in names(defunct) ){
        expect_error(
            get(fun)(tse), paste0("'", fun, "' is defunct. Use '",
                defunct[[fun]], "'"),
            class = "defunctError")
    }
    expect_error(relabundance(tse) <- 1, class = "defunctError")
})

test_that("Defunct argument names give an error naming the new argument", {
    data(GlobalPatterns, package = "mia")
    tse <- GlobalPatterns
    expect_error(
        getPrevalence(tse, assay_name = "counts"),
        "'assay_name' is defunct. Use 'assay.type' instead.")
    expect_error(
        getPrevalent(tse, include_lowest = TRUE),
        "'include_lowest' is defunct. Use 'include.lowest' instead.")
    expect_error(
        getDominant(tse, rank = "Genus"),
        "'rank' is defunct. Use 'group' instead.")
    expect_error(
        agglomerateByVariable(tse, by = "rows", f = "Phylum"),
        "'f' is defunct. Use 'group' instead.")
    # The new names give the same result as before
    expect_equal(
        getPrevalence(tse, assay.type = "counts"),
        getPrevalence(assay(tse, "counts")))
    expect_equal(
        getDominant(tse, group = "Genus"),
        getDominant(agglomerateByRank(tse, rank = "Genus")))
    expect_equal(
        agglomerateByVariable(tse, by = "rows", group = "Phylum"),
        agglomerateByVariable(tse, by = "rows", rowData(tse)$Phylum))
})

test_that("calculateDMN is deprecated for matrices and defunct for SE", {
    data(GlobalPatterns, package = "mia")
    tse <- agglomerateByRank(GlobalPatterns, rank = "Phylum")
    expect_error(calculateDMN(tse), "defunct")
    mat <- t(assay(tse, "counts"))
    expect_warning(res <- calculateDMN(mat, k = 1:2, seed = 1), "deprecated")
    expect_length(res, 2L)
    expect_s4_class(res[[1L]], "DMN")
})
