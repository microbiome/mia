context("PERMANOVA")
# Sample data setup
data(GlobalPatterns, package = "mia")
tse <- GlobalPatterns
tse <- transformAssay(tse, method = "relabundance")

test_that("getPERMANOVA works on SummarizedExperiment", {
    # Test basic PERMANOVA without homogeneity check
    res <- getPERMANOVA(
        tse, assay.type = "relabundance",
        formula = x ~ SampleType,
        test.homogeneity = FALSE,
        permutations = 99
    )
    expect_s3_class(res, "anova.cca")
    
    # Test PERMANOVA with homogeneity check enabled
    res <- getPERMANOVA(
        tse, assay.type = "relabundance",
        formula = x ~ SampleType,
        test.homogeneity = TRUE,
        permutations = 99
    )
    expect_type(res, "list")
    expect_true("permanova" %in% names(res))
    expect_true("homogeneity" %in% names(res))
    expect_s3_class(res$permanova, "anova.cca")
    expect_s3_class(res$homogeneity, "data.frame")
    
    # Test full results with nested structure validation for detailed outputs
    res <- getPERMANOVA(
        tse, assay.type = "relabundance",
        formula = x ~ SampleType,
        test.homogeneity = TRUE,
        full = TRUE,
        permutations = 99
    )
    expect_type(res, "list")
    expect_s3_class(res[[2]][[2]][[1]][[1]], "betadisper")
    expect_s3_class(res[[2]][[2]][[1]][[2]], "permutest.betadisper")
})

test_that("addPERMANOVA stores results in metadata", {
    # Test that addPERMANOVA saves results in metadata as expected
    tse_with_meta <- addPERMANOVA(
        tse, assay.type = "relabundance",
        formula = x ~ SampleType,
        permutations = 99,
        name = "permanova_test"
    )
    metadata_res <- metadata(tse_with_meta)[["permanova_test"]]
    expect_type(metadata_res, "list")
    expect_true("permanova" %in% names(metadata_res))
    expect_true("homogeneity" %in% names(metadata_res))
    expect_s3_class(metadata_res$permanova, "anova.cca")
})

test_that("getPERMANOVA input validations", {
    # Check error handling for a nonexistent assay type
    expect_error(
        getPERMANOVA(tse, assay.type = "nonexistent", formula = x ~ SampleType)
    )
    # Check that an incorrect formula type raises an error
    expect_error(
        getPERMANOVA(tse, assay.type = "relabundance", formula = "SampleType"),
        "'formula' must be formula or NULL."
    )
    # Check that dimension mismatches between matrix and covariate data raise an error
    expect_error(
        getPERMANOVA(
            assay(tse[, 1:10]), formula = x ~ SampleType, data = colData(tse)
        ),
        "Number of columns in 'x' should match with number of rows in 'data'"
    )
    # Check that invalid homogeneity test settings raise an error
    expect_error(
        getPERMANOVA(
            tse, assay.type = "relabundance",
            formula = x ~ SampleType,
            test.homogeneity = "invalid"
        ),
        "'test.homogeneity' must be TRUE or FALSE"
    )
    
    # Invalid assay type (non-existent type)
    expect_error(
        getPERMANOVA(tse, assay.type = "nonexistent", formula = x ~ SampleType)
    )
    
    # Invalid variable
    expect_error(
        getPERMANOVA(tse, assay.type = "relabundance", col.var = "invalid")
    )
    
    # Incorrect permutations
    expect_error(
        getPERMANOVA(
            tse, assay.type = "relabundance", formula = x ~ SampleType, 
            permutations = -10
        )
    )
    
    # Incorrect 'by' parameter
    expect_error(
        getPERMANOVA(
            tse, assay.type = "relabundance", formula = x ~ SampleType, 
            by = "invalid_option"
        )
    )
    
    # Missing or incorrect 'homogeneity.test' option
    expect_error(
        getPERMANOVA(
            tse, assay.type = "relabundance", formula = x ~ SampleType, 
            homogeneity.test = "unsupported_test"
        )
    )
})

test_that("getPERMANOVA works on SingleCellExperiment", {
    # Convert tse to SingleCellExperiment format and test getPERMANOVA functionality
    sce <- as(tse, "SingleCellExperiment")
    res <- getPERMANOVA(
        sce, assay.type = "relabundance",
        formula = x ~ SampleType,
        permutations = 99
    )
    expect_s3_class(res$permanova, "anova.cca")
    expect_s3_class(res$homogeneity, "data.frame")
})

test_that("getPERMANOVA 'by' and 'homogeneity.test' options", {
    # Test 'by' parameter with 'terms' for PERMANOVA
    res_by_terms <- getPERMANOVA(
        tse, assay.type = "relabundance",
        formula = x ~ SampleType, by = "terms"
    )
    expect_s3_class(res_by_terms$permanova, "anova.cca")
    
    # Test homogeneity test with ANOVA option
    res_anova <- getPERMANOVA(
        tse, assay.type = "relabundance",
        formula = x ~ SampleType,
        homogeneity.test = "anova"
    )
    expect_s3_class(res_anova$homogeneity, "data.frame")
    
    # Test homogeneity test with Tukey HSD and full results
    res_tukey <- getPERMANOVA(
        tse, assay.type = "relabundance",
        formula = x ~ SampleType,
        homogeneity.test = "tukeyhsd",
        full = TRUE
    )
    expect_s3_class(res_tukey[[2]][[2]][[1]][[1]], "betadisper")
    expect_s3_class(res_tukey[[2]][[2]][[1]][[2]], "TukeyHSD")
})

test_that("getPERMANOVA handles edge cases", {
    # Test handling of a missing formula (default behavior)
    expect_no_error(getPERMANOVA(tse, assay.type = "relabundance"))
    
    # Test for error when permutations count is zero
    expect_error(
        getPERMANOVA(
            tse, assay.type = "relabundance",
            formula = x ~ SampleType, permutations = 0
        )
    )
    
    # Test for error when input matrix is NULL
    expect_error(getPERMANOVA(NULL, formula = x ~ SampleType))
    
    # Test for warning when only one level exists in the factor variable
    tse_subset <- tse[, tse$SampleType == "Soil"]
    expect_warning(
        getPERMANOVA(
            tse_subset, assay.type = "relabundance",
            formula = x ~ SampleType
        )
    )
})

test_that("getPERMANOVA matches direct calculations", {
    # Perform direct calculations with vegan package for comparison
    permanova_direct <- vegan::adonis2(
        t(assay(tse, "relabundance")) ~ SampleType,
        data = colData(tse),
        permutations = 99
    )
    homogeneity_direct <- vegan::betadisper(
        vegan::vegdist(t(assay(tse, "relabundance"))),
        group = tse$SampleType
    )
    
    # Run the getPERMANOVA function and compare results
    res <- getPERMANOVA(
        tse, assay.type = "relabundance",
        formula = x ~ SampleType,
        test.homogeneity = TRUE,
        full = TRUE,
        permutations = 99
    )
    
    # Verify permanova results match
    expect_equal(res$permanova$aov.tab, permanova_direct$aov.tab)
    
    # Verify homogeneity results match
    expect_equal(
        res[[2]][[2]][[1]][[1]]$distances, homogeneity_direct$distances)
})

test_that("getPERMANOVA works with UniFrac and dist object directly", {
    # Test on-the-fly UniFrac calculation on TreeSummarizedExperiment
    set.seed(42)
    res_unifrac <- getPERMANOVA(
        tse, assay.type = "counts",
        method = "unifrac",
        formula = x ~ SampleType,
        test.homogeneity = TRUE,
        permutations = 99
    )
    expect_type(res_unifrac, "list")
    expect_s3_class(res_unifrac$permanova, "anova.cca")
    expect_s3_class(res_unifrac$homogeneity, "data.frame")

    # Test passing a dist object directly
    d <- getDissimilarity(tse, method = "unifrac")
    expect_s3_class(d, "dist")
    setseed_wrap <- set.seed(42)
    res_dist <- suppressWarnings(getPERMANOVA(
        d, formula = x ~ SampleType, data = colData(tse),
        test.homogeneity = TRUE, permutations = 99
    ))
    expect_type(res_dist, "list")
    expect_s3_class(res_dist$permanova, "anova.cca")
    expect_equal(res_unifrac$permanova$SumOfSqs, res_dist$permanova$SumOfSqs)
    expect_equal(res_unifrac$permanova$R2, res_dist$permanova$R2)
    expect_equal(res_unifrac$permanova$F, res_dist$permanova$F)
})

test_that("getPERMANOVA and addPERMANOVA work with pre-calculated dissimilarity (dis.name)", {
    # Pre-calculate dissimilarity and store in metadata
    tse_diss <- addDissimilarity(tse, method = "unifrac", name = "unifrac_dist")
    expect_true("unifrac_dist" %in% names(metadata(tse_diss)))

    # Run getPERMANOVA using dis.name
    set.seed(42)
    res_precalc <- getPERMANOVA(
        tse_diss, dis.name = "unifrac_dist",
        formula = x ~ SampleType, permutations = 99
    )
    expect_s3_class(res_precalc$permanova, "anova.cca")

    d <- getDissimilarity(tse, method = "unifrac")
    set.seed(42)
    res_dist <- suppressWarnings(getPERMANOVA(
        d, formula = x ~ SampleType, data = colData(tse), permutations = 99
    ))
    expect_equal(res_precalc$permanova$SumOfSqs, res_dist$permanova$SumOfSqs)
    expect_equal(res_precalc$permanova$F, res_dist$permanova$F)

    # Test addPERMANOVA with dis.name
    tse_meta <- addPERMANOVA(
        tse_diss, dis.name = "unifrac_dist",
        formula = x ~ SampleType, name = "perm_precalc", permutations = 99
    )
    expect_true("perm_precalc" %in% names(metadata(tse_meta)))
    expect_s3_class(metadata(tse_meta)[["perm_precalc"]]$permanova, "anova.cca")

    # Test addPERMANOVA with method = unifrac on-the-fly
    tse_meta2 <- addPERMANOVA(
        tse, method = "unifrac",
        formula = x ~ SampleType, name = "perm_unifrac", permutations = 99
    )
    expect_true("perm_unifrac" %in% names(metadata(tse_meta2)))
    expect_s3_class(metadata(tse_meta2)[["perm_unifrac"]]$permanova, "anova.cca")
})

test_that("getPERMANOVA input validations for dis.name and dist objects", {
    # Nonexistent dis.name in metadata
    expect_error(
        getPERMANOVA(tse, dis.name = "nonexistent", formula = x ~ SampleType),
        "'dis.name' must be a valid name of metadata"
    )

    # Invalid dis.name type
    expect_error(
        getPERMANOVA(tse, dis.name = 123, formula = x ~ SampleType),
        "'dis.name' must be a single character value"
    )

    # Dist object with sample size mismatch
    d <- getDissimilarity(tse, method = "bray")
    expect_error(
        getPERMANOVA(d, formula = x ~ SampleType, data = colData(tse)[1:10, ]),
        "Number of samples in 'x' should match"
    )

    # Invalid x (neither matrix nor dist)
    expect_error(
        getPERMANOVA(1:10, formula = x ~ SampleType, data = colData(tse)),
        "'x' must be a matrix or a 'dist' object"
    )
})

