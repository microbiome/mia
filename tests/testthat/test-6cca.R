context("CCA")
test_that("CCA", {
    # .remove_special_functions_from_terms
    expect_error(mia:::.remove_special_functions_from_terms(),
                 'argument "terms" is missing, with no default')
    expect_equal(mia:::.remove_special_functions_from_terms("abc"),
                 c(abc = "abc"))
    expect_equal(mia:::.remove_special_functions_from_terms("Condition(abc)"),
                 c("Condition(abc)" = "abc"))
    expect_equal(mia:::.remove_special_functions_from_terms(c("abc","def")),
                 c(abc = "abc", def = "def"))
    expect_equal(mia:::.remove_special_functions_from_terms(c("Condition(abc)","def")),
                 c("Condition(abc)" = "abc", def = "def"))
    #
    skip_if_not(requireNamespace("vegan", quietly = TRUE))
    data(dune, dune.env, package = "vegan")
    sce <- SingleCellExperiment(assays = list(counts = t(dune)),
                               colData = DataFrame(dune.env))
    # .get_variables_from_data_and_formula
    expect_null(mia:::.get_variables_from_data_and_formula())
    form <- dune ~ Condition(Management) + Manure + A1
    expect_error(mia:::.get_variables_from_data_and_formula(formula = form),
                 'argument "x" is missing, with no default')
    actual <- mia:::.get_variables_from_data_and_formula(sce, form)
    expect_s4_class(actual, "DataFrame")
    expect_named(actual, c("Management", "Manure", "A1"))
    # .get_dependent_var_name
    expect_error(mia:::.get_dependent_var_name(),
                 'argument "formula" is missing, with no default')
    actual <- mia:::.get_dependent_var_name(form)
    expect_equal(actual, "dune")
    #
    mcca <- vegan::cca(form, dune.env, scale = TRUE)
    mrda <- vegan::rda(form, dune.env, scale = FALSE)
    sce <- runCCA(sce, form)
    actual <- reducedDim(sce,"CCA")
    expect_equal(as.vector(actual), as.vector(mcca$CCA$u))
    sce <- runRDA(sce, form)
    actual <- reducedDim(sce,"RDA")
    expect_equal(abs( as.vector(actual) ), abs( as.vector(mrda$CCA$u) ))
    sce <- runRDA(sce, form, distance = "bray", name = "rda_bray")
    actual <- reducedDim(sce,"rda_bray")
    rda_bray <- vegan::dbrda(form, dune.env, distance = "bray")
    expect_equal(abs( as.vector(actual) ), abs( as.vector(rda_bray$CCA$u) ))
    #
    sce <- runRDA(sce)
    test <- reducedDim(sce,"RDA")
    # Test that eigenvalues match
    test <- attr(test, "rda")$CA$eig
    res <- vegan::rda(t(assay(sce)))$CA$eig
    expect_equal(unname(test), unname(res))
    data(GlobalPatterns, package="mia")
    GlobalPatterns <- estimateDiversity(GlobalPatterns, index = "shannon")
    expect_error(calculateRDA(GlobalPatterns, variables = c("Primer", "test")))
    res1 <- calculateRDA(GlobalPatterns, variables = c("shannon", "SampleType"))
    res1 <- attr(res1, "rda")$CCA
    res2 <- calculateRDA(GlobalPatterns, formula = data ~ shannon + SampleType)
    res2 <- attr(res2, "rda")$CCA
    expect_equal(res1, res2)
})
