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
    # Check that input check of scores is working
    expect_error(runCCA(sce, form, scores = 1))
    expect_error(runCCA(sce, form, scores = TRUE))
    expect_error(runCCA(sce, form, scores = NULL))
    expect_error(runCCA(sce, form, scores = c("wa", "u")))
    # Remove "v" ption from parameters, it is already in rotation atrribute.
    # 
    mcca <- vegan::cca(form, dune.env)
    sce <- runCCA(sce, form)
    actual <- reducedDim(sce,"CCA")
    expect_equal(as.vector(actual), as.vector(mcca$CCA$wa))
    sce <- runCCA(sce, form, scores = "v") # With tes.signif = F this works...
    actual <- reducedDim(sce,"CCA")
    expect_equal(as.vector(actual), as.vector(mcca$CCA$v))
    #
    mcca <- vegan::cca(form, dune.env, scale = TRUE)
    mrda <- vegan::rda(form, dune.env, scale = FALSE)
    
    sce <- runCCA(sce, form, scores = "u")
    actual <- reducedDim(sce,"CCA")
    expect_equal(as.vector(actual), as.vector(mcca$CCA$u))
    # Check that test.signif works
    expect_error( runCCA(sce, test.signif = 1) )
    expect_error( runCCA(sce, test.signif = "TRUE") )
    expect_error( runCCA(sce, test.signif = NULL) )
    expect_error( runCCA(sce, test.signif = c(TRUE, TRUE)) )
    mat <- calculateRDA(sce, scores = "u", test.signif = FALSE)
    expect_true(is.null(attributes(mat)$significance))
    # Check that significance calculations are correct
    set.seed(46)
    sce <- runCCA(sce, variables = "Manure", full = TRUE)
    actual <- reducedDim(sce,"CCA")
    res <- attributes(actual)$significance
    # Create a function that calculates significances
    calc_signif <- function(obj, assay, betadisp_group){
        permanova <- vegan::anova.cca(obj, permutations = 999)
        permanova2 <- vegan::anova.cca(obj, permutations = 999, by = "margin")
        dist <- vegan::vegdist(t(assay), method = "euclidean")
        betadisp <- vegan::betadisper(dist, group = betadisp_group)
        betadisp_permanova <- vegan::permutest(betadisp, permutations = 999)
        betadisp_anova <- anova(betadisp)
        betadisp_tukeyhsd <- TukeyHSD(betadisp)
        res <- list(
            permanova = permanova,
            permanova_variables = permanova2,
            betadisper = betadisp,
            betadisper_permanova = betadisp_permanova,
            betadisper_anova = betadisp_anova,
            betadisper_tukeyhsd = betadisp_tukeyhsd
        )
        return(res)
    }
    set.seed(46)
    test <- calc_signif(attributes(actual)$cca, assay(sce), colData(sce)[["Manure"]])
    # Permanova
    expect_equal(res$permanova$model, test$permanova)
    expect_equal(res$permanova$variables, test$permanova_variables)
    # Betadisper (homogeneity of groups)
    expect_equal(res$homogeneity$variables$Manure$betadisper$vectors, test$betadisper$vectors)
    # Significance of betadisper with different permanova, anova and tukeyhsd
    expect_equal(res$homogeneity$variables$Manure$permanova, test$betadisper_permanova)
    set.seed(46)
    sce <- runCCA(sce, form, full = TRUE, homogeneity.test = "anova")
    actual <- reducedDim(sce,"CCA")
    res <- attributes(actual)$significance
    expect_equal(res$homogeneity$variables$Manure$anova, test$betadisper_anova)
    set.seed(46)
    sce <- runCCA(sce, form, full = TRUE, homogeneity.test = "tukeyhsd")
    actual <- reducedDim(sce,"CCA")
    res <- attributes(actual)$significance
    expect_equal(res$homogeneity$variables$Manure$tukeyhsd, test$betadisper_tukeyhsd)
    # Test that data is subsetted correctly
    data("enterotype", package = "mia")
    variable_names <- c("ClinicalStatus", "Gender", "Age")
    res <- runRDA(enterotype, variables = variable_names, na.action = na.exclude)
    expect_equal(colnames(res), colnames(enterotype))
    res <- runRDA(enterotype, variables = variable_names, na.action = na.exclude, subset_result = TRUE)
    enterotype <- enterotype[, complete.cases(colData(enterotype)[, variable_names])]
    expect_equal(colnames(res), colnames(enterotype))
    #
    sce <- runRDA(sce, form, scores = "u")
    actual <- reducedDim(sce,"RDA")
    expect_equal(abs( as.vector(actual) ), abs( as.vector(mrda$CCA$u) ))
    sce <- runRDA(sce, form, distance = "bray", name = "rda_bray", scores = "u")
    actual <- reducedDim(sce,"rda_bray")
    rda_bray <- vegan::dbrda(form, dune.env, distance = "bray")
    expect_equal(abs( as.vector(actual) ), abs( as.vector(rda_bray$CCA$u) ))
    #
    sce <- runRDA(sce, scores = "u")
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
