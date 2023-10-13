test_that("Estimate Alpha Diversity Indices with Rarefaction", {
    data(GlobalPatterns, package="mia")
    tse <- GlobalPatterns
    # Calculate the default Shannon index with 1 rarefaction round
    tse <- estimateAlphaWithRarefaction(tse)
    expect_true(any(grepl("shannon", colnames(colData(tse)))))

    # Calculate the default observed richness with 10 rarefaction rounds
    tse <- estimateAlphaWithRarefaction(tse, nrounds=10,
     FUN=mia::estimateRichness, args.fun=list(index="observed"))
    expect_true(any(grepl("observed", colnames(colData(tse)))))
})