
context("estimates")
test_that("estimates", {
    data("esophagus", package = "MicrobiomeExperiment")
    FUN_list <- list(breakaway = estimateBreakaway,
                     chao1 = estimateChao1,
                     chao_bunge = estimateChaoBunge,
                     kemp = estimateKemp,
                     poisson = estimatePoissonModel,
                     WLRM = estimateWLRMtransformed,
                     WLRMun = estimateWLRMuntransformed,
                     InvSimpson = estimateInvSimpson,
                     Simpson = estimateSimpson,
                     Shannon = estimateShannon,
                     ShannonE = estimateShannonE,
                     Richness = estimateRichness)
    for(FUN in FUN_list){
        suppressWarnings(esophagus <- do.call(FUN,list(esophagus)))
    }
    cd <- colData(esophagus)
    expect_true(all(names(FUN_list) %in% colnames(cd)))
    expect_equal(unname(round(cd$breakaway, 4)), c(38.2520, 46.5959, 184.5630))
    expect_equal(unname(round(cd$chao1, 4)), c(42.0833, 38.7857, 78.3333))
    expect_equal(unname(round(cd$chao_bunge, 3)), c(38.252, 43.387, 634.750))
    expect_equal(unname(round(cd$kemp, 3)), c(NA, NA, 184.563))
    expect_equal(unname(round(cd$poisson, 4)), c(33.3355, 34.4299, 45.3274))
    expect_equal(unname(round(cd$WLRM, 4)), c(NA, 41.4148, 108.2979))
    expect_equal(unname(round(cd$WLRMun, 4)), c(NA, 46.5959, NA))
    expect_equal(unname(round(cd$InvSimpson, 5)), c(5.93021, 10.34606, 2.99177))
    expect_equal(unname(round(cd$Simpson, 7)), c(0.1686282, 0.0966551, 0.3342507))
    expect_equal(unname(round(cd$Shannon, 5)), c(2.24937, 2.76239, 2.03249))
    expect_equal(unname(round(cd$ShannonE, 6)), c(0.675039, 0.790042, 0.558748))
    expect_equal(unname(cd$Richness), c(28, 33, 38))
})
