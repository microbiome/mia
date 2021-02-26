context("estimateEvenness")

test_that("estimateEvenness", {

    test_internal_estimateEvenness <- function(mat){
        expect_equal(round(as.vector(mia:::.run_evenness(mat, index = "camargo")), 7),
                     round(c(0.6942294, 0.6230541, 0.8010094)), 7)

        expect_equal(round(as.vector(mia:::.run_evenness(mat, index="pielou")),7),
                     round(c(0.6750387, 0.7900423, 0.5587478),7))

        expect_equal(round(as.vector(mia:::.run_evenness(mat, index="simpson")), 7),
                     round(c(0.21179306, 0.31351703, 0.07873068), 7))

        expect_equal(round(as.vector(mia:::.run_evenness(mat, index="evar")), 7),
                     round(c(0.3723086, 0.4073989, 0.4820153), 7))

        expect_equal(round(as.vector(mia:::.run_evenness(mat, index="bulla")), 7),
                     round(c(0.3627075, 0.4897059, 0.3519684), 7))
    }

    # TSE object
    data("esophagus")
    tse <- esophagus
    test_internal_estimateEvenness(assay(esophagus,"counts"))
    assay(tse,"counts") <- DelayedArray(assay(tse,"counts"))
    test_internal_estimateEvenness(assay(tse,"counts"))
})
