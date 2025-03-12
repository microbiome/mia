context("decontam")
test_that("decontam", {
    data(esophagus, package="mia")
    # setup of some mock data
    colData(esophagus)$concentration <- c(1,2,3)
    colData(esophagus)$control <- c(FALSE,FALSE,TRUE)
    #
    tse <- addContaminantQC(
        esophagus, method = "frequency", concentration = "concentration") |>
        expect_warning()
    expect_is(rowData(tse)[,"p"], "numeric")
    #
    tse <- addNotContaminantQC(
        esophagus, control = "control", detailed = FALSE) |>
        expect_warning()
    expect_is(rowData(tse)[, "not_contaminant"], "logical")
    #
    tse <- addNotContaminantQC(
        esophagus, control = "control", detailed = TRUE) |>
        expect_warning()
    expect_is(rowData(tse)[, "p"], "numeric")
    #
    res <- isContaminant(esophagus, control = "control") |> expect_warning()
    expect_s4_class(res, "DataFrame")
    #
    res <- isNotContaminant(esophagus, control = "control", detailed = FALSE) |>
        expect_warning()
    expect_is(res, "logical")
})
