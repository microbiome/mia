context("decontam")
test_that("decontam", {
    data(esophagus)
    # setup of some mock data
    colData(esophagus)$concentration <- c(1,2,3)
    colData(esophagus)$control <- c(FALSE,FALSE,TRUE)
    expect_warning(esophagus <- isContaminant(esophagus,
                                              method = "frequency",
                                              concentration = "concentration"))
    expect_s4_class(rowData(esophagus)[,"isContaminant"],"DataFrame")
    expect_warning(esophagus <- isNotContaminant(esophagus, control = "control"))
    expect_type(rowData(esophagus)[,"isNotContaminant"],"logical")
    expect_warning(esophagus <- isNotContaminant(esophagus, control = "control",
                                                 detailed = TRUE))
    expect_s4_class(rowData(esophagus)[,"isNotContaminant"],"DataFrame")
})
