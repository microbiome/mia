test_that("getDysbiosis and addDysbiosis", {
    data(Tengeler2020, package = "mia")
    tse <- Tengeler2020
    res <- getDysbiosis(tse, col.var = "patient_status", reference = "Control")
    expect_s4_class(res, "DataFrame")
    expect_identical(rownames(res), colnames(tse))
    expect_identical(colnames(res), "median_distance")
    # Reference: median Bray-Curtis to controls, excluding the sample itself
    diss <- as.matrix(vegan::vegdist(t(assay(tse, "counts")), "bray"))
    controls <- colnames(tse)[tse[["patient_status"]] == "Control"]
    ref <- vapply(colnames(tse), function(s){
        median(diss[s, setdiff(controls, s)])
    }, numeric(1))
    expect_equal(res[["median_distance"]], unname(ref))
    # Samples of the same group are excluded
    res <- getDysbiosis(
        tse, col.var = "patient_status", reference = "Control",
        group = "cohort")
    ref <- vapply(colnames(tse), function(s){
        same <- colnames(tse)[tse[["cohort"]] == tse[, s][["cohort"]]]
        median(diss[s, setdiff(controls, same)])
    }, numeric(1))
    expect_equal(res[["median_distance"]], unname(ref))
    # By default, the first level is used as controls
    tse[["patient_status"]] <- relevel(
        factor(tse[["patient_status"]]), ref = "Control")
    expect_equal(
        getDysbiosis(tse, col.var = "patient_status")[["median_distance"]],
        getDysbiosis(
            tse, col.var = "patient_status",
            reference = "Control")[["median_distance"]])
    # Scores are added to colData
    tse <- addDysbiosis(tse, col.var = "patient_status", name = "dysbiosis")
    expect_true("dysbiosis" %in% colnames(colData(tse)))
    # Scores are calculated and added for altExp
    altExp(tse, "genus") <- agglomerateByRank(Tengeler2020, rank = "Genus")
    tse <- addDysbiosis(
        tse, col.var = "patient_status", reference = "Control",
        altexp = "genus")
    expect_equal(
        altExp(tse, "genus")[["median_distance"]],
        getDysbiosis(
            altExp(tse, "genus"), col.var = "patient_status",
            reference = "Control")[["median_distance"]])
    # Input check
    expect_error(getDysbiosis(tse, col.var = "test"))
    expect_error(getDysbiosis(tse, col.var = "patient_status", reference = "x"))
    expect_error(getDysbiosis(tse[, tse[["patient_status"]] == "ADHD"],
        col.var = "patient_status"))
    expect_error(getDysbiosis(tse, col.var = "patient_status", index = "test"))
    expect_error(getDysbiosis(tse, col.var = "patient_status", group = "test"))
})
