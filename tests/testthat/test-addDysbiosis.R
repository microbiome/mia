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

test_that("getDysbiosis scores adapted from dysbiosisR", {
    data(Tengeler2020, package = "mia")
    tse <- Tengeler2020
    ctrl <- tse[["patient_status"]] == "Control"
    index <- c("centroid_distance", "shannon_jsd", "reference_plane", "cloud")
    res <- getDysbiosis(tse, col.var = "patient_status", reference = "Control",
        index = index, k = 3)
    expect_identical(colnames(res), index)
    # Centroid distance equals the distance to group means for Euclidean
    # distances
    clr <- transformAssay(tse, method = "clr", pseudocount = 1)
    X <- t(assay(clr, "clr"))
    to_mean <- function(rows) sqrt(rowSums(sweep(X, 2, colMeans(X[rows, ]))^2))
    expect_equal(
        getDysbiosis(clr, col.var = "patient_status", reference = "Control",
            assay.type = "clr", index = "centroid_distance",
            method = "euclidean")[["centroid_distance"]],
        unname(to_mean(ctrl) - to_mean(!ctrl)))
    # Shannon-JSD: mean Shannon difference times mean JSD to other controls
    rel <- t(t(assay(tse)) / colSums(assay(tse)))
    shannon <- apply(rel, 2, function(p) -sum(p[p > 0] * log(p[p > 0])))
    jsd <- as.matrix(getDissimilarity(tse, method = "jsd"))
    ref <- vapply(seq_len(ncol(tse)), function(i){
        r <- ctrl
        r[i] <- FALSE
        mean(shannon[i] - shannon[r]) * mean(jsd[i, r])
    }, numeric(1))
    expect_equal(res[["shannon_jsd"]], ref)
    # CLOUD: mean distance to the 3 nearest controls relative to controls
    d <- as.matrix(vegan::vegdist(t(assay(tse)), "bray"))
    diag(d) <- NA
    nn <- apply(d[, ctrl], 1, function(v) mean(sort(v)[1:3]))
    expect_equal(res[["cloud"]], unname(nn / mean(nn[ctrl])))
    expect_equal(mean(res[["cloud"]][ctrl]), 1)
    # The reference plane distance is never negative
    expect_true(all(res[["reference_plane"]] >= 0))
    # Input check
    expect_error(getDysbiosis(tse, col.var = "patient_status",
        index = "cloud", k = 0))
})
