context("calculateDistance")
test_that("calculateDistance", {
    library(scater)
    mat <- matrix(1:60, nrow = 6)
    df <- DataFrame(n = c(1:6))
    se <- SummarizedExperiment(assays = list(counts = mat),
                               rowData = df)
    # default
    actual <- calculateDistance(se)
    actual2 <- dist(t(assay(se,"counts")))
    expect_true(all(actual == actual2))
    #
    data(enterotype)
    data(esophagus)
    # default
    actual <- calculateDistance(enterotype)
    expect_equal(actual[1], 0.22382452, tolerance = .0000001)
    expect_equal(actual[length(actual)], 0.3316593, tolerance = .0000001)
    # JSD
    # actual <- calculateJSD(enterotype)
    # expect_equal(actual[1], 0.1282428, tolerance = .0000001)
    # expect_equal(actual[length(actual)], 0.06595946, tolerance = .0000001)
    # UniFrac
    actual <- calculateUniFrac(esophagus, weighted = FALSE, normalized = TRUE)
    expect_equal(actual[1], 0.5175550, tolerance = .0000001)
    expect_equal(actual[length(actual)], 0.5422394, tolerance = .0000001)
    actual <- calculateUniFrac(esophagus, weighted = TRUE, normalized = TRUE)
    expect_equal(actual[1], 0.2035424, tolerance = .0000001)
    expect_equal(actual[length(actual)], 0.2477016, tolerance = .0000001)
    actual <- calculateUniFrac(esophagus, weighted = TRUE, normalized = FALSE)
    expect_equal(actual[1], 0.1050480, tolerance = .0000001)
    expect_equal(actual[length(actual)], 0.1422409, tolerance = .0000001)
    #
    enterotype <- runMDS(enterotype, FUN = calculateJSD, name = "JSD",
                         exprs_values = "counts")
    expect_named(reducedDims(enterotype),"JSD")
    jsd <- matrix(c(-0.2993825,-0.1876141,-0.1761079,-0.1508764,-0.2822249,
                    -0.3040909,-0.14425004,-0.13490508,-0.02830419,-0.06346108,
                    -0.11871798,-0.11856633), ncol=2)
    rownames(jsd) <- colnames(enterotype)[seq_len(nrow(jsd))]
    expect_equal(round(head(reducedDim(enterotype,"JSD")),7),
                 round(jsd,7))
    esophagus <- runMDS(esophagus, FUN = calculateUniFrac, name = "UniFrac",
                        exprs_values = "counts",
                        tree = rowTree(esophagus))
    expect_named(reducedDims(esophagus),"UniFrac")
    # Test overlap
    tse <- enterotype
    overlap <- calculateOverlap(tse, detection = 3)
    expect_equal(class(overlap), "dist")
    overlap <- as.matrix(overlap)
    expect_equal(nrow(overlap), ncol(assay(tse)))
    expect_equal(nrow(overlap), ncol(overlap))
    expect_error(calculateOverlap(tse, detection = "test"))
    expect_error(calculateOverlap(tse, detection = TRUE))
    expect_error(calculateOverlap(tse, abund_values = "test"))
    esophagus <- transformSamples(esophagus, method = "relabundance")
    overlap <- as.matrix(calculateOverlap(esophagus, abund_values = "relabundance", detection = 0.05))
    dimnames(overlap) <- NULL
    # Reference values from microbiome::overlap(esophagus, detection = 0.05)
    col1 <- c(0.0000000, 0.6572394, 0.5434127)
    col2 <- c(0.6572394, 0.0000000, 0.4849852)
    col3 <- c(0.5434127, 0.4849852, 0.0000000)
    test <- matrix(c(col1, col2, col3), ncol = 3)
    expect_equal(round(overlap,7), round(test,7))
})
