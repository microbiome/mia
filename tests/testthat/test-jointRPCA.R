test_that("single-view .joint_rpca returns the expected structure", {
    set.seed(123)
    #features x samples
    X <- abs(matrix(rnorm(20 * 12), nrow = 20, ncol = 12))
    rownames(X) <- paste0("f", 1:20)
    colnames(X) <- paste0("s", 1:12)
    
    fit <- mia:::.joint_rpca(
        tables = list(assay1 = X),
        n.components = 3,
        max.iterations = 2,
        rclr.transform.tables = FALSE,  
        n.test.samples = 4
    )
    
    #top-level structure
    expect_type(fit, "list")
    expect_true(all(c("ord.res", "dist", "cv.stats", "rclr.tables") %in% names(fit)))
    
    #ordination structure
    OR <- fit$ord.res
    expect_true(all(c("eigvals", "samples", "features", "proportion.explained") %in% names(OR)))
    
    #samples: should include both train + projected test = all original samples
    expect_equal(nrow(OR$samples), ncol(X))
    expect_equal(colnames(OR$samples), paste0("PC", 1:3))
    expect_true(all(colnames(OR$samples) %in% names(OR$proportion.explained)))
    
    #features: now per-view list, even for single-view
    expect_true(is.list(OR$features))
    expect_true("assay1" %in% names(OR$features))
    
    F1 <- OR$features$assay1
    if (is.data.frame(F1)) F1 <- as.matrix(F1)
    expect_equal(nrow(F1), nrow(X))
    expect_equal(colnames(F1), paste0("PC", 1:3))
    
    #cv stats & dist
    expect_true(is.data.frame(fit$cv.stats))
    expect_true(all(c("mean_CV", "std_CV", "run", "iteration") %in% names(fit$cv.stats)))
    
    expect_s3_class(fit$dist, "DistanceMatrix")
    expect_true(is.matrix(fit$dist$data))
    expect_equal(nrow(fit$dist$data), ncol(X))
    expect_equal(colnames(fit$dist$data), colnames(X))
})

test_that("multi-view .joint_rpca preserves per-view feature loadings", {
    set.seed(42)
    #two views, same samples
    S <- paste0("s", 1:10)
    A <- abs(matrix(rnorm(25 * 10), 25, 10,
                    dimnames = list(paste0("a", 1:25), S)))
    B <- abs(matrix(rnorm(30 * 10), 30, 10,
                    dimnames = list(paste0("b", 1:30), S)))
    
    fit <- mia:::.joint_rpca(
        tables = list(MGX = A, MTX = B),
        n.components = 2,
        max.iterations = 2,
        rclr.transform.tables = FALSE,
        n.test.samples = 3
    )
    
    OR <- fit$ord.res
    expect_true(is.list(OR$features))
    expect_true(all(c("MGX", "MTX") %in% names(OR$features)))
    
    #check each view's loading matrix dimensions & rownames
    MGXv <- OR$features$MGX; MTXv <- OR$features$MTX
    if (is.data.frame(MGXv)) MGXv <- as.matrix(MGXv)
    if (is.data.frame(MTXv)) MTXv <- as.matrix(MTXv)
    
    expect_equal(nrow(MGXv), nrow(A))
    expect_equal(nrow(MTXv), nrow(B))
    expect_equal(colnames(MGXv), paste0("PC", 1:2))
    expect_equal(colnames(MTXv), paste0("PC", 1:2))
    expect_true(all(rownames(MGXv) %in% rownames(A)))
    expect_true(all(rownames(MTXv) %in% rownames(B)))
    
    #samples include all training + projected test
    expect_equal(nrow(OR$samples), length(S))
    expect_equal(colnames(OR$samples), paste0("PC", 1:2))
})

test_that("unshared samples are dropped with a warning and alignment is correct", {
    set.seed(7)
    S1 <- paste0("s", 1:8)
    S2 <- c(paste0("s", 3:10)) 
    A <- abs(matrix(rnorm(20 * length(S1)), 20, length(S1),
                    dimnames = list(paste0("a", 1:20), S1)))
    B <- abs(matrix(rnorm(15 * length(S2)), 15, length(S2),
                    dimnames = list(paste0("b", 1:15), S2)))
    
    expect_warning(
        fit <- mia:::.joint_rpca(
            tables = list(MGX = A, MTX = B),
            n.components = 2,
            max.iterations = 2,
            rclr.transform.tables = FALSE,
            n.test.samples = 2
        ),
        regexp = "Removing.*sample\\(s\\).*overlap",
        all = FALSE
    )
    
    OR <- fit$ord.res
    expect_equal(nrow(OR$samples), 6)
})

test_that("projection of new samples via .transform appends rows and keeps component names", {
    set.seed(99)
    S_train <- paste0("s", 1:8)
    featsA <- paste0("a", 1:18)
    featsB <- paste0("b", 1:22)
    
    A <- abs(matrix(rnorm(length(featsA) * length(S_train)), length(featsA), length(S_train),
                    dimnames = list(featsA, S_train)))
    B <- abs(matrix(rnorm(length(featsB) * length(S_train)), length(featsB), length(S_train),
                    dimnames = list(featsB, S_train)))
    
    fit <- mia:::.joint_rpca(
        tables = list(MGX = A, MTX = B),
        n.components = 3,
        max.iterations = 2,
        rclr.transform.tables = FALSE,
        n.test.samples = 3
    )
    
    OR <- fit$ord.res
    n_before <- nrow(OR$samples)
    
    #create new samples (same features, new sample IDs)
    S_new <- c("s9", "s10")
    A_new <- abs(matrix(rnorm(length(featsA) * length(S_new)), length(featsA), length(S_new),
                        dimnames = list(featsA, S_new)))
    B_new <- abs(matrix(rnorm(length(featsB) * length(S_new)), length(featsB), length(S_new),
                        dimnames = list(featsB, S_new)))
    
    #project new samples
    OR2 <- mia:::.transform(
        ordination = OR,
        tables = list(MGX = A_new, MTX = B_new),
        apply.rclr = FALSE
    )
    
    expect_true(is.list(OR2))
    expect_true(all(c("samples", "features", "eigvals", "proportion.explained") %in% names(OR2)))
    expect_equal(colnames(OR2$samples), colnames(OR$samples))          
    expect_true(all(S_new %in% rownames(OR2$samples)))                 
    expect_equal(nrow(OR2$samples), n_before + length(S_new))          
})

test_that("errors surface for duplicated sample IDs during preprocessing", {
    set.seed(101)
    X <- abs(matrix(rpois(24, lambda = 5), 6, 4))
    rownames(X) <- paste0("f", 1:6)
    colnames(X) <- c("s1", "s1", "s2", "s3")  
    
    expect_error(
        mia:::.joint_rpca(
            tables = list(assay1 = X),
            n.components = 2,
            max.iterations = 2,
            rclr.transform.tables = FALSE,
            n.test.samples = 2
        ),
        regexp = "duplicate sample \\(column\\) IDs",
        ignore.case = TRUE
    )
})

test_that("jointRPCAuniversal works on MultiAssayExperiment", {
    set.seed(2025)
    
    #two views, same samples
    S <- paste0("s", 1:8)
    MGX_mat <- abs(matrix(
        rnorm(15 * length(S)),
        nrow = 15, ncol = length(S),
        dimnames = list(paste0("mgx", 1:15), S)
    ))
    MTX_mat <- abs(matrix(
        rnorm(10 * length(S)),
        nrow = 10, ncol = length(S),
        dimnames = list(paste0("mtx", 1:10), S)
    ))
    
    se_mgx <- SummarizedExperiment::SummarizedExperiment(
        assays = list(counts = MGX_mat)
    )
    se_mtx <- SummarizedExperiment::SummarizedExperiment(
        assays = list(counts = MTX_mat)
    )
    
    mae2 <- MultiAssayExperiment::MultiAssayExperiment(
        experiments = list(MGX = se_mgx, MTX = se_mtx)
    )
    
    fit <- mia:::jointRPCAuniversal(
        x                     = mae2,
        experiments           = c("MGX", "MTX"),
        n.components          = 2,
        max.iterations        = 2,
        rclr.transform.tables = FALSE,
        n.test.samples        = 3
    )
    
    #basic structure
    expect_type(fit, "list")
    expect_true(all(c("ord.res", "dist", "cv.stats", "rclr.tables") %in% names(fit)))
    
    #rclr.tables should be named by experiment
    expect_true(is.list(fit$rclr.tables))
    expect_equal(names(fit$rclr.tables), c("MGX", "MTX"))
    
    #per-view feature loadings should also carry MGX/MTX names
    OR <- fit$ord.res
    expect_true(is.list(OR$features))
    expect_true(all(c("MGX", "MTX") %in% names(OR$features)))
    
    MGXv <- OR$features$MGX
    MTXv <- OR$features$MTX
    if (is.data.frame(MGXv)) MGXv <- as.matrix(MGXv)
    if (is.data.frame(MTXv)) MTXv <- as.matrix(MTXv)
    
    expect_equal(nrow(MGXv), nrow(MGX_mat))
    expect_equal(nrow(MTXv), nrow(MTX_mat))
    expect_equal(colnames(MGXv), paste0("PC", 1:2))
    expect_equal(colnames(MTXv), paste0("PC", 1:2))
})

test_that("jointRPCAuniversal uses default assays per experiment in a MAE and records them", {
    skip_if_not_installed("MultiAssayExperiment")
    skip_if_not_installed("SummarizedExperiment")
    
    S <- paste0("s", 1:6)
    A1_counts <- matrix(abs(rnorm(10 * 6)), nrow = 10,
                        dimnames = list(paste0("a", 1:10), S))
    A1_alt    <- A1_counts * 2  
    
    A2_counts <- matrix(abs(rnorm(8 * 6)), nrow = 8,
                        dimnames = list(paste0("b", 1:8), S))
    A2_alt    <- A2_counts * 3
    
    se1 <- SummarizedExperiment::SummarizedExperiment(
        assays = list(counts = A1_counts, alt = A1_alt)
    )
    se2 <- SummarizedExperiment::SummarizedExperiment(
        assays = list(counts = A2_counts, alt = A2_alt)
    )
    
    mae <- MultiAssayExperiment::MultiAssayExperiment(
        experiments = list(MGX = se1, MTX = se2)
    )
    
    fit <- mia:::jointRPCAuniversal(
        x = mae,
        experiments = c("MGX", "MTX"),
        n.components = 2,
        max.iterations = 2,
        rclr.transform.tables = FALSE,
        n.test.samples = 2
    )
    
    #basic structure
    expect_true(is.list(fit))
    expect_true(all(c("MGX", "MTX") %in% names(fit$rclr.tables)))
    
    expect_true(all(c("experiment.names", "assay.names.used") %in% names(fit)))
    expect_equal(fit$experiment.names, c("MGX", "MTX"))
    expect_named(fit$assay.names.used, c("MGX", "MTX"))
    expect_equal(fit$assay.names.used[["MGX"]], "counts")
    expect_equal(fit$assay.names.used[["MTX"]], "counts")
})