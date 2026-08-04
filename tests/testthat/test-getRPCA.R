
test_that(".calculate_pca returns valid PCA results", {
    set.seed(1)

    mat <- matrix(rnorm(50), nrow = 10)
    rownames(mat) <- paste0("S", seq_len(nrow(mat)))
    colnames(mat) <- paste0("F", seq_len(ncol(mat)))

    res <- .calculate_pca(mat, ncomponents = 3)

    expect_type(res, "list")
    expect_named(
        res,
        c("sample_scores", "varExplained", "rotation", "center")
    )

    expect_equal(dim(res$sample_scores), c(10, 3))
    expect_equal(dim(res$rotation), c(5, 3))
    expect_length(res$varExplained, 3)

    expect_equal(
        rownames(res$sample_scores),
        rownames(mat)
    )
    expect_equal(
        rownames(res$rotation),
        colnames(mat)
    )
})

test_that(".construct_rpca_result stores expected attributes", {

    pca <- list(
        sample_scores = matrix(1, 4, 2),
        varExplained = c(4, 2),
        rotation = matrix(1, 3, 2),
        center = c(1, 2, 3)
    )

    opt <- list(
        matrix = matrix(0, 4, 3),
        raw = list(
            X = matrix(1, 4, 2),
            Y = matrix(1, 3, 2),
            S = diag(2)
        )
    )

    d <- dist(matrix(rnorm(8), 4))

    res <- .construct_rpca_result(pca, opt, d)

    expect_equal(dim(res), c(4, 2))
    expect_true(is.numeric(attr(res, "percentVar")))
    expect_equal(
        sum(attr(res, "percentVar")),
        100,
        tolerance = 1e-10
    )

    expect_true(inherits(attr(res, "distance"), "dist"))
    expect_equal(attr(res, "lower_dim"), opt$matrix)
})

test_that(".determine_test_set_for_rpca validates arguments", {

    mat <- matrix(rnorm(100), 20)
    rownames(mat) <- seq_len(nrow(mat))

    expect_error(
        .determine_test_set_for_rpca(
            mat,
            n.test.samples = -1
        ),
        "n.test.samples"
    )

    expect_error(
        .determine_test_set_for_rpca(
            mat,
            test.ratio = 2
        ),
        "test.ratio"
    )

    expect_error(
        .determine_test_set_for_rpca(
            mat,
            test.ratio = 0
        ),
        "test.ratio"
    )
})

test_that(".determine_test_set_for_rpca returns requested number of samples", {

    set.seed(1)

    mat <- matrix(rnorm(100), 20)
    rownames(mat) <- paste0("S", seq_len(nrow(mat)))

    idx <- .determine_test_set_for_rpca(
        mat,
        n.test.samples = 5
    )

    expect_length(idx, 5)
    expect_true(all(idx %in% seq_len(nrow(mat))))
})

test_that(".get_lower_rank_mat preserves dimensions and dimnames", {

    set.seed(1)

    mat <- matrix(rnorm(50), 10)
    rownames(mat) <- paste0("S", 1:10)
    colnames(mat) <- paste0("F", 1:5)

    res <- .get_lower_rank_mat(mat)

    expect_named(res, c("matrix", "raw"))

    expect_equal(dim(res$matrix), dim(mat))
    expect_equal(rownames(res$matrix), rownames(mat))
    expect_equal(colnames(res$matrix), colnames(mat))
})

test_that(".project_test_set_to_rpca returns projected coordinates", {

    set.seed(1)

    pca <- matrix(rnorm(10), 5, 2)
    attr(pca, "rotation") <- matrix(rnorm(8), 4, 2)
    attr(pca, "varExplained") <- c(3, 2)

    mat <- matrix(rnorm(20), 5, 4)

    proj <- .project_test_set_to_rpca(pca, mat)

    expect_equal(dim(proj), c(5, 2))
    expect_false(any(is.na(proj)))
})

test_that(".project_test_set_to_rpca handles missing values", {

    set.seed(1)

    pca <- matrix(rnorm(10), 5, 2)
    attr(pca, "rotation") <- matrix(rnorm(8), 4, 2)
    attr(pca, "varExplained") <- c(3, 2)

    mat <- matrix(rnorm(20), 5, 4)
    mat[1, 1] <- NA

    proj <- .project_test_set_to_rpca(pca, mat)

    expect_false(any(is.na(proj)))
})

test_that(".run_joint_rpca_analysis validates test.set", {

    expect_error(
        .run_joint_rpca_analysis(
            list(matrix(1, 2, 2)),
            test.set = 1
        ),
        "test.set"
    )
})

test_that(".run_joint_rpca_analysis returns JointRPCA object", {

    set.seed(1)

    mat1 <- matrix(rnorm(120), 20, 6)
    mat2 <- matrix(rnorm(80), 20, 4)

    rownames(mat1) <- rownames(mat2) <- paste0("S", 1:20)

    res <- .run_joint_rpca_analysis(
        list(
            microbiota = mat1,
            metabolites = mat2
        ),
        ncomponents = 2
    )

    expect_s3_class(res, "JointRPCA")
    expect_equal(nrow(res), 20)
    expect_equal(ncol(res), 2)

    expect_true(!is.null(attr(res, "rotation")))
    expect_true(!is.null(attr(res, "percentVar")))
    expect_true(!is.null(attr(res, "reconstruct_error")))
    expect_true(!is.null(attr(res, "n_features")))
})

test_that("invalid inputs fail", {

    expect_error(
        .joint_optspace(matrix(1)),
        "'x'"
    )

    expect_error(
        .joint_optspace(
            list(
                matrix(1,2,2),
                matrix(1,3,2)
            )
        ),
        "equal number of samples"
    )

    expect_error(
        .joint_optspace(
            list(matrix(Inf,2,2))
        ),
        "Infinite"
    )

    expect_error(
        .joint_optspace(
            list(matrix(1,5,5)),
            ropt = 10
        ),
        "ropt"
    )
})

test_that("cv error is finite", {

    x <- matrix(rnorm(50), 10, 5)

    res <- .calculate_optspace_cv_error(
        x,
        diag(3),
        matrix(rnorm(15),5,3)
    )

    expect_true(is.finite(res))
    expect_gt(res, 0)
})

test_that("double centering produces zero row and column means", {

    x <- matrix(rnorm(100), 10)

    y <- .apply_double_centering(x)

    expect_equal(
        rowMeans(y),
        rep(0, nrow(y)),
        tolerance = 1e-10
    )

    expect_equal(
        colMeans(y),
        rep(0, ncol(y)),
        tolerance = 1e-10
    )
})

test_that("reconstruction error is non-negative", {

    x <- list(
        matrix(rnorm(50),10,5),
        matrix(rnorm(40),10,4)
    )

    res <- .run_joint_rpca_analysis(x)

    err <- attr(res, "reconstruct_error")

    expect_true(all(err >= 0))
    expect_true(all(is.finite(err)))
})

test_that("explained variance sums to 100%", {

    x <- list(
        matrix(rnorm(50),10,5),
        matrix(rnorm(40),10,4)
    )

    res <- .run_joint_rpca_analysis(x)

    expect_equal(
        sum(attr(res, "percentVar")),
        100,
        tolerance = 1e-8
    )
})

test_that("algorithm handles missing values", {

    set.seed(1)

    x <- list(
        matrix(rnorm(50),10,5),
        matrix(rnorm(40),10,4)
    )

    x[[1]][1,1] <- NA
    x[[2]][3,2] <- NA

    expect_no_error(
        .run_joint_rpca_analysis(x)
    )
})

