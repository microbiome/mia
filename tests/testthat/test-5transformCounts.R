context("transformAssay")

# Small deterministic data with zeros. Small data keeps the tests fast; the
# numerical details of the transformations are tested in the packages that
# implement them (e.g. vegan).
.mat <- matrix(
    (0:39) %% 7, nrow = 8, ncol = 5,
    dimnames = list(paste0("f", 1:8), paste0("s", 1:5)))
.se <- SummarizedExperiment(assays = list(counts = .mat))

test_that("transformAssay input checks", {
    se <- .se
    expect_error(transformAssay(se))
    expect_error(transformAssay(se, method = "test"))
    expect_error(
        transformAssay(se, method = c("log10", "pa"), name = "test"),
        "'method' must be a non-empty single character value")
    expect_error(
        transformAssay(se, method = "relabundance", assay.type = "test"),
        "'assay.type' must be a valid name of assays")
    expect_error(
        transformAssay(se, method = "relabundance", name = FALSE),
        "'name' must be a non-empty single character value")
    expect_error(
        transformAssay(se, method = "relabundance", name = c("a", "b")))
    expect_error(transformAssay(se, method = "relabundance", MARGIN = 3))
    expect_error(
        transformAssay(se, method = "log10", pseudocount = "pseudocount"))
    expect_error(transformAssay(se, method = "log10", pseudocount = -1))
    expect_error(
        transformAssay(se, method = "alr", reference = c(1, 2)),
        "'reference' must be a single value")
    # Deprecated argument
    expect_warning(res <- transformAssay(
        as(se, "TreeSummarizedExperiment"), method = "pa",
        assay_name = "counts"))
    expect_true("pa" %in% assayNames(res))
})

test_that("transformAssay adds assay and keeps the input intact", {
    res <- transformAssay(.se, method = "relabundance")
    expect_named(assays(res), c("counts", "relabundance"))
    expect_equal(assay(res, "counts"), .mat)
    expect_identical(dimnames(assay(res, "relabundance")), dimnames(.mat))
    res <- transformAssay(.se, method = "pa", name = "test")
    expect_true("test" %in% assayNames(res))
    # Matrix without dimnames
    se <- SummarizedExperiment(assays = list(counts = unname(.mat)))
    res <- transformAssay(se, method = "relabundance")
    expect_equal(
        assay(res, "relabundance"), t(t(.mat) / colSums(.mat)),
        check.attributes = FALSE)
    expect_null(dimnames(assay(res, "relabundance")))
    # Internal helper used by other functions
    expect_equal(
        .calc_rel_abund(.mat), t(t(.mat) / colSums(.mat)),
        check.attributes = FALSE)
})

test_that("transformAssay equals vegan::decostand", {
    # mia method name = vegan method name
    methods <- c(
        relabundance = "total", total = "total", frequency = "frequency",
        hellinger = "hellinger", log_zero_replace = "log", max = "max",
        normalize = "normalize", pa = "pa", range = "range", rank = "rank",
        rrank = "rrank", standardize = "standardize",
        chi.square = "chi.square", rclr = "rclr")
    for( MARGIN in c(1L, 2L) ){
        for( method in names(methods) ){
            test <- transformAssay(.se, method = method, MARGIN = MARGIN)
            test <- assay(test, method)
            ref <- vegan::decostand(
                .mat, method = methods[[method]], MARGIN = MARGIN)
            # chi.square returns a transposed table
            if( identical(dim(ref), rev(dim(.mat))) ){
                ref <- t(ref)
            }
            expect_equal(
                unclass(test)[, ], unclass(ref)[, ],
                check.attributes = FALSE,
                info = paste0(method, ", MARGIN = ", MARGIN))
        }
    }
    # Additional arguments are passed to vegan
    test <- transformAssay(.se, method = "log_zero_replace", logbase = 10)
    ref <- vegan::decostand(.mat, method = "log", logbase = 10, MARGIN = 2)
    expect_equal(
        assay(test, "log_zero_replace"), ref, check.attributes = FALSE)
    # clr requires pseudocount with zeroes
    expect_error(transformAssay(.se, method = "clr"))
    test <- transformAssay(.se, method = "clr", pseudocount = 1)
    ref <- vegan::decostand(.mat, method = "clr", pseudocount = 1, MARGIN = 2)
    expect_equal(assay(test, "clr"), ref, check.attributes = FALSE)
    # alr: vegan drops the reference feature; mia adds it back and fills it
    # with 'reference'
    test <- transformAssay(.se, method = "alr", pseudocount = 1)
    test <- assay(test, "alr")
    ref <- vegan::decostand(.mat, method = "alr", pseudocount = 1, MARGIN = 2)
    expect_identical(dimnames(test), dimnames(.mat))
    expect_equal(
        test[rownames(ref), ], ref, check.attributes = FALSE)
    expect_true(all(is.na(test[setdiff(rownames(.mat), rownames(ref)), ])))
    test <- transformAssay(.se, method = "alr", pseudocount = 1, reference = 0)
    expect_true(all(assay(test, "alr")[1, ] == 0))
    # alr feature-wise; vegan drops the reference sample
    test <- assay(
        transformAssay(.se, method = "alr", pseudocount = 1, MARGIN = 1L),
        "alr")
    ref <- vegan::decostand(.mat, method = "alr", pseudocount = 1, MARGIN = 1)
    expect_identical(dimnames(test), dimnames(.mat))
    expect_equal(test[, colnames(ref)], ref, check.attributes = FALSE)
    expect_true(all(is.na(test[, setdiff(colnames(.mat), colnames(ref))])))
    # Deprecated alias
    expect_warning(
        test <- transformAssay(.se, method = "z", MARGIN = 1L, name = "z"))
    ref <- vegan::decostand(.mat, method = "standardize", MARGIN = 1)
    expect_equal(assay(test, "z"), ref, check.attributes = FALSE)
})

test_that("transformAssay pseudocount", {
    min_half <- min(.mat[.mat > 0]) / 2
    # TRUE adds half of the minimum positive value
    res <- transformAssay(.se, method = "pseudocount")
    expect_equal(assay(res, "pseudocount"), .mat + min_half,
        check.attributes = FALSE)
    res <- transformAssay(.se, method = "relabundance", pseudocount = TRUE)
    expect_equal(
        attr(assay(res, "relabundance"), "parameters")[["pseudocount"]],
        min_half)
    # Pseudocount is added before transformation
    ref <- transformAssay(.se, method = "relabundance", pseudocount = 2)
    se <- .se
    assay(se, "counts") <- .mat + 2
    test <- transformAssay(se, method = "relabundance")
    expect_equal(
        assay(test, "relabundance"), assay(ref, "relabundance"),
        check.attributes = FALSE)
    # FALSE equals 0
    a <- transformAssay(.se, method = "relabundance", pseudocount = FALSE)
    b <- transformAssay(.se, method = "relabundance", pseudocount = 0)
    expect_equal(assay(a, "relabundance"), assay(b, "relabundance"))
    # Negative values with TRUE
    se <- .se
    assay(se, "neg") <- .mat - 2
    expect_error(transformAssay(
        se, assay.type = "neg", method = "relabundance", pseudocount = TRUE))
    # NAs are ignored when calculating the pseudocount
    assay(se, "na") <- .mat
    assay(se, "na")[1, 2] <- NA
    expect_warning(res <- transformAssay(
        se, assay.type = "na", method = "pseudocount"))
    expect_equal(
        attr(assay(res, "pseudocount"), "parameters")[["pseudocount"]],
        min_half)
    # Only positive values: no pseudocount added with TRUE
    assay(se, "pos") <- .mat + 1
    expect_message(res <- transformAssay(
        se, assay.type = "pos", method = "relabundance", pseudocount = TRUE))
    expect_equal(
        attr(assay(res, "relabundance"), "parameters")[["pseudocount"]], 0)
    # ... but warning if added with numeric value
    expect_warning(transformAssay(
        se, assay.type = "pos", method = "relabundance", pseudocount = 1))
    # Warning if numeric pseudocount is added to negative values
    expect_warning(
        res <- transformAssay(
            se, assay.type = "neg", method = "pseudocount", pseudocount = 2),
        "negative values")
    expect_equal(assay(res, "pseudocount"), .mat, check.attributes = FALSE)
})

test_that("transformAssay log", {
    for( method in c("log", "log2", "log10") ){
        # Error with zeroes
        expect_error(transformAssay(.se, method = method), "zeroes")
        FUN <- match.fun(method)
        for( MARGIN in c(1L, 2L) ){
            test <- transformAssay(
                .se, method = method, MARGIN = MARGIN, pseudocount = 1)
            expect_equal(
                assay(test, method), FUN(.mat + 1), check.attributes = FALSE)
        }
    }
    # Error with negative values
    se <- .se
    assay(se, "neg") <- .mat - 2
    expect_error(
        transformAssay(se, assay.type = "neg", method = "log"),
        "negative values")
})

test_that("transformAssay css", {
    test <- assay(transformAssay(.se, method = "css", percentile = 0.5), "css")
    expect_identical(dim(test), dim(.mat))
    # Parameters affect the result
    test2 <- transformAssay(.se, method = "css", percentile = 0.75)
    expect_false(isTRUE(all.equal(test, assay(test2, "css"))))
    test2 <- transformAssay(.se, method = "css", percentile = 0.5,
        scaling = 2000)
    expect_equal(assay(test2, "css"), test * 2, check.attributes = FALSE)
    expect_length(.calc_scaling_factors(.mat, 0.5), ncol(.mat))
    expect_error(transformAssay(.se, method = "css", percentile = "a"))
    expect_error(transformAssay(.se, method = "css", threshold = 2))
    expect_error(transformAssay(.se, method = "css", scaling = "a"))
    # Percentile cannot be estimated if a sample has 1 or less features
    se <- .se
    assay(se, "counts")[-1, 1] <- 0
    expect_error(transformAssay(se, method = "css"), "1 or less found")
    # Compare to values from metagenomeSeq (CSS is implemented in mia)
    data(GlobalPatterns, package = "mia")
    expect_message(tse <- transformAssay(GlobalPatterns, method = "css"))
    ref <- matrix(
        c(1.358696, 317.846608, 4.076087, 1.474926), nrow = 2, byrow = TRUE)
    expect_equal(
        assay(tse, "css")[19:20, 17:18], ref,
        tolerance = 1e-6, check.attributes = FALSE)
})

test_that("transformAssay binning", {
    tse <- transformAssay(.se, method = "binning", nbins = 3)
    res <- assay(tse, "binning")
    expect_true(all(res >= 0 & res <= 3))
    expect_true(all(res[.mat == 0] == 0))
    expect_true(all(res[.mat != 0] > 0))
    # Fewer non-zero values than bins
    se <- SummarizedExperiment(assays = list(
        counts = matrix(c(10, 5, 0, 0), ncol = 1)))
    se <- transformAssay(se, method = "binning", nbins = 4)
    expect_equal(as.vector(assay(se, "binning")), c(4, 1, 0, 0))
    # Feature-wise
    se <- SummarizedExperiment(assays = list(
        counts = matrix(c(10, 5, 0, 20, 0, 10), nrow = 2, byrow = TRUE)))
    se <- transformAssay(se, method = "binning", nbins = 3, MARGIN = 1L)
    res <- assay(se, "binning")
    expect_equal(as.vector(res[1, ]), c(3, 1, 0))
    expect_equal(as.vector(res[2, ]), c(3, 0, 1))
    # Single feature
    se <- SummarizedExperiment(assays = list(
        counts = matrix(c(10, 5, 0), nrow = 1)))
    se <- transformAssay(se, method = "binning", nbins = 2)
    expect_equal(dim(assay(se, "binning")), c(1L, 3L))
    # Errors
    expect_error(transformAssay(.se, method = "binning", nbins = "a"))
    expect_error(transformAssay(.se, method = "binning", nbins = 0))
    se <- transformAssay(.se, method = "clr", pseudocount = 1)
    expect_error(transformAssay(se, method = "binning", assay.type = "clr"))
})

test_that("transformAssay cutoff", {
    expect_identical(
        .apply_cutoff(c(0, 1, 2, 3, 0.5), 1), c(NA, NA, 2, 3, NA))
    expect_identical(
        .apply_cutoff(c(0, 1, 2, 3, 0.5), 2, value = 999),
        c(999, 999, 999, 3, 999))
    expect_error(
        .apply_cutoff(c(1, 2, 3), 1, value = c(1, 2)),
        "'value' must be a single numeric value or NA")
    expect_error(
        transformAssay(.se, method = "cutoff", threshold = "a"),
        "'threshold' must be a single numeric value")
    test <- transformAssay(.se, method = "cutoff", threshold = 2, value = 0)
    expect_equal(
        assay(test, "cutoff"), ifelse(.mat <= 2, 0, .mat),
        check.attributes = FALSE)
})

test_that("transformAssay difference and division", {
    mat <- .mat[1:4, ]
    se <- SummarizedExperiment(assays = list(counts = mat))
    pairs <- combn(rownames(mat), 2)
    # Difference
    res <- transformAssay(se, method = "difference", MARGIN = 1L)
    expect_true("difference" %in% altExpNames(res))
    res <- assay(altExp(res, "difference"))
    expect_s4_class(res, "dgCMatrix")
    expect_equal(attr(res, "mia"), "difference")
    expect_equal(dim(res), c(choose(nrow(mat), 2), ncol(mat)))
    ref <- mat[pairs[1, ], ] - mat[pairs[2, ], ]
    expect_equal(
        as.matrix(res)[paste0(pairs[1, ], "-", pairs[2, ]), ], ref,
        check.attributes = FALSE)
    # Alias
    res2 <- transformAssay(se, method = "-", MARGIN = 1L, name = "difference")
    expect_equal(assay(altExp(res2, "difference")), res)
    # Division; zero division gives 0
    expect_warning(res <- transformAssay(se, method = "division", MARGIN = 1L))
    res <- assay(altExp(res, "division"))
    expect_equal(attr(res, "mia"), "division")
    ref <- mat[pairs[1, ], ] / mat[pairs[2, ], ]
    ref[!is.finite(ref)] <- 0
    expect_equal(
        as.matrix(res)[paste0(pairs[1, ], "-", pairs[2, ]), ], ref,
        check.attributes = FALSE)
    # Missing names are generated
    se <- SummarizedExperiment(assays = list(counts = unname(mat)))
    expect_warning(
        res <- transformAssay(se, method = "difference", MARGIN = 1L),
        "No names found")
    expect_equal(nrow(altExp(res, "difference")), choose(nrow(mat), 2))
    # Warning with many features
    expect_warning(
        .apply_transformation_difference_or_division(
            matrix(1, nrow = 1, ncol = 1001), "difference", 2L),
        "over 1000 features")
})

test_that("transformAssay invnorm", {
    mat <- .mat
    mat[2, 3] <- NA
    se <- SummarizedExperiment(assays = list(counts = mat))
    res <- assay(transformAssay(
        se, method = "invnorm", ties.method = "average", offset = 0.5),
        "invnorm")
    expect_identical(dimnames(res), dimnames(mat))
    expect_identical(is.na(res), is.na(mat))
    expect_identical(attr(res, "mia"), "invnorm")
    expect_identical(attr(res, "parameters")$ties.method, "average")
    expect_identical(attr(res, "parameters")$offset, 0.5)
    ref <- apply(mat, 2, function(x){
        ok <- !is.na(x)
        x[ok] <- qnorm((rank(x[ok]) - 0.5) / sum(ok))
        x
    })
    expect_equal(res, ref, check.attributes = FALSE)
    # ties.method changes result
    res2 <- transformAssay(se, method = "invnorm", ties.method = "max")
    expect_false(isTRUE(all.equal(assay(res2, "invnorm"), res)))
    expect_error(transformAssay(se, method = "invnorm", offset = 0.75))
    expect_error(transformAssay(se, method = "invnorm", offset = -0.01))
    expect_error(transformAssay(se, method = "invnorm", ties.method = "nope"))
})

test_that("transformAssay philr", {
    skip_if_not_installed("philr")
    mat <- .mat[1:5, 1:4]
    tree <- ape::read.tree(text = "((f1:1,f2:1):1,((f3:1,f4:1):1,f5:1):1);")
    tse <- TreeSummarizedExperiment(
        assays = list(counts = mat), rowTree = tree)
    # The default is columns, and colTree is not present
    expect_error(transformAssay(tse, method = "philr", pseudocount = 1))
    # There must not be 0s
    expect_error(transformAssay(tse, method = "philr", MARGIN = 1L))
    # Tree must match
    wrong <- ape::read.tree(text = "((a:1,b:1):1,(c:1,d:1):1);")
    expect_error(transformAssay(
        tse, method = "philr", MARGIN = 1L, pseudocount = 1, tree = wrong))
    expect_message(tse <- transformAssay(
        tse, method = "philr", pseudocount = 1, MARGIN = 1L))
    expect_true("philr" %in% altExpNames(tse))
    ref <- t(philr::philr(t(mat + 1), tree))
    expect_equal(
        assay(altExp(tse, "philr"), "philr"), ref, check.attributes = FALSE)
    expect_equal(colData(tse), colData(altExp(tse, "philr")))
    # return.all
    tse <- TreeSummarizedExperiment(
        assays = list(counts = mat), rowTree = tree)
    expect_error(transformAssay(
        tse, method = "philr", MARGIN = 1L, pseudocount = 1,
        return.all = "a"), "'return.all' must be TRUE or FALSE")
    res <- transformAssay(
        tse, method = "philr", MARGIN = 1L, pseudocount = 1,
        return.all = TRUE) |> expect_message()
    res <- assay(altExp(res, "philr"), "philr")
    expect_equal(res, ref, check.attributes = FALSE)
    expect_false(is.null(attr(res, "V")))
    # Tree must be provided for SE
    se <- SummarizedExperiment(assays = list(counts = mat))
    expect_error(
        transformAssay(se, method = "philr", MARGIN = 1L, pseudocount = 1),
        "'tree' must be provided")
    # Tree must have branch lengths
    no_len <- ape::read.tree(text = "((f1,f2),((f3,f4),f5));")
    expect_error(
        transformAssay(se, method = "philr", MARGIN = 1L, pseudocount = 1,
            tree = no_len),
        "does not have any branches")
    # node.label must match the features
    expect_error(
        transformAssay(se, method = "philr", MARGIN = 1L, pseudocount = 1,
            tree = tree, node.label = "f1"),
        "'node.label' must be NULL")
    # Features without node label are dropped
    expect_warning(
        res <- transformAssay(
            se, method = "philr", MARGIN = 1L, pseudocount = 1, tree = tree,
            node.label = c(paste0("f", 1:4), NA)),
        "subsetted") |> suppressMessages()
    expect_equal(ncol(altExp(res, "philr")), ncol(mat))
    ref <- t(philr::philr(t(mat[1:4, ] + 1), ape::keep.tip(tree, 1:4)))
    expect_equal(
        assay(altExp(res, "philr"), "philr"), ref, check.attributes = FALSE)
    # Features must have names
    expect_error(
        transformAssay(
            SummarizedExperiment(assays = list(counts = unname(mat))),
            method = "philr", MARGIN = 1L, pseudocount = 1, tree = tree),
        "'x' must have rownames")
    # Warning if tree is not binary
    multi <- ape::read.tree(text = "(f1:1,f2:1,f3:1,f4:1,f5:1);")
    expect_warning(
        try(transformAssay(se, method = "philr", MARGIN = 1L,
            pseudocount = 1, tree = multi), silent = TRUE) |>
            suppressMessages(),
        "rooted and binary")
    # Sample-wise with colTree: dimensions change, so only the transformed
    # data is returned
    tse <- TreeSummarizedExperiment(
        assays = list(counts = t(mat)), colTree = tree)
    expect_warning(
        res <- transformAssay(tse, method = "philr", pseudocount = 1),
        "columns of the transformed data do not match")
    expect_equal(
        assay(res, "philr"), philr::philr(t(mat) + 1, tree),
        check.attributes = FALSE)
    # Sample-wise with colTree; samples without node label are dropped
    w <- (res <- transformAssay(
        tse, method = "philr", pseudocount = 1, tree = tree,
        node.label = c(paste0("f", 1:4), NA))) |>
        capture_warnings()
    grepl("subsetted", w) |> any() |> expect_true()
    ref <- philr::philr(t(mat[1:4, ]) + 1, ape::keep.tip(tree, 1:4))
    expect_equal(assay(res, "philr"), ref, check.attributes = FALSE)
})

test_that("transformAssay works with DelayedArray", {
    se <- .se
    assay(se, "counts") <- DelayedArray(.mat)
    args <- list(
        list(method = "relabundance"), list(method = "clr", pseudocount = 1),
        list(method = "log10", pseudocount = 1), list(method = "rank"),
        list(method = "pa", MARGIN = 1L),
        list(method = "css", percentile = 0.5),
        list(method = "binning"), list(method = "invnorm"))
    for( a in args ){
        test <- do.call(transformAssay, c(list(se), a))
        ref <- do.call(transformAssay, c(list(.se), a))
        expect_equal(
            as.matrix(assay(test, a$method)), as.matrix(assay(ref, a$method)),
            check.attributes = FALSE, info = a$method)
    }
})

test_that("transformAssay altexp", {
    tse <- as(.se, "TreeSummarizedExperiment")
    altExp(tse, "sub1") <- tse[1:3, ]
    altExp(tse, "sub2") <- tse[4:6, ]
    expect_error(transformAssay(tse, altexp = "test", method = "pa"))
    expect_error(transformAssay(tse, altexp = 3, method = "pa"))
    expect_error(transformAssay(tse, altexp = TRUE, method = "pa"))
    expect_error(transformAssay(tse, altexp = NA, method = "pa"))
    res <- transformAssay(tse, altexp = "sub1", method = "relabundance")
    ref <- transformAssay(altExp(tse, "sub1"), method = "relabundance")
    expect_equal(
        assay(altExp(res, "sub1"), "relabundance"), assay(ref, "relabundance"))
    expect_false("relabundance" %in% assayNames(altExp(res, "sub2")))
    res <- transformAssay(tse, altexp = altExpNames(tse), method = "pa")
    expect_true(all(vapply(
        altExps(res), function(y) "pa" %in% assayNames(y), logical(1))))
    # Errors in altExps are reported
    assay(tse, "test") <- .mat
    expect_error(
        transformAssay(tse, assay.type = "test", altexp = "sub1",
            method = "pa"),
        "Transforming altExps")
})

test_that("transformAssay and rarefyAssay warn when overwriting an assay", {
    se <- transformAssay(.se, method = "relabundance")
    expect_warning(
        se <- transformAssay(se, method = "relabundance"),
        "already present in `assays`")
    # The assay used as input can be overwritten
    ref <- transformAssay(se, method = "relabundance", name = "test")
    expect_warning(
        res <- transformAssay(se, method = "relabundance", name = "counts"),
        "already present in `assays`")
    expect_equal(assay(res, "counts"), assay(ref, "test"))
    expect_warning(
        rarefyAssay(se, name = "counts", sample = 5, verbose = FALSE),
        "already present in `assays`")
})

test_that("Multiplicative replacement of zeros", {
    data(GlobalPatterns, package = "mia")
    tse <- agglomerateByRank(GlobalPatterns, rank = "Phylum")
    tse <- transformAssay(tse, method = "multiplicative_replacement")
    x <- assay(tse, "counts")
    y <- assay(tse, "multiplicative_replacement")
    # Samples sum to one, zeros take 0.65 of one read
    expect_equal(unname(colSums(y)), rep(1, ncol(y)))
    total <- colSums(x)
    expect_equal(y[x == 0], (matrix(0.65 / total, nrow(x), ncol(x),
        byrow = TRUE))[x == 0])
    # Ratios between non-zero parts are kept
    nz <- which(x[, 1] > 0)[1:2]
    expect_equal(y[nz[1], 1] / y[nz[2], 1], x[nz[1], 1] / x[nz[2], 1])
    # Proportions: half the smallest positive proportion of each sample
    tse <- transformAssay(tse, method = "relabundance")
    tse <- transformAssay(tse, assay.type = "relabundance",
        method = "multiplicative_replacement", name = "mr")
    p <- assay(tse, "relabundance")
    dl <- apply(p, 2, function(v) min(v[v > 0]) / 2)
    expect_equal(assay(tse, "mr")[p == 0], (matrix(0.65 * dl, nrow(p),
        ncol(p), byrow = TRUE))[p == 0])
    # The clr of the replaced proportions
    tse <- transformAssay(tse, assay.type = "multiplicative_replacement",
        method = "clr")
    expect_false(anyNA(assay(tse, "clr")))
    expect_error(transformAssay(tse, method = "multiplicative_replacement",
        fraction = 2))
    # Zeros that would replace the whole sample
    expect_error(transformAssay(GlobalPatterns[1:200, colSums(assay(
        GlobalPatterns, "counts")[1:200, ]) > 0],
        method = "multiplicative_replacement"), "sum to one")
    expect_error(transformAssay(tse, method = "multiplicative_replacement",
        detection = -1))
})
