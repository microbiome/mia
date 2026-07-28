#' @name
#' getRPCA
#'
#' @title
#' Run (joint) robust principal component analysis (RPCA)
#'
#' @description
#' These functions implement robust principal component analysis (RPCA)
#' for single tables and joint RPCA for multiple tables.
#' \code{*RPCA} functions run RPCA for single table. \code{*JointRPCA} runs the
#' analysis jointly for multiple tables. \code{get*} functions return the
#' results of the analysis while \code{add*} adds them to the input object.
#'
#' @details
#' These functions perform robust principal component analysis (RPCA) using a
#' low-rank matrix approximation followed by principal component analysis.
#' Missing values are handled using matrix completion based on the OptSpace
#' algorithm.
#'
#' \strong{Single-table RPCA}
#'
#' For a single assay, the workflow is:
#' \enumerate{
#'   \item Extract the selected assay matrix
#'   \item Estimate a low-rank approximation of the matrix using the OptSpace
#'   algorithm, which reconstructs missing values and reduces noise.
#'   \item Apply double centering (row and column centering).
#'   \item Apply PCA (via singular value decomposition) to the centered matrix.
#'   \item Return the sample coordinates in PCA space together with additional
#'   information such as loadings, explained variance, and pairwise distances.
#' }
#'
#' \strong{Joint RPCA}
#'
#' When multiple assays are provided, joint RPCA estimates a shared low-rank
#' representation across all tables. Each table may contain different features,
#' but they must share the same samples.
#'
#' The workflow is:
#' \enumerate{
#'   \item Extract the selected experiments and assays.
#'   \item Estimate a joint low-rank representation using a joint OptSpace
#'   algorithm that learns shared sample factors while allowing
#'   table-specific feature loadings.
#'   \item Apply double centering (row and column centering).
#'   \item Perform PCA on the centered matrix. to obtain the final
#'   sample coordinates.
#'   \item Return the sample coordinates in PCA space together with additional
#'   information such as loadings, explained variance, and pairwise distances.
#' }
#'
#' To assess generalization, the joint RPCA procedure optionally splits the
#' samples into training and test sets. The low-rank model is learned using the
#' training data and the test samples are projected into the resulting PCA
#' space. Reconstruction error for the test set is returned as a measure of
#' model fit.
#'
#' The returned object contains the PCA sample scores as the main result,
#' with additional information (e.g., rotation matrix, explained variance,
#' reconstructed matrix, and reconstruction error) stored in attributes.
#'
#'
#' @return
#' \code{matrix}, \code{TreeSummarizedExperiment} or \code{MultiAssayExperiment}
#' object.
#'
#' @inheritParams addAlpha
#'
#' @param assay.type \code{Character scalar}. Specifies the name of assay
#' used in calculation. (Default: \code{"counts"})
#'
#' @param experiments \code{Character vector} or \code{integer scalar}. Names
#' or indices of experiments selected from \code{x}.
#'
#' @param assay.types \code{Character vector}. Names assays selected from
#' \code{experiments}.
#'
#'@param name \code{Character scalar}. Name of results stored in \code{x}.
#'(Default: \code{"RPCA"} or \code{"JointRPCA"} depending on the method)
#'
#' @param ... additional arguments:
#' \itemize{
#'   \item \code{ncomponents}:\code{Integer scalar}. The number of components
#'   to estimate. (Default: \code{3})
#'   \item \code{max.iterations}:\code{Integer scalar}. The number of iterations
#'   run in OptSpace algorithm. (Default: \code{3})
#'   \item \code{tolerance}:\code{Numeric scalar}. Accepted error between
#'   lower rank representation obtained by OptSpace algorithm and original
#'   table. (Default: \code{3})
#' }
#'
#' @examples
#' data("ibdmdb")
#' mae <- ibdmdb
#'
#' # Apply filtering
#' mae[[1]] <- filterRPCAInput(mae[[1]], assay.type = "mgx")
#' mae[[2]] <- filterRPCAInput(mae[[2]], assay.type = "mtx")
#'
#' # Apply data transformations. With impute=FALSE, missing values are preserved
#' # and not imputed.
#' mae[[1]] <- transformAssay(
#'     mae[[1]], assay.type = "mgx", method = "rclr", impute = FALSE)
#' mae[[2]] <- transformAssay(
#'     mae[[2]], assay.type = "mtx", method = "rclr", impute = FALSE)
#'
#' # Run joint-RPCA
#' res <- getJointRPCA(
#'     mae,
#'     experiments = c(1, 2),
#'     assay.types = c("rclr", "rclr")
#' )
#'
#' # Run RPCA for single experiment
#' res <- getRPCA(mae[[1]], assay.type = "rclr")
#'
#' @seealso
#' \code{\link[scater:runPCA]{scater::runPCA}}
#'
#' @author
#' The RPCA method is reported in Martino et al. (2020) and the
#' R/Bioconductor implementation utilizes the robust Aitchison
#' distance from \code{\link[vegan:decostand]{vegan::decostand}}.
#'
#' The Joint-RPCA method was adapted from the original
#' Python-based implementation in biocore/Gemelli by
#' Bianca Cordazzo Vargas, Liat Shenhav, and Cameron Martino.
#' The R/Bioconductor implementation was subsequently prepared by
#' Aituar Bektanov, Sabuj Bhowmick, Tuomas Borman, and Leo Lahti.
#'
#' @references
#'
#' Martino, C. and Shenhav, L. et al. (2020)
#' Context-aware dimensionality reduction deconvolutes gut microbial community
#' dynamics.
#' _Nat. Biotechnol._ doi:10.1038/s41587-020-0660-7
#'
NULL

#' @rdname getRPCA
#' @export
setMethod("getRPCA", signature = c(x = "SingleCellExperiment"),
    function(x, ...){
        x <- .check_and_get_altExp(x, ...)
        res <- callNextMethod(x, ...)
        return(res)
    }
)

#' @rdname getRPCA
#' @export
setMethod("getRPCA", signature = c(x = "SummarizedExperiment"),
    function(x, assay.type = "counts", ...){
        .check_assay_present(assay.type, x)
        mat <- assay(x, assay.type) |> t()
        res <- .calculate_rpca(mat, ...)
        return(res)
    }
)

#' @rdname getRPCA
#' @export
setMethod("addRPCA", signature = c(x = "SummarizedExperiment"),
    function(x, name = "RPCA", ...){
        if( !.is_a_string(name) ){
            stop("'name' must be a single character value.", call. = FALSE)
        }
        # Hiddenly support altExp
        x <- .check_and_get_altExp(x, ...)
        # Calculate indices
        args <- c(list(x = x), list(...))
        args <- args[ !names(args) %in% c("altexp") ]
        res <- do.call(getRPCA, args)
        # Add object to reducedDim
        x <- .add_object_to_reduceddim(x, res, name = name, ...)
        return(x)
    }
)

#' @rdname getRPCA

#' @export
setMethod("getJointRPCA", signature = c(x = "MultiAssayExperiment"),
    function(x, experiments, assay.types, ...){
        if( !(length(experiments) == length(assay.types) &&
                length(experiments) > 1L) ){
            stop("The lengths of 'experiments' and 'assay.types' must match ",
                "and there must be multiple experiments selected.",
                call. = FALSE)
        }
        mat_list <- .prepare_mae_for_joint_rpca(x, experiments, assay.types)
        res <- .run_joint_rpca_analysis(mat_list, ...)
        return(res)
    }
)

#' @rdname getRPCA
#' @export
setMethod("getJointRPCA", signature = c(x = "SingleCellExperiment"),
    function(x, experiments, assay.types, ...){
        if( !(length(experiments) == length(assay.types) &&
                length(experiments) > 1L) ){
            stop("The lengths of 'experiments' and 'assay.types' must match ",
                "and there must be multiple experiments selected.",
                call. = FALSE)
        }
        mat_list <- .prepare_tse_for_joint_rpca(x, experiments, assay.types)
        res <- .run_joint_rpca_analysis(mat_list, ...)
        return(res)
    }
)

#' @rdname getRPCA
#' @export
setMethod("addJointRPCA", signature = c(x = "SingleCellExperiment"),
    function(x, name = "JointRPCA", ...){
        if( !.is_a_string(name) ){
            stop("'name' must be a single character value.", call. = FALSE)
        }
        res <- getJointRPCA(x, ...)
        x <- .add_object_to_reduceddim(x, res, name = name, ...)
        return(x)
    }
)

#' @rdname getRPCA
#' @export
setMethod("addJointRPCA", signature = c(x = "MultiAssayExperiment"),
    function(x, name = "JointRPCA", ...){
        if( !.is_a_string(name) ){
            stop("'name' must be a single character value.", call. = FALSE)
        }
        res <- getJointRPCA(x, ...)
        x <- .add_values_to_mae_metadata(x, res, names = name, ...)
        return(x)
    }
)

################################ HELP FUNCTIONS ################################

# This function retrieves specific tables from MAE
#' @importFrom MultiAssayExperiment intersectColumns
.prepare_mae_for_joint_rpca <- function(x, experiments, assay.types, ...){
    # Select experiments from MAE
    x <- .select_experiments(x, experiments)
    # Select samples that are shared between experiments
    x <- intersectColumns(x)
    # Select alternative experiments
    x <- .select_altexps(x, ...)
    # Select assays of each experiment
    x <- .select_assays_from_mae(x, assay.types)
    # Get tables as a list
    mat_list <- MultiAssayExperiment::assays(x)
    # Change orientation so that samples are in rows
    mat_list <- lapply(mat_list, t)
    return(mat_list)
}

# This function retrieves specific tables from TreeSE
.prepare_tse_for_joint_rpca <- function(x, experiments, assay.types){
    # Select experiments from TreeSE
    tse_list <- .select_experiments_from_tse(x, experiments)
    # Select assays of each experiment
    tse_list <- .select_assays_from_tse_list(tse_list, assay.types)
    # Get tables as a list
    mat_list <- lapply(tse_list, assay)
    # Change orientation so that samples are in rows
    mat_list <- lapply(mat_list, t)
    return(mat_list)
}

# This function calculates Joint-RPCA for multiple tables. It runs the whole
# analysis from test/train set split to projecting the results to test set.
.run_joint_rpca_analysis <- function(mat_list, test.set = NULL, ...){
    if( !(is.null(test.set) || is.character(test.set)) ){
        stop("'test.set' must specify sample names for test set or be NULL.",
            call. = FALSE)
    }
    # Determine train/test split. User can define test set samples with vector
    # or then we can select representative samples based on RPCA of first table.
    all_samples <- mat_list[[1L]] |> rownames()
    if( !is.null(test.set) ){
        test_samples <- which( all_samples %in% test.set )
    } else{
        test_samples <- .determine_test_set_for_rpca(mat_list[[1L]], ...)
    }
    all_index <- all_samples |> length() |> seq_len()
    train_samples <- all_index[ !all_index %in% test_samples ]

    # Split tables to train and test sets
    train_set <- lapply(mat_list, function(x){
        x[train_samples, , drop = FALSE]
    })
    test_set <- lapply(mat_list, function(x){
        x[test_samples, , drop = FALSE]
    })
    # If user did not specify test set, use training set as test set
    if( all(lengths(test_set) == 0L) ){
        test_set <- train_set
    }

    # Calculate RPCA
    res <- .calculate_joint_rpca(train_set, test_set, ...)
    # Project test samples to PCA space determined by train set
    test_mat <- do.call(cbind, test_set)
    projected <- .project_test_set_to_rpca(res, test_mat)
    # Add test samples to pca results
    attr_list <- attributes(res)
    res <- rbind(res, projected)
    # Sort back to original order
    res <- res[
        order(c(train_samples, test_samples)), , drop = FALSE]
    # Add additional info back
    attr_list <- c(attributes(res), attr_list)
    attr_list <- attr_list[ !duplicated(names(attr_list)) ]
    attributes(res) <- attr_list

    # Calculate error between low rank and original test set
    num_features <- vapply(mat_list, ncol, numeric(1L))
    reconstruct_error <- .calculate_reconstruct_error(
        res, test_set, num_features)
    attributes(res)[["reconstruct_error"]] <- reconstruct_error

    return(res)
}

# This function runs RPCA to single table
#' @importFrom stats dist
.calculate_rpca <- function(mat, ncomponents = 3L, ...){
    # Get lower rank representation of the data
    opt_results <- .get_lower_rank_mat(mat, ncomponents = ncomponents, ...)
    # The result might have lower number of columns if they were not able to be
    # estimated.
    ncomponents <- opt_results[["raw"]][["S"]] |> ncol()
    # Apply pca to lower rank representation
    pca_results <- .calculate_pca(
        opt_results[["matrix"]], ncomponents = ncomponents)
    # Calculate distance in PCA space
    distance <- pca_results[["sample_scores"]] |> dist()
    # Create a final results to return to user
    res <- .construct_rpca_result(pca_results, opt_results, distance)
    return(res)
}

# This function runs Joint-RPCA. The only difference to .calculate_rpca is that
# the lower dimension matrix is estimated by optimizing feature loadings
# separately (sample loadings and singular values are estimated jointly).
#' @importFrom stats dist
.calculate_joint_rpca <- function(mat_list, test_set, ncomponents = 3L, ...){
    # Get lower rank representation of the data
    opt_results <- .get_lower_rank_joint_mat(
        mat_list, test_set, ncomponents = ncomponents, ...)
    # The result might have lower number of columns if they were not able to be
    # estimated.
    ncomponents <- opt_results[["raw"]][["S"]] |> ncol()
    # Apply pca to lower rank representation
    pca_results <- .calculate_pca(
        opt_results[["matrix"]], ncomponents = ncomponents)
    # Calculate distance in PCA space
    distance <- pca_results[["sample_scores"]] |> dist()
    # Create a final results to return to user
    pca_result <- .construct_rpca_result(pca_results, opt_results, distance)
    return(pca_result)
}

# This function constructs a lower rank representation from the data. The idea
# is to extract the essential from the data and to remove noise.
#' @importFrom vegan optspace
.get_lower_rank_mat <- function(
        mat,
        ncomponents = pmin(3L, nrow(mat), ncol(mat)),
        max.iterations = 5L,
        tolerance = 1e-5,
        ...){
    # Create lower rank representation
    opt_result <- optspace(
        x = mat,
        ropt = ncomponents,
        niter = max.iterations,
        tol = tolerance,
        verbose = FALSE
    )

    # Reconstruct the matrix
    X_hat <- opt_result[["X"]] %*% opt_result[["S"]] %*% t(opt_result[["Y"]])
    # Add old feature and sample names as they are dropped off
    dimnames(X_hat) <- dimnames(mat)

    # Create a result list
    res <- list(
        matrix = X_hat,
        raw = opt_result
    )

    return(res)
}

# Construct lower rank representation from a list of matrices.
.get_lower_rank_joint_mat <- function(
        mat_list, mat_list_test,
        ncomponents = pmin(3L, nrow(mat_list[[1L]]), ncol(mat_list[[1L]])),
        max.iterations = 5L,
        ...){
    # Create lower rank representation
    opt_result <- .joint_optspace(
        x = mat_list,
        x_test = mat_list_test,
        ropt = ncomponents,
        niter = max.iterations
    )

    # Reconstruct the matrix
    X_hat <- opt_result[["X"]] %*% opt_result[["S"]] %*% t(opt_result[["Y"]])
    # Add old feature and sample names as they are dropped off
    nams <- do.call(cbind, mat_list) |> dimnames()
    dimnames(X_hat) <- nams

    # Create a result list
    res <- list(
        matrix = X_hat,
        raw = opt_result
    )

    return(res)
}

# This function applies PCA to the data.
.calculate_pca <- function(mat, ncomponents){
    # Double center the data
    mat <- .apply_double_centering(mat)

    # Run PCA
    svd_result <- mat |> svd()
    u <- svd_result[["u"]]
    s <- svd_result[["d"]]
    v <- svd_result[["v"]]

    # We multiply U by singular values so that the sample coordinates reflect
    # actual variance magnitude rather than just orthonormal directions.
    # u <- u %*% diag(s)

    # Subset. There might be more components than requested.
    u <- u[ , seq_len(ncomponents), drop = FALSE]
    s <- s[ seq_len(ncomponents) ]
    v <- v[ , seq_len(ncomponents), drop = FALSE]

    # Adjust dimnames
    names(s) <- paste0("PC", s |> length() |> seq_len())
    dimnames(u) <- list(rownames(mat), names(s))
    dimnames(v) <- list(colnames(mat), names(s))

    # Create a result list
    pca_results <- list(
        sample_scores = u,
        varExplained = s,
        rotation = v,
        center = attributes(mat)[["center"]]
    )

    return(pca_results)
}

# This function constructs a final result to user. It does not calculate, but
# re-structures the results to returned format.
.construct_rpca_result <- function(pca_results, opt_results, distance){
    # We return the PCA sample scores as main results
    mat <- pca_results[["sample_scores"]]
    attr_list <- pca_results[ c("varExplained", "rotation", "center", "scale") ]

    # Calculate the explained variance in percentages
    percent_var <- pca_results[["varExplained"]]^2 /
        sum(pca_results[["varExplained"]]^2) * 100
    attr_list[["percentVar"]] <- percent_var

    # Add the lower rank representation and distance to result list
    attr_list[["lower_dim"]] <- opt_results[["matrix"]]
    attr_list <- c(attr_list, opt_results[["raw"]])
    attr_list[["distance"]] <- distance

    # Add additional results to main result matrix
    attr_list <- c(attributes(mat), attr_list)
    attributes(mat) <- attr_list

    return(mat)
}

# This function selects a representative set of samples to test set. This is
# done by applying RPCA for the first table and selecting samples from PC1 with
# a highest variance.
.determine_test_set_for_rpca <- function(
        mat, n.test.samples = NULL, test.ratio = 0.2, ...){
    if( !(is.null(n.test.samples) ||
            (.is_an_integer(n.test.samples) && n.test.samples > 0)) ){
        stop("'n.test.samples' must be a single positive integer value.",
            call. = FALSE)
    }
    if( !(.is_a_numeric(test.ratio) && test.ratio > 0 && test.ratio < 1) ){
        stop("'test.ratio' must be a numeric value in the range [0, 1].",
            call. = FALSE)
    }
    # Select number of samples
    if( is.null(n.test.samples) ){
        n.test.samples <- ceiling(test.ratio * nrow(mat))
    }
    test_samples <- c()
    if( n.test.samples > 0L ){
        # Calculate RPCA
        pca_result <- .calculate_rpca(mat, ...)
        # Select N samples so that they span over the PC1 axis. Idea is that
        # these samples should represent the dataset the best as there are
        # maximally different samples.
        first_component <- pca_result[, 1] |> sort()
        test_samples <- seq(
            1, length(first_component), length.out = n.test.samples) |>
            round()
        test_samples <- names(first_component)[ test_samples ]
        test_samples <- match(test_samples, rownames(pca_result))
    }
    return(test_samples)
}

# This function projects test set samples to PCA space that were obtained with
# train set.
.project_test_set_to_rpca <- function(pca_result, mat){
    # Extract PCA components
    feature_scores <- attributes(pca_result)[["rotation"]]
    singular_values <- attributes(pca_result)[["varExplained"]]

    # Apply double-centering
    mat <- mat |> .apply_double_centering()

    # NAs are set to zero during matrix multiplication so they do not contribute
    # to the projection. Otherwise, the NAs propagate and the projection step
    # fails.
    mat[ is.na(mat) ] <- 0

    # Project into PCA space
    projected <- mat %*% feature_scores

    # Normalize based on singular values
    projected <- projected / sqrt(sum(singular_values^2))

    return(projected)
}

# This function calculates error between raw test set values and their
# lower rank representation. The lower rank representation is calculated based
# on parameters learned from train set. The idea is to assess, how well the
# lower rank representation learns the generic, generalizable patterns from the
# data.
.calculate_reconstruct_error <- function(res, test_set, num_features){
    # Get learned parameters
    u_shared <- attributes(res)[["X"]]
    s_shared <- attributes(res)[["S"]]
    y_shared <- attributes(res)[["Y"]]

    # Split feature loadings by table
    ends   <- cumsum(num_features)
    starts <- c(1, head(ends, -1) + 1)
    y_individual <- mapply(function(s, e) {
        y_shared[s:e, ]
    }, starts, ends)

    # Calculate error separately for each table
    errors_per_set <- vapply(seq_len(length(test_set)), function(i){
        test_mat <- test_mat_zeroed <- test_set[[i]]
        # NAs are set to zero during matrix multiplication so they do not
        # contribute to the projection.
        test_mat_zeroed[ is.na(test_mat_zeroed) ] <- 0
        # Create projection for test set
        u_test <- test_mat_zeroed %*% y_individual[[i]]
        u_test <- sweep(u_test, 2, diag(s_shared), "/")
        recon_test <- u_test %*% s_shared %*% t(y_individual[[i]])

        # Calculate error between actual values and lower rank representation
        # Note: we don't use "test_mat_zeroed" here since we only compute the
        # error on observed entries (these are non-NA, non-zero entries)
        error <- test_mat - recon_test
        error[is.na(error)] <- 0
        error <- norm(error, "F") / sqrt(sum(!is.na(test_mat)))
        return(error)
    }, numeric(1L))
    names(errors_per_set) <- names(num_features)

    return(errors_per_set)
}

# -----------------------------------------------------------------------------
# .joint_optspace
#
# Run the *joint OptSpace* matrix completion algorithm on multiple tables that
# share the same rows (samples) but may contain different feature sets.
#
# Goal
# ----
# Estimate a shared low-rank representation across several matrices with
# missing values. Each matrix shares the same sample space (rows) but has
# its own features (columns). The algorithm learns:
#
#   X (U_shared)        : shared sample scores
#   S (S_shared)        : singular values (component scaling)
#   Y (V_individual)    : feature loadings for each matrix (table-specific)
#
# Model approximation:
#
#   M_i ≈ U_shared %*% S_shared %*% t(V_i)
#
# where M_i is table i and V_i are feature loadings specific to that table.
#
# References used for implementation guidance:
#   - gemelli implementation: https://github.com/biocore/gemelli/blob/00e3993f7358006d99a481ad281ee5801c6e159e/gemelli/optspace.py#L200
#   - vegan::optspace: https://github.com/vegandevs/vegan/blob/5f992e14f1d644d1734195a30042f18028d090b1/R/optspace.R#L12
#
# Returns
# -------
# list(
#   X        = sample scores (U_shared)
#   S        = singular value matrix
#   Y        = stacked feature loadings
#   cv_error = cross-validation error
# )
# -----------------------------------------------------------------------------
#
.joint_optspace <- function(x, x_test = x, ropt = 3, niter = 5){
    # Validate input
    if( !(is.list(x) && all(vapply(
            x, function(mat) is.matrix(mat) || is.data.frame(mat),
            logical(1L))) ) ){
        stop("'x' must be a list of matrices.", call. = FALSE)
    }
    if( vapply(x, nrow, integer(1L)) |> unique() |> length() != 1L ){
        stop("All tables must have equal number of samples (rows).",
            call. = FALSE)
    }
    if( do.call(cbind, x) |> is.infinite() |> any() ){
        stop("Infinite values are not allowed.", call. = FALSE)
    }
    if( !.is_an_integer(niter) ){
        stop("'niter' must be a single integer value.", call. = FALSE)
    }
    if( !.is_an_integer(ropt) ){
        stop("'ropt' must be a single integer value.", call. = FALSE)
    }
    #
    # Input tables can be matrices also, so make sure that they are now
    # matrices.
    x <- lapply(x, as.matrix)
    x_test <- lapply(x_test, as.matrix)
    # Number of input tables
    n_tables <- x |> length()
    # Number of samples (rows). Guaranteed equal across tables.
    n_samples <- x[[1]] |> nrow()
    # Number of features per table
    n_features <- vapply(x, ncol, integer(1L))
    # Total number of features after stacking tables
    total_n_features <- n_features |> sum()
    # Smallest feature count among tables
    # Used for validating maximum rank.
    min_feat <- n_features |> min()
    # ropt cannot exceed number of samples or smallest feature dimension
    if( (ropt < 1) || (ropt > min_feat) || (ropt > n_samples) ){
        stop("'ropt' must be integer in [1, ",
            min(n_samples, min_feat),
            "]", call. = FALSE)
    }

    # Prepare tables for optspace. In first table convert all NA values to 0.
    observed_list <- lapply(x, function(mat){
        mat[ is.na(mat) ] <- 0
        return(mat)
    })
    # The second table shows which cells included a value. Zeroes are also
    # treated as missing.
    mask_list <- lapply(x, function(mat){
        mat[ is.na(mat) ] <- 0
        mask <- abs(mat) > 0
        storage.mode(mask) <- "integer"
        return(mask)
    })

    # Stack all tables column-wise into a single matrix. We will use these
    # tables in initialization and calculating
    observed_stacked <- do.call(cbind, observed_list)
    mask_stacked <- do.call(cbind, mask_list)

    # eps is used to rescale the observed matrix to account for missing entries.
    total_non_zeroes <- mask_stacked |> sum()
    eps <- total_non_zeroes / sqrt(total_n_features * n_samples)
    # rho is gradient scaling factor. Helps stabilize gradient descent when
    # matrices are sparse.
    rho <- eps * n_samples

    # When many entries are missing, the magnitude of observed entries tends
    # to underestimate the magnitude of the full matrix. We compensate
    # for this by rescaling the observed matrix before optimization.
    #
    # The scale factor is reversed at the end of the algorithm.
    rescale_param <- sum(mask_stacked) * ropt
    rescale_param <- sqrt(rescale_param / (norm(observed_stacked, "F")^2))
    observed_stacked <- rescale_param * observed_stacked

    # We make initial guess for U_shared, S_shared and V. SVD with stacked data
    # gives sensible starting point for gradient descent.
    init_res <- .initialize_joint_optspace(
        observed_stacked, observed_list, mask_stacked, mask_list, ropt, eps,
        n_samples, n_features)
    U_shared <- init_res[["U_shared"]]
    S_shared <- init_res[["S_shared"]]
    V_list <- init_res[["V_list"]]

    # Now we start iteratively refine U_shared, S_shared and V
    cv_errors <- data.frame(mean = numeric(0L), sd = numeric(0L))
    for( i in seq_len(niter) ){
        sample_loadings <- vector("list", n_tables)
        cv_iter <- c()

        # Iterate over tables
        for( table_i in seq_len(n_tables) ){
            # Perform gradient update for current table
            res <- .gradient_update_joint_optspace(
                table_i, observed_list, mask_list, U_shared,
                S_shared, V_list, rho)
            # Store table-specific updates
            # Proposal for shared U
            sample_loadings[[table_i]] <- res[["U_i"]]
            # Table specific, updated V
            V_list[[table_i]] <- res[["V_i"]]

            # Calculate cross-validation error for test set
            cv_error <- .calculate_optspace_cv_error(
                x_test[[table_i]],
                V = res[["V_i"]],
                S = res[["S_i"]]
            )
            cv_iter <- c(cv_iter, cv_error)
        }

        # Combine CV error. We will create a table of errors where each row is
        # single iteration.
        # NumPy in Python implementation uses population SD (ddof=0) by default,
        # while R’s sd() uses sample SD (ddof=1).
        cv_iter <- data.frame(
            mean_CV = mean(cv_iter), std_CV = .population_sd(cv_iter))
        cv_errors <- rbind(cv_errors, cv_iter)

        # Update the shared sample factors (U_shared)
        # - Each table proposes its own update of U (sample loadings)
        # - Take the average across all tables to get a consensus
        U_shared <- Reduce("+", sample_loadings) / n_tables

        # Update the shared singular values (S_shared)
        # - First compute the average covariance of the updated U's
        # - Perform SVD on X_U to extract the principal singular values
        # - Normalize by Frobenius norm
        X_U <- Reduce(
            "+", lapply(sample_loadings, function(u) u %*% t(u))) / n_tables
        svd_res <- svd(X_U)

        # svd() can return fewer than `ropt` singular values when X_U is
        # rank-deficient (common here: X_U has rank <= ropt by construction,
        # and small singular values may be dropped by some LAPACK builds,
        # notably on Windows). Pad with zeros to keep length == ropt.
        d <- svd_res[["d"]]
        if (length(d) < ropt) {
            d <- c(d, rep(0, ropt - length(d)))
        }
        d <- d[seq_len(ropt)]

        # R's svd() returns singular values in descending order;
        # switch to ascending to match the Gemelli reference
        # implementation. Use nrow = ropt so that the ropt == 1 case
        # still yields a 1x1 matrix (diag(scalar) would mis-interpret
        # the scalar as a dimension).
        S_shared <- diag(rev(d), nrow = ropt)
        S_shared <- S_shared / norm(S_shared, "F")

        # Align table-specific loadings with updated S_shared for consistent
        # reconstruction
        V_list <- lapply(V_list, function(V){
            t(S_shared %*% t(V))
        })
    }

    # Reverse the earlier rescaling so singular values match the scale
    # of the original input data.
    S_shared <- S_shared / rescale_param

    # Ensure components are ordered by decreasing singular value.
    index_order <- order(diag(S_shared), decreasing = TRUE)
    U_shared <- U_shared[, index_order, drop = FALSE]
    S_shared <- S_shared[index_order, index_order, drop = FALSE]
    # Combine feature loadings into one matrix
    V_stacked <- do.call(rbind, V_list)
    V_stacked <- V_stacked[, index_order, drop = FALSE]

    # Return U, S, and V with CV errors
    res <- list(
        X = U_shared,
        S = S_shared,
        Y = V_stacked,
        cv_error = cv_errors
    )
    return(res)
}

# Initialize U, S and V by doing SVD with stacked data.
.initialize_joint_optspace <- function(
        observed_stacked, observed_list, mask_stacked, mask_list, ropt, eps,
        n_samples, n_features){
    # Run SVD. Our initial first guess are the loadings generated
    # by the traditional SVD.
    svd_res <- svd(observed_stacked)

    # Pad singular values to length `ropt` if svd() returned fewer
    # (rank-deficient case, BLAS-dependent).
    d <- svd_res[["d"]]
    if (length(d) < ropt) {
        d <- c(d, rep(0, ropt - length(d)))
    }
    d <- d[seq_len(ropt)]

    U_shared <- svd_res[["u"]][, seq_len(ropt), drop = FALSE]
    U_shared <- U_shared[
        , U_shared |> ncol() |> seq_len() |> rev(), drop = FALSE]
    S_shared <- diag(rev(d), nrow = ropt)
    V_shared <- svd_res[["v"]][, seq_len(ropt), drop = FALSE]
    V_shared <- V_shared[
        , V_shared |> ncol() |> seq_len() |> rev(), drop = FALSE]

    # The shape and number of non-zero values
    # can set the input parameters for the gradient
    # decent.
    U_shared <- U_shared * sqrt(n_samples)
    V_shared <- V_shared * sqrt(sum(n_features))
    S_shared <- S_shared / eps

    # Generate the new singular values from
    # the initialization of U and V
    S_shared <- .aux_getoptS(
        U_shared, V_shared, observed_stacked, mask_stacked)

    # Split feature loadings by table
    ends   <- cumsum(n_features)
    starts <- c(1, head(ends, -1) + 1)
    V_list <- mapply(function(s, e) {
        V_shared[seq(s, e), , drop = FALSE]
    }, starts, ends)

    res <- list(
        U_shared = U_shared,
        S_shared = S_shared,
        V_list = V_list
    )
    return(res)
}

# Perform one gradient descent update for a single table. For each table, this
# produces proposal for U_shared and S_shared which are then later averaged.
.gradient_update_joint_optspace <- function(
        table_i, observed_list, mask_list, U_shared, S_shared, V_list,
        rho, step.size = 1e4, sign.correction = -1){
    if( !(.is_an_integer(step.size) && step.size > 0) ){
        stop("'step.size' must be a single positive integer.", call. = FALSE)
    }
    if( !(.is_an_integer(sign.correction) && sign.correction %in% c(-1, 1)) ){
        stop("'sign.correction' must be either -1 or 1.", call. = FALSE)
    }
    obs  <- observed_list[[table_i]]
    mask <- mask_list[[table_i]]
    V_i  <- V_list[[table_i]]

    # Compute gradient for table i
    grad_res <- .aux_gradF_t(
        U_shared,
        V_i,
        S_shared,
        obs,
        mask,
        step.size,
        rho
    )
    U_update <- grad_res[["W"]]
    V_update <- grad_res[["Z"]]

    # Line search for optimal step
    step <- .aux_getoptT(
        U_shared,
        U_update,
        V_i,
        V_update,
        S_shared,
        obs,
        mask,
        step.size,
        rho
    )

    # Update table-specific matrices
    U_i <- U_shared - sign.correction * step * U_update
    V_i <- V_i - sign.correction * step * V_update

    # Recompute singular values for this table
    S_i <- .aux_getoptS(
        U_shared,
        V_i,
        obs,
        mask
    )

    # Return updated parameters
    res <- list(
        U_i = U_i,
        S_i = S_i,
        V_i = V_i

    )
    return(res)
}

# This function calculates error between original table and reconstructed table.
.calculate_optspace_cv_error <- function(x_test, S, V, ...){
    # Get mask, i.e., info on which values are NA
    mask <- x_test |> is.na()
    n_obs <- sum(!mask)
    # For the test tables, we double center them. The idea is that CV focuses
    # on latent structure reconstruction and not mean shifts.
    x_test <- x_test |> .apply_double_centering()
    # Treat missing entries as zero during multiplication
    x_filled <- x_test
    x_filled[mask] <- 0

    # Project test data onto feature loadings
    U_test <- x_filled %*% V
    U_test <- sweep(U_test, 2, diag(S), "/")

    # Reconstruct test matrix from projected U and singular values
    reconstruct_test <- U_test %*% S %*% t(V)

    # Keep NA structure like masked arrays
    reconstruct_test[mask] <- NA

    # Compute error on observed entries
    obs_error <- x_test - reconstruct_test
    obs_error[is.na(obs_error)] <- 0

    # Scale error by the number of observed entries
    error <- norm(obs_error, "F") / sqrt(n_obs)

    return(error)
}

# This function applies double centering of the data, i.e., it centers columns
# and rows.
.apply_double_centering <- function(mat, na.rm = TRUE, ...){
    mat <- sweep(mat, 1L, rowMeans(mat, na.rm = na.rm), "-")
    mat <- sweep(mat, 2L, colMeans(mat, na.rm = na.rm), "-")
    return(mat)
}

# Calculate population standard deviation
.population_sd <- function(x){
    (x - mean(x))^2 |> mean() |> sqrt()
}
