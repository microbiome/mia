#' @name
#' getRPCA
#'
#' @title
#' Add here
#'
#' @description
#' Add here
#'
#' @details
#' Add here.
#'
#' @return
#' \code{SummarizedExperiment} object.
#'
#' @inheritParams addAlpha
#'
#' @param assay.type \code{Character scalar}. Specifies the name of assay
#' used in calculation. (Default: \code{"counts"})
#'
#' @param ... additional arguments.
#'
#' @examples
#'
#' data(GlobalPatterns)
#' tse <- GlobalPatterns
#'
#' @seealso
#' \code{\link[scater::runPCA]{runPCA}}
#'
#' @references
#'
#' Add here
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
        res <- .calculate_joint_rpca(mat_list, ...)
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
        res <- .calculate_joint_rpca(mat_list, ...)
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
        x <- .add_values_to_mae_metadata(x, res, name = name, ...)
        return(x)
    }
)
##########################

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

# This function runs RPCA to single table
.calculate_rpca <- function(mat, ncomponents = 3L, ...){
    # Get lower rank representation of the data
    opt_results <- .get_lower_rank_mat(mat, ncomponents = ncomponents, ...)
    # The result might have lower number of columns if they were not able to be
    # estimated.
    ncomponents <- opt_results[["raw"]][["S"]] |> ncol()
    # Apply pca to lower rank representation
    pca_results <- .calculate_pca(
        opt_results[["matrix"]], ncomponents = ncomponents, ...)
    # Calculate distance in PCA space
    distance <- pca_results[["sample_scores"]] |> dist()
    # Create a final results to return to user
    res <- .construct_rpca_result(pca_results, opt_results, distance)
    return(res)
}

# This function calculates Joint-RPCA for multiple tables
.calculate_joint_rpca <- function(mat_list, test.set = NULL, ...){
    if( !(is.null(test.set) || is.character(test.set)) ){
        stop(".")
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

    # Create stacked tables for train and test sets
    train_mat <- do.call(cbind, train_set)
    test_mat <- do.call(cbind, test_set)

    # Calculate RPCA for stacked train set. The results are calculated as
    # features were from the same table.
    pca_result <- .calculate_rpca(train_mat, ...)

    # Project test samples to PCA space determined by train set
    projected <- .project_test_set_to_rpca(pca_result, test_mat)
    # Add test samples to pca results
    attr_list <- attributes(pca_result)
    pca_result <- rbind(pca_result, projected)
    # Sort back to original order
    pca_result <- pca_result[
        order(c(train_samples, test_samples)), , drop = FALSE]
    # Add additional info back
    attr_list <- c(attributes(pca_result), attr_list)
    attr_list <- attr_list[ !duplicated(names(attr_list)) ]
    attributes(pca_result) <- attr_list

    # Calculate error between low rank and original test set
    num_features <- vapply(mat_list, ncol, numeric(1L))
    reconstruct_error <- .calculate_reconstruct_error(
        pca_result, test_set, num_features)
    attributes(pca_result)[["reconstruct_error"]] <- reconstruct_error

    return(pca_result)
}

# This function constructs a lower rank representation from the data. The idea
# is to extract the essential from the data and to remove noise.
#' @importFrom vegan optspace
.get_lower_rank_mat <- function(
        mat,
        ncomponents = pmin(3L, nrow(mat), ncol(mat)),
        max.iterations = 5L,
        tol = 1e-5,
        ...){
    # Add input check

    # Create lower rank representation
    opt_result <- optspace(
        x = mat,
        ropt = ncomponents,
        niter = max.iterations,
        tol = tol,
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

# This function applies PCA to the data.
.calculate_pca <- function(
        mat, ncomponents, center = TRUE, scale = FALSE, ...){

    if( center ){
        # Row and column means
        row_means <- rowMeans(mat)
        col_means <- colMeans(mat)
        grand_mean <- mean(mat)

        # Do double centering. Add overall mean so that we do not substract the
        # data effectively 2 times. The result is a matrix that has row and
        # column means in zero.
        mat <- sweep(mat, 1L, row_means, "-")
        mat <- sweep(mat, 2L, col_means, "-")
        mat <- mat + grand_mean
    } else {
        row_means <- NULL
        col_means <- NULL
        grand_mean <- NULL
    }

    # Optional column scaling (PCA on correlation matrix)
    if( scale ){
        col_sds <- colSds(mat)
        col_sds[col_sds == 0] <- 1
        mat <- sweep(mat, 2L, col_sds, "/")
    } else {
        col_sds <- NULL
    }

    # Run PCA
    svd_result <- mat |> svd()
    u <- svd_result[["u"]]
    s <- svd_result[["d"]]
    v <- svd_result[["v"]]

    # We multiply U by singular values so that the sample coordinates reflect
    # actual variance magnitude rather than just orthonormal directions.
    # u <- u %*% diag(s) ###################################################### UNCOMMENT

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
        center = c(
            row = row_means,
            col = col_means,
            grand = grand_mean
        ),
        scale = col_sds
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
        mat, test.ratio = 0.2, ...){
    # Select number of samples
    test_n <- ceiling(test.ratio * nrow(mat))
    test_samples <- c()
    if( test_n > 0L ){
        # Calculate RPCA
        pca_result <- .calculate_rpca(mat, ...)
        # Select N samples so that they span over the PC1 axis. Idea is that
        # these samples should represent the dataset the best as there are
        # maximally different samples.
        first_component <- pca_result[, 1] |> sort()
        test_samples <- seq(1, length(first_component), length.out = test_n) |>
            round()
        test_samples <- names(first_component)[ test_samples ]
        test_samples <- match(test_samples, rownames(pca_result))
    }
    return(test_samples)
}

# This function projects test set samples to PCA space that were obtained with
# train set.
.project_test_set_to_rpca <- function(pca_result, test_mat){
    # Get results
    sample_scores <- pca_result
    feature_scores <- attributes(pca_result)[["rotation"]]
    eigenvalues <- attributes(pca_result)[["varExplained"]]

    # Calculate projection
    projected <- test_mat %*% feature_scores %*% diag(1/eigenvalues)

    return(projected)
}

# This function calculates error between raw test set values and their
# lower rank representation. The lower rank representation is calculated based
# on parameters learned from train set. The idea is to assess, how well the
# lower rank representation learns the generic, generalizable patterns from the
# data.
.calculate_reconstruct_error <- function(pca_result, test_set, num_features){
    # Get learned parametes
    u_shared <- attributes(pca_result)[["X"]]
    s_shared <- attributes(pca_result)[["S"]]
    y_shared <- attributes(pca_result)[["Y"]]

    # Split feature loadings by table
    ends   <- cumsum(num_features)
    starts <- c(1, head(ends, -1) + 1)
    y_individual <- mapply(function(s, e) {
        y_shared[s:e, ]
    }, starts, ends)

    # Calculate error separately for each table
    errors_per_set <- vapply(seq_len(length(test_set)), function(i){
        # Calculate lower rank representation
        test_mat <- test_set[[i]]
        u_test <- test_mat %*% y_individual[[i]]
        u_test <- sweep(u_test, 2, diag(s_shared), "/")
        recon_test <- u_test %*% s_shared %*% t(y_individual[[i]])
        ########################################################################### REMOVE THESE LINES
        # Center for consistency
        recon_test <- scale(recon_test, center = TRUE, scale = FALSE)
        recon_test <- t(scale(t(recon_test), center = TRUE, scale = FALSE))
        # Calculate error between actual values and lower rank representation
        error <- test_mat - recon_test
        error[is.na(error)] <- 0
        error <- norm(error, "F") / sqrt(sum(!is.na(test_mat)))
        return(error)
    }, numeric(1L))
    names(errors_per_set) <- names(num_features)

    return(errors_per_set)
}
