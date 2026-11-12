#' Merge SE objects into single SE object.
#'
#' @inheritParams rarefyAssay
#' @inheritParams getDominant
#'
#' @param y a \code{\link[SummarizedExperiment]{SummarizedExperiment}} object
#' when \code{x} is a \code{\link[SummarizedExperiment]{SummarizedExperiment}}
#' object. Must be missing when \code{x} is a list.
#'
#' @param join \code{Character scalar}. A value for selecting the joining
#' method. Must be 'full', 'inner', 'left', or 'right'. 'left' and 'right' are
#' disabled when more than two objects are being merged.
#' (Default: \code{"full"})
#'
#' @param missing.values \code{NA}, \code{0} or \code{Character scalar}.
#' Specifies the notation of missing values. (By default: \code{NA})
#'
#' @param missing_values Deprecated. Use \code{missing.values} instead.
#'
#' @param collapse.cols \code{Logical scalar}. Determines whether to collapse
#' identically named samples to one. (Default: \code{FALSE})
#'
#' @param collapse_samples Deprecated. Use \code{collapse.cols} instead.
#'
#' @param collapse.rows \code{Logical scalar}. Determines whether to collapse
#' features that have identical names and taxonomy to one.
#' (Default: \code{TRUE})
#'
#' @param collapse_features Deprecated. Use \code{collapse.rows} instead.
#'
#' @param ... optional arguments (not used).
#'
#' @return A single \code{SummarizedExperiment} object.
#'
#' @details
#' This function merges multiple \code{SummarizedExperiment} objects. It
#' combines \code{assays}, \code{rowData}, \code{colData}, and
#' \code{metadata} so that the output includes each unique row and column.
#' The merging is done based on \code{rownames} and \code{colnames}.
#'
#' Rows are matched based on \code{rownames} and taxonomy ranks of
#' \code{rowData}, i.e., rows are interpreted as equal if they have the same
#' name and taxonomy (unless \code{collapse.rows = FALSE}). Rows that have
#' the same name but differently annotated taxonomy are kept as separate rows.
#' To merge data at certain taxonomy level, agglomerate the objects first with
#' \code{\link{agglomerateByRank}}. Duplicated
#' names within an object are made unique. For samples, collapsing is disabled by default meaning that
#' equally named samples that are stored in different objects are interpreted
#' as unique. Collapsing can be enabled with \code{collapse.cols = TRUE} when
#' equally named samples describe the same sample. Equally named rows or
#' columns that are not collapsed are made unique with
#' \code{\link[base:make.unique]{make.unique}}. When the same value is found
#' from multiple objects, the value from the first object is used. If equally
#' named \code{rowData} or \code{colData} variables have different classes,
#' they are combined as character. Taxonomy rank columns of \code{rowData}
#' are matched case-insensitively, e.g., \code{"genus"} and \code{"Genus"}.
#'
#' If, for example, all rows are not shared with individual objects, there are
#' missing values in \code{assays}. The notation of missing can be specified
#' with the \code{missing.values} argument.
#'
#' If input consists of \code{TreeSummarizedExperiment} objects, also
#' \code{referenceSeq}, \code{rowTree}, and \code{colTree} are preserved if
#' possible. Reference sequences are merged if a sequence is available for
#' every row. Trees are not merged; a tree is added only if a single input
#' tree includes all the rows (or columns) of the merged data. Otherwise, the
#' tree is dropped with a warning.
#'
#' Compared to \code{cbind} and \code{rbind} \code{mergeSEs}
#' allows more freely merging since \code{cbind} and \code{rbind} expect
#' that rows and columns are matching, respectively.
#'
#' You can choose joining methods from \code{'full'}, \code{'inner'},
#'  \code{'left'}, and  \code{'right'}. In all the methods, all the samples are
#'  included in the result object. However, with different methods, it is
#'  possible
#'  to choose which rows are included.
#'
#' \itemize{
#'   \item{\code{full} -- all unique features}
#'   \item{\code{inner} -- all shared features}
#'   \item{\code{left} -- all the features of the first object}
#'   \item{\code{right} -- all the features of the second object}
#' }
#'
#' The output depends on the input. If the input contains
#' \code{SummarizedExperiment}
#' object, then the output will be \code{SummarizedExperiment}. When all the
#' input
#' objects belong to \code{TreeSummarizedExperiment}, the output will be
#' \code{TreeSummarizedExperiment}.
#'
#' @seealso
#' \itemize{
#'   \item{\code{TreeSummarizedExperiment::cbind}}
#'   \item{\code{TreeSummarizedExperiment::rbind}}
#'   \item{\code{\link[dplyr:full_join]{full_join}}}
#'   \item{\code{\link[dplyr:inner_join]{inner_join}}}
#'   \item{\code{\link[dplyr:left_join]{left_join}}}
#'   \item{\code{\link[dplyr:right_join]{right_join}}}
#' }
#'
#' @name mergeSEs
#' @export
#'
#' @examples
#' data(GlobalPatterns)
#' data(esophagus)
#' data(enterotype)
#'
#' # Take only subsets so that it wont take so long
#' tse1 <- GlobalPatterns[1:100, ]
#' tse2 <- esophagus
#' tse3 <- enterotype[1:100, ]
#'
#' # Merge two TreeSEs
#' tse <- mergeSEs(tse1, tse2)
#'
#' # Merge a list of TreeSEs
#' list <- SimpleList(tse1, tse2, tse3)
#' tse <- mergeSEs(list, assay.type = "counts", missing.values = 0)
#' tse
#'
#' # With 'join', it is possible to specify the merging method. Subsets are used
#' # here just to show the functionality
#' tse_temp <- mergeSEs(tse[1:10, 1:10], tse[5:100, 11:20], join = "left")
#' tse_temp
#'
#' # If your objects contain samples that describe one and same sample,
#' # you can collapse equally named samples to one by specifying 'collapse.cols'
#' tse_temp <- mergeSEs(list(tse[1:10, 1], tse[1:20, 1], tse[1:5, 1]),
#'                        collapse.cols = TRUE,
#'                        join = "inner")
#' tse_temp
#'
#' # Merge all available assays
#' tse <- transformAssay(tse, method="relabundance")
#' ts1 <- transformAssay(tse1, method="relabundance")
#' tse_temp <- mergeSEs(tse, tse1, assay.type = assayNames(tse))
#'
NULL

###################### Function for SimpleList of TreeSEs ######################

#' @rdname mergeSEs
#' @export
setMethod("mergeSEs", signature = c(x = "SimpleList", y = "missing"),
    function(x, y, assay.type = assay_name, assay_name = "counts",
            join = "full", missing.values = missing_values, missing_values = NA,
            collapse.cols = collapse_samples, collapse_samples = FALSE,
            collapse.rows = collapse_features, collapse_features = TRUE, ...){
        ################## Input check ##################
        x <- x |> as.list() |> unname()
        if( !all(vapply(x, is, logical(1L), "SummarizedExperiment")) ){
            stop("Input includes an object that is not ",
                "'SummarizedExperiment'.", call. = FALSE)
        }
        if( !(.is_a_string(join) &&
                join %in% c("full", "inner", "left", "right")) ){
            stop("'join' must be 'full', 'inner', 'left', or 'right'.",
                call. = FALSE)
        }
        if( length(x) > 2L && join %in% c("left", "right") ){
            stop("Joining method 'left' and 'right' are not available ",
                "when more than two objects are being merged.", call. = FALSE)
        }
        if( !(length(missing.values) == 1L && (is.na(missing.values) ||
                .is_a_string(missing.values) || .is_zero(missing.values))) ){
            stop("'missing.values' must be 0, NA, or a single character ",
                "value.", call. = FALSE)
        }
        if( !.is_a_bool(collapse.cols) ){
            stop("'collapse.cols' must be TRUE or FALSE.", call. = FALSE)
        }
        if( !.is_a_bool(collapse.rows) ){
            stop("'collapse.rows' must be TRUE or FALSE.", call. = FALSE)
        }
        # Keep assays that are found from all the objects
        found <- vapply(assay.type, function(a)
            all(vapply(x, function(obj) a %in% assayNames(obj), logical(1L))),
            logical(1L))
        if( !any(found) ){
            stop("'assay.type' must specify an assay that is found from all ",
                "the objects.", call. = FALSE)
        }
        if( !all(found) ){
            warning("The following assay(s) was not found from all the ",
                "objects so it is dropped from the output: ",
                paste0("'", assay.type[!found], "'", collapse = ", "),
                call. = FALSE)
        }
        ################ Input check end ################
        .merge_SEs(x, assay.type[found], join, missing.values, collapse.cols,
            collapse.rows, ...)
    }
)

########################### Function for two TreeSEs ###########################

#' @rdname mergeSEs
#' @export
setMethod("mergeSEs",
    signature = c(x = "SummarizedExperiment", y = "SummarizedExperiment"),
    function(x, y, ...){
        mergeSEs(SimpleList(x, y), ...)
    }
)

########################### Function for list TreeSEs ##########################

#' @rdname mergeSEs
#' @export
setMethod("mergeSEs", signature = c(x = "list", y = "missing"),
    function(x, y, ...){
        mergeSEs(SimpleList(x), ...)
    }
)

################################ HELP FUNCTIONS ################################

# This function merges a list of SE objects into one SE. The final rows and
# columns are determined first, and then each component is filled in based on
# the position of each row and column in the output. When the same value is
# found from multiple objects, the value from the first object is used.
#' @importFrom SingleCellExperiment SingleCellExperiment
.merge_SEs <- function(
        x, assay.type, join, missing.values, collapse.cols, collapse.rows, ...){
    rows_list <- lapply(x, rownames)
    cols_list <- lapply(x, colnames)
    if( any(lengths(c(rows_list, cols_list)) == 0L) ){
        stop("Each object must have rownames and colnames.", call. = FALSE)
    }
    # Make names unique within each object
    rows_list <- lapply(rows_list, make.unique)
    cols_list <- lapply(cols_list, make.unique)
    # rowData and colData are handled as lists for speed
    rds <- lapply(x, function(y) as.list(rowData(y))) |>
        .standardize_rank_names()
    cds <- lapply(x, function(y) as.list(colData(y)))
    # Rows are matched based on rownames and taxonomy ranks
    keys_list <- .get_row_keys(rds, rows_list)
    # Make names unique between objects if they should not be collapsed
    if( !collapse.rows ){
        keys_list <- .make_names_unique(keys_list)
    }
    if( !collapse.cols ){
        cols_list <- .make_names_unique(cols_list)
    }
    keys <- switch(join,
        full = keys_list |> unlist() |> unique(),
        inner = Reduce(intersect, keys_list),
        left = keys_list[[1L]],
        right = keys_list[[2L]])
    cols <- cols_list |> unlist() |> unique()
    # Position of each row and column in the output
    ri_list <- keys_list |> unlist() |> match(keys) |> .split_as(keys_list)
    ci_list <- cols_list |> unlist() |> match(cols) |> .split_as(cols_list)
    # Rows are named by the original rownames. Rows that have equal rownames
    # but different taxonomy are made unique.
    rows <- unlist(rows_list)[ match(keys, unlist(keys_list)) ] |>
        make.unique()
    # Merge data
    assays <- lapply(assay.type, function(a)
        lapply(x, assay, a, withDimnames = FALSE) |>
            .merge_assays(ri_list, ci_list, rows, cols, missing.values, ...))
    names(assays) <- assay.type
    args <- list(
        assays = assays,
        rowData = .merge_tables(rds, ri_list, rows),
        colData = .merge_tables(cds, ci_list, cols),
        metadata = lapply(x, metadata) |> do.call(what = c))
    # The output class is the most specific class shared by all the objects
    if( all(vapply(x, is, logical(1L), "TreeSummarizedExperiment")) ){
        args[["referenceSeq"]] <- .merge_reference_seqs(x, ri_list, rows)
        tse <- do.call(TreeSummarizedExperiment, args)
        # Trees are added with changeTree since the constructor would copy the
        # assays.
        tse <- tse |>
            .add_tree(x, ri_list, "row") |>
            .add_tree(x, ci_list, "col")
    } else if( all(vapply(x, is, logical(1L), "SingleCellExperiment")) ){
        tse <- do.call(SingleCellExperiment, args)
    } else{
        tse <- do.call(SummarizedExperiment, args)
    }
    if( length(unique(lapply(x, class))) > 1L ){
        warning("The input consists of multiple classes. The output is '",
            class(tse), "'.", call. = FALSE)
    }
    return(tse)
}

# This function combines rownames with taxonomy ranks of rowData, similarly
# to agglomerateByRank. Missing ranks are treated as empty.
.get_row_keys <- function(rds, rows_list){
    ranks <- lapply(rds, names) |> unlist() |>
        intersect(x = .capitalize(TAXONOMY_RANKS))
    if( length(ranks) == 0L ){
        return(rows_list)
    }
    mapply(function(rd, nams){
        tax <- lapply(ranks, function(rank){
            v <- as.character(rd[[rank]])
            v[ is.na(v) ] <- ""
            v
        })
        do.call(paste, c(list(nams), tax, sep = "_"))
    }, rds, rows_list, SIMPLIFY = FALSE)
}

# This function makes names unique across a list of name vectors.
.make_names_unique <- function(nams){
    nams |> unlist() |> make.unique() |> .split_as(nams)
}

# This function splits a vector into a list that has the same lengths as the
# given list.
.split_as <- function(v, lst){
    split(v, rep(seq_along(lst), lengths(lst))) |> unname()
}

# This function merges assays into a single matrix. Objects are looped in
# reverse order so that values of the first object are kept.
#' @importFrom DelayedArray type
.merge_assays <- function(
        mats, ri_list, ci_list, rows, cols, missing.values, dense = 1e8, ...){
    # Sparse matrix is used if the dense matrix would be too large, i.e., it
    # would have more cells than 'dense'. It is possible only if missing values
    # are zeros.
    too_large <- as.numeric(length(rows)) * length(cols) > dense
    if( too_large && .is_zero(missing.values) ){
        return(.merge_assays_sparse(mats, ri_list, ci_list, rows, cols))
    }
    # Initialize with the type of data so that matrix is not copied when filled
    init <- c(vector(type(mats[[1L]])), missing.values)
    res <- matrix(init, length(rows), length(cols), dimnames = list(rows, cols))
    for( i in rev(seq_along(mats)) ){
        keep <- !is.na(ri_list[[i]])
        res[ri_list[[i]][keep], ci_list[[i]]] <-
            as.matrix(mats[[i]][keep, , drop = FALSE])
    }
    return(res)
}

# This function merges assays into a sparse matrix. Non-zero values are
# collected from each object, and the first value of each cell is kept.
.merge_assays_sparse <- function(mats, ri_list, ci_list, rows, cols){
    .require_package("Matrix")
    vals <- lapply(seq_along(mats), function(k){
        m <- as(mats[[k]], "CsparseMatrix") |> Matrix::summary()
        list(i = ri_list[[k]][ m[["i"]] ], j = ci_list[[k]][ m[["j"]] ],
            x = m[["x"]])
    })
    i <- lapply(vals, `[[`, "i") |> unlist()
    j <- lapply(vals, `[[`, "j") |> unlist()
    x <- lapply(vals, `[[`, "x") |> unlist()
    keep <- !is.na(i) & !duplicated(i + (j - 1) * length(rows))
    Matrix::sparseMatrix(
        i = i[keep], j = j[keep], x = x[keep],
        dims = c(length(rows), length(cols)), dimnames = list(rows, cols))
}

# Checks whether value is numeric zero
.is_zero <- function(x){
    isTRUE(is.numeric(x) && x == 0)
}

# This function merges rowData or colData tables. Tables are looped in reverse
# order so that non-missing values of the first table are kept.
#' @importFrom S4Vectors DataFrame
.merge_tables <- function(dfs, ind_list, ids){
    col_names <- lapply(dfs, names) |> unlist() |> unique()
    res <- lapply(col_names, function(col){
        has_col <- vapply(dfs, function(df) col %in% names(df), logical(1L))
        vals <- lapply(dfs[has_col], `[[`, col) |> .harmonize_classes()
        inds <- ind_list[has_col]
        # Missing values with the same class as the data
        out <- vals[[1L]][ rep(NA_integer_, length(ids)) ] |> unname()
        for( i in rev(seq_along(vals)) ){
            keep <- !is.na(inds[[i]]) & !is.na(vals[[i]])
            out[ inds[[i]][keep] ] <- vals[[i]][keep]
        }
        out
    })
    names(res) <- col_names
    res <- DataFrame(res, check.names = FALSE)
    if( length(col_names) == 0L ){
        res <- DataFrame(row.names = ids)
    }
    rownames(res) <- ids
    return(res)
}

# This function renames taxonomy rank columns with the spelling of
# TAXONOMY_RANKS so that they are matched case-insensitively.
.standardize_rank_names <- function(dfs){
    lapply(dfs, function(df){
        ind <- names(df) |> tolower() |> match(TAXONOMY_RANKS)
        names(df)[ !is.na(ind) ] <- TAXONOMY_RANKS[ ind[!is.na(ind)] ] |>
            .capitalize()
        df
    })
}

# This function converts columns into a class that can hold all of them.
.harmonize_classes <- function(vals){
    classes <- vapply(vals, function(v) class(v)[[1L]], character(1L))
    if( all(classes == "factor") ){
        lvls <- lapply(vals, levels) |> unlist() |> unique()
        vals <- lapply(vals, factor, levels = lvls)
    } else if( any(classes != classes[[1L]]) ){
        numeric <- all(vapply(vals, is.numeric, logical(1L)))
        vals <- lapply(vals, if( numeric ) as.numeric else as.character)
    }
    return(vals)
}

# Trees are not merged. A tree is added only if a single input tree includes
# all the rows (or columns) of the merged data.
#' @importFrom TreeSummarizedExperiment rowLinks colLinks rowTree colTree
#'     changeTree
.add_tree <- function(tse, x, ind_list, MARGIN){
    links_FUN <- switch(MARGIN, row = rowLinks, col = colLinks)
    tree_FUN <- switch(MARGIN, row = rowTree, col = colTree)
    slot_name <- paste0(MARGIN, "Tree")
    n <- switch(MARGIN, row = nrow(tse), col = ncol(tse))
    # Get all the trees as a list
    trees <- lapply(x, tree_FUN, NULL) |> unlist(recursive = FALSE) |>
        unique()
    if( length(trees) == 0L || n == 0L ){
        return(tse)
    }
    # Get node label of each row from the first object that has the link
    labels <- rep(NA_character_, n)
    for( i in rev(seq_along(x)) ){
        labs <- links_FUN(x[[i]])[["nodeLab"]]
        keep <- !is.na(ind_list[[i]]) & !is.na(labs)
        labels[ ind_list[[i]][keep] ] <- labs[keep]
    }
    for( tree in trees ){
        if( all(labels %in% c(tree[["tip.label"]], tree[["node.label"]])) ){
            args <- list(tse, tree, labels)
            names(args) <- c("x", slot_name, paste0(MARGIN, "NodeLab"))
            # Multiple rows can be linked to the same node. TreeSE warns about
            # it, but it is expected here.
            tse <- do.call(changeTree, args)
            return(tse)
        }
    }
    warning(slot_name, " is dropped since none of the input trees includes ",
        "all the ", MARGIN, "s. Trees are not merged.", call. = FALSE)
    return(tse)
}

# This function merges reference sequences. They are added only if every row
# has a sequence.
.merge_reference_seqs <- function(x, ri_list, rows){
    seqs <- lapply(x, referenceSeq)
    has_seqs <- vapply(seqs, is, logical(1L), "DNAStringSet")
    if( !any(has_seqs) || length(rows) == 0L ){
        return(NULL)
    }
    # Position of first sequence of each row
    ind <- ri_list[has_seqs] |> unlist() |> match(x = seq_along(rows))
    if( anyNA(ind) ){
        warning("referenceSeq is dropped since sequences are not available ",
            "for all the rows.", call. = FALSE)
        return(NULL)
    }
    seqs <- seqs[has_seqs] |> unname() |> do.call(what = c)
    seqs <- seqs[ind]
    names(seqs) <- rows
    return(seqs)
}
