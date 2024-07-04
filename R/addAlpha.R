#' Estimate alpha diversity indices
#' 
#' The function estimates alpha diversity indices optionally using rarefaction,
#' then stores results in \code{\link{colData}}.
#' 
#' @param x a \code{\link{SummarizedExperiment}} object.
#' 
#' @param assay.type the name of the assay used for calculation of the
#'   sample-wise estimates (Default: \code{"counts"}).
#'   
#' @param index a \code{character} vector, specifying the alpha diversity 
#'   indices to be calculated.
#'   
#' @param name a name for the column(s) of the colData the results should be
#'   stored in. By default this will use the original names of the calculated
#'   indices(Default: \code{index}).
#' 
#' @param n.iter \code{NULL} or a single \code{integer} value for the number of
#'   rarefaction rounds. Rarefaction is not applied when \code{n.iter=NULL}
#'   (see @details section). (Default: \code{NULL}).
#'   
#' @param ... optional arguments passed to mia::rarefyAssay():
#' \itemize{
#'   \item a \code{numeric} value specifying the rarefaction depth i.e. the
#'   sample size drawn from samples.
#'   (Default: \code{min(colSums2(assay(x, assay.type)))})
#' }
#' 
#' @return \code{x} with additional \code{\link{colData}} named after the index 
#'    used.
#' 
#' @examples
#' 
#' data("GlobalPatterns")
#' tse <- GlobalPatterns
#' 
#' # Calculate the default Shannon index with no rarefaction
#' tse <- addAlpha(mae[[1]], index = c("shannon", "observed_richness"))
#' 
#' # Shows the estimated Shannon index
#' tse$shannon
#'
#' # Calculate observed richness with 10 rarefaction rounds
#' tse <- addAlpha(tse,
#'    assay.type = "counts",
#'    index = "observed_richness",
#'    sample=min(colSums(assay(tse, "counts")), na.rm = TRUE),
#'    n.iter=10)
#' 
#' # Shows the estimated observed richness
#' tse$observed_richness
#' 
#' @name addAlpha
#' @rdname addAlpha
#' @export
NULL

#' @rdname addAlpha
#' @export
setGeneric(
    "addAlpha", signature = c("x"),
    function(
        x, assay.type = "counts", 
        index = c(
            "coverage_diversity", "fisher_diversity", "faith_diversity",
            "gini_simpson_diversity", "inverse_simpson_diversity",
            "log_modulo_skewness_diversity", "shannon_diversity",
            "absolute_dominance", "dbp_dominance",
            "core_abundance_dominance", "gini_dominance",
            "dmn_dominance", "relative_dominance",
            "simpson_lambda_dominance", "camargo_evenness",
            "pielou_evenness", "simpson_evenness",
            "evar_evenness", "bulla_evenness", "ace_richness",
            "chao1_richness", "hill_richness", "observed_richness"),
        name = index, n.iter = NULL, ...)
    standardGeneric("addAlpha"))

#' @rdname addAlpha
#' @export
setMethod("addAlpha", signature = c(x = "SummarizedExperiment"),
    function(
        x, assay.type = "counts", 
        index = c(
            "coverage_diversity", "fisher_diversity", "faith_diversity",
            "gini_simpson_diversity", "inverse_simpson_diversity",
            "log_modulo_skewness_diversity", "shannon_diversity",
            "absolute_dominance", "dbp_dominance",
            "core_abundance_dominance", "gini_dominance",
            "dmn_dominance", "relative_dominance",
            "simpson_lambda_dominance", "camargo_evenness",
            "pielou_evenness", "simpson_evenness",
            "evar_evenness", "bulla_evenness", "ace_richness",
            "chao1_richness", "hill_richness", "observed_richness"),
        name = index, n.iter = NULL, ...){
        ############################## Input check #############################
        # Check that index is a character vector
        if( !.is_non_empty_character(index) ){
            stop("'index' should be a character vector.", call. = FALSE)
        }
        # if multiple indices to be estimated, name should a vector of
        # same length
        if( !.is_non_empty_character(name) || length(name) != length(index) ){
            stop(
                "'name' must be a non-empty character value and have the ",
                "same length than 'index'.",
                call. = FALSE)
        }
        # Check n.tier
        if( !(is.null(n.iter) || (.is_an_integer(n.iter) && n.iter >= 0)) ){
            stop("'n.iter' must be NULL or an integer.", call. = FALSE)
        }
        # Check if index exists. For each index input, detect it and get
        # information (e.g. internal function) to calculate the index.
        index <- .get_indices(index, name, x)
        ############################ Input check end ###########################
        # Looping over the vector of indices to be estimated
        for( i in seq_len(nrow(index)) ){
            # Performing rarefaction if sample is specified
            if( !is.null(n.iter) && n.iter > 0 ){
                x <- .alpha_rarefaction(
                    x, assay.type = assay.type, n.iter = n.iter,
                    FUN = index[i, "FUN"], index = index[i, "index"],
                    name = index[i, "name"], ...)
            } else {
                # Estimate index without rarefaction
                args <- c(
                    list(x, assay.type = assay.type, index = index[i, "index"],
                        name = index[i, "name"]), list(...))
                x <- do.call(index[i, "FUN"], args = args)
            }
        }
        return(x)
    }
)

################################ HELP FUNCTIONS ################################

# Search alpha diversity index that user wants to calculate.
.get_indices <- function(index, name, x){
    # Initialize list for supported indices
    supported <- list()
    # Supported diversity indices
    temp <- c(
        "coverage", "faith", "fisher", "gini_simpson", "inverse_simpson",
        "log_modulo_skewness", "shannon")
    temp <- data.frame(index = temp)
    temp[["measure"]] <- "diversity"
    temp[["index_long"]] <- paste0(temp[["index"]], "_", temp[["measure"]])
    temp[["FUN"]] <- ".estimate_diversity"
    supported[["diversity"]] <- temp
    # Supported dominance indices
    temp <- c(
        "absolute", "dbp", "core_abundance", "gini", "dmn", "relative",
        "simpson_lambda")
    temp <- data.frame(index = temp)
    temp[["measure"]] <- "dominance"
    temp[["index_long"]] <- paste0(temp[["index"]], "_", temp[["measure"]])
    temp[["FUN"]] <- ".estimate_dominance"
    supported[["dominance"]] <- temp
    # Supported eveness indices
    temp <- c(
        "camargo", "pielou", "simpson", "evar", "bulla")
    temp <- data.frame(index = temp)
    temp[["measure"]] <- "evenness"
    temp[["index_long"]] <- paste0(temp[["index"]], "_", temp[["measure"]])
    temp[["FUN"]] <- ".estimate_evenness"
    supported[["eveness"]] <- temp
    # Supported richness indices
    temp <- c(
        "ace", "chao1", "hill", "observed")
    temp <- data.frame(index = temp)
    temp[["measure"]] <- "richness"
    temp[["index_long"]] <- paste0(temp[["index"]], "_", temp[["measure"]])
    temp[["FUN"]] <- ".estimate_richness"
    supported[["richness"]] <- temp
    # Combine
    supported <- do.call(rbind, supported)
    # Find the index that user wants to calculate
    ind <- match(tolower(index), supported[["index_long"]])
    ind_short <- match(tolower(index), supported[["index"]])
    ind[ is.na(ind) ] <- ind_short[ is.na(ind) ]
    detected <- supported[ind, ]
    # Add the index that was searched
    detected[["search"]] <- index
    # Add names that user wants to use when storing results to colData
    detected[["name"]] <- name
    # Check if there are indices that were not detected
    if( any(is.na(detected[["index"]])) ){
        not_detected <- paste0(
            detected[is.na(detected[["index"]]), "search"], collapse = "', '")
        not_detected <- paste0("'", not_detected, "'")
        stop(
            "'index' is corresponding to none of the alpha diversity ",
            "indices. The following 'index' was not detected: ", not_detected,
            call. = FALSE)
    }
    # Faith index is available only for TreeSE with rowTree
    if( "faith" %in% detected[["index"]] &&
            !(is(x, "TreeSummarizedExperiment") && !is.null(rowTree(x))) ){
        # Drop faith index from indices being calculated
        detected <- detected[!detected[["index"]] %in% c("faith"), ]
        # If there are still other indices being calculated, give warning.
        # Otherwise, give error if faith was the only index that user wants to
        # calculate.
        FUN <- if( nrow(detected) == 0 ) stop else warning
        FUN("'faith' index can be calculated only for TreeSE with rowTree(x) ",
            "populated.", call. = FALSE)
    }
    # Check if there are indices left
    return(detected)
}

# This function rarifies the data n.iter of times and calculates index for the
# rarified data. The result is a mean of the iterations.
#' @importFrom DelayedMatrixStats colMeans2
.alpha_rarefaction <- function(
        x, assay.type, n.iter, FUN, index, name, ...){
    # Calculating the mean of the subsampled alpha estimates ans storing them
    res <- lapply(seq(n.iter), function(i){
        # Subsampling the counts from the original TreeSE object.
        x_sub <- rarefyAssay(x, assay.type = assay.type, verbose = FALSE, ...)
        # Calculating the diversity indices on the subsampled object
        x_sub <- do.call(FUN, args = list(
            x_sub, assay.type = assay.type, index = index,
            name = "rarefaction_temp_result", list(...)))
        # Get results as a vector from colData
        temp <- colData(x_sub)[["rarefaction_temp_result"]]
        return(temp)
    })
    # Combine list of vectors from multiple iterations
    res <- do.call(rbind, res)
    # Calculate mean of iterations
    res <- colMeans2(res)
    # Give warning about missing samples. Same might have been dropped during
    # rarefaction.
    if( !all(colnames(x) %in% names(res)) ){
        warning(
            "Some samples were dropped during rarefaction leading to missing ",
            "diversity values. Consider lower 'sample'.", call. = FALSE)
    }
    # It might be that certain samples were dropped off if they have lower
    # abundance than rarefaction  depth --> order so that data includes all the
    # samples
    res <- res[match(colnames(x), names(res))]
    res <- unname(res)
    # Add to original data
    colData(x)[[name]] <- res
    return(x)
}
