#' Estimate alpha diversity indices.
#' 
#' The function estimates alpha diversity indices optionally using of rarefaction,
#' then stores results at \code{\link{colData}}.
#' 
#' @param x a \code{\link{SummarizedExperiment}} object.
#' 
#' @param assay.type the name of the assay used for
#'   calculation of the sample-wise estimates (default: \code{assay.type = "counts"}).
#'   
#' @param index a \code{character} vector, specifying the alpha diversity indices
#'   to be calculated.
#'   
#' @param name a name for the column(s) of the colData the results should be
#'   stored in. By default this will use the original names of the calculated
#'   indices(By default: \code{name = index}).
#'   
#' @param ... optional arguments.
#' 
#' @param n.iter a single \code{integer} value for the number of rarefaction
#' rounds (By default: \code{n.iter = 10}).
#' 
#' @param rarefaction.depth a \code{double} value as for the minimim size or 
#' rarefaction.depth. (By default: \code{rarefaction.depth = NULL})
#' 
#' @return \code{x} with additional \code{\link{colData}} named after the index 
#' used.
#' 
#' @examples
#' 
#' data("GlobalPatterns")
#' tse <- GlobalPatterns
#' 
#' # Calculate the default Shannon index with no rarefaction
#' tse <- addAlpha(tse, assay.type = "counts", index = "shannon")
#' 
#' # Shows the estimated Shannon index
#' tse$shannon
#'
#' # Calculate observed richness with 10 rarefaction rounds
#' tse <- addAlpha(tse, assay.type = "counts", index = "observed_richness",
#' rarefaction.depth=min(colSums(assay(tse, "counts")), na.rm = TRUE), n.iter=10)
#' 
#' # Shows the estimated observed richness
#' tse$observed_richness
#' 

#' @name addAlpha
#' @rdname addAlpha
#' @export
setGeneric(
    "addAlpha", signature = c("x"), function(
        x, assay.type = "counts", 
        index = c(
            "coverage_diversity", "fisher_diversity", "faith_diversity",
            "gini_simpson_diversity", "inverse_simpson_diversity",
            "log_modulo_skewness_diversity", "shannon_diversity",
            "absolute_dominance", "dbp_dominance", "core_abundance_dominance",
            "gini_dominance", "dmn_dominance", "relative_dominance",
            "simpson_lambda_dominance", "camargo_evenness", "pielou_evenness",
            "simpson_evenness", "evar_evenness", "bulla_evenness",
            "ace_richness", "chao1_richness", "hill_richness",
            "observed_richness"),
        name = index, n.iter = 10, rarefaction.depth = NULL, ...)
    standardGeneric("addAlpha"))

#' @rdname addAlpha
#' @export
setMethod(
    "addAlpha", signature = c(x = "SummarizedExperiment"), function(
        x, assay.type = "counts", 
        index = c(
            "coverage_diversity", "fisher_diversity", "faith_diversity",
            "gini_simpson_diversity", "inverse_simpson_diversity",
            "log_modulo_skewness_diversity", "shannon_diversity",
            "absolute_dominance", "dbp_dominance", "core_abundance_dominance",
            "gini_dominance", "dmn_dominance", "relative_dominance",
            "simpson_lambda_dominance", "camargo_evenness", "pielou_evenness",
            "simpson_evenness", "evar_evenness", "bulla_evenness",
            "ace_richness", "chao1_richness", "hill_richness",
            "observed_richness"),
        name = index, n.iter = 10, rarefaction.depth = NULL, ...){
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
        if( !.is_an_integer(n.iter) ) {
            stop("'n.iter' must be an integer.", call. = FALSE)
        }
        # Check that rarefaction.depth is a numeric > 0
        if( !is.null(rarefaction.depth) && 
            !(is.numeric(rarefaction.depth) && rarefaction.depth > 0)) {
            stop("'rarefaction.depth' must be a non-zero positive double.",
                 call. = FALSE)
        }
        # Check if index exists
        index <- lapply(index, .get_indices)
        index <- do.call(rbind, index)
        index[["name"]] <- name
        if( any(is.na(index[["index"]])) ){
            stop(
                "'index' is corresponding to none of the alpha diversity ",
                "indices. The following 'index' was not detected: ",
                paste0(
                    index[is.na(index[["index"]]), "search"], collapse = ", "),
                call. = FALSE)
        }
        ############################ Input check end ###########################
        # Looping over the vector of indices to be estimated
        for( i in seq_len(nrow(index)) ){
            # Performing rarefaction if rarefaction.depth is specified
            if( !is.null(rarefaction.depth) ){
                x <- .alpha_rarefaction(
                    x, assay.type = assay.type, n.iter = n.iter,
                    rarefaction.depth = rarefaction.depth,
                    FUN = index[i, "FUN"], index = index[i, "index"],
                    name = index[i, "name"], ...)
            } else {
                # Estimate index without rarefaction
                x <- do.call(
                    index[i, "FUN"], args = c(
                        list(x, assay.type = assay.type,
                             index = index[i, "index"],
                             name = index[i, "name"]),
                        list(...)))
            }
        }
        return(x)
    }
)

################################ HELP FUNCTIONS ################################

# Search index that user wants to calculate.
.get_indices <- function(index) {
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
    # Find the index that user wanst to calculate
    ind <- index == supported[["index"]] | index == supported[["index_long"]]
    detected <- supported[ind, ]
    # If not found, create an empty vector
    if( nrow(detected) == 0 ){
        detected <- rep(NA, ncol(supported))
        names(detected) <- c("index", "measure", "index_long", "FUN")
    }
    # Add the index that was searched
    detected[["search"]] <- index
    return(detected)
}

# This function rarifies the data n.iter of times and calculates index for the
# rarified data. The result is a mean of the iterations.
#' @importFrom DelayedMatrixStats colMeans2
.alpha_rarefaction <- function(
        x, assay.type, n.iter, rarefaction.depth, FUN, index, name, ...){
    # Calculating the mean of the subsampled alpha estimates ans storing them
    res <- lapply(seq(n.iter), function(i){
        # Subsampling the counts from the original tse object
        x_sub <- subsampleCounts(
            x, assay.type = assay.type, min_size = rarefaction.depth,
            verbose = FALSE)
        # Calculating the diversity indices on the subsampled object
        x_sub <- do.call(FUN, args = list(
            x_sub, assay.type = assay.type, index = index,
            name = "rarefaction_temp_result", list(...)))
        # Get results
        res <- x_sub[["rarefaction_temp_result"]]
        names(res) <- colnames(x_sub)
        return(res)
    })
    # Combine results from multiple iterations
    res <- do.call(rbind, res)
    # Calculate mean of iterations
    res <- colMeans2(res)
    # It might be that certain samples were dropped off if they have lower
    # abundance than rarefaction  depth --> order so that data includes all the
    # samples
    res <- res[match(colnames(x), names(res))]
    res <- unname(res)
    # Add to original data
    colData(x)[[name]] <- res
    return(x)
}
