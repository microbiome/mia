#' Estimate alpha diversity indices.
#' 
#' The function estimates alpha diversity indices optionally using of rarefaction,
#' then stores results at \code{\link{colData}}.
#' 
#' @param x a \code{\link{SummarizedExperiment}} object.
#' 
#' @param assay.type the name of the assay used for
#'   calculation of the sample-wise estimates.
#'   
#' @param index a \code{character} vector, specifying the alpha diversity indices
#'   to be calculated
#'   
#' @param name a name for the column(s) of the colData the results should be
#'   stored in. By default this will use the original names of the calculated
#'   indices.
#'   
#' @param ... optional arguments.
#' 
#' @param n.iter a single \code{integer} value for the number of rarefaction
#' rounds.
#' 
#' @param rarefaction.depth a \code{double} value as for the minimim size or 
#' rarefaction.depth. (default: \code{min(colSums(assay(x, "counts")), na.rm = TRUE)})
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
#' tse <- estimateAlpha(tse, assay.type = "counts", index = "shannon")
#' 
#' # Shows the estimated Shannon index
#' tse$shannon_diversity
#'
#'# Calculate observed richness with 10 rarefaction rounds
#' tse <- estimateAlpha(tse, assay.type = "counts", index = "observed_richness",
#' rarefaction.depth=min(colSums(assay(tse, "counts")), na.rm = TRUE), n.iter=10)
#' 
#' # Shows the estimated observed richness
#' tse$observed_richness
#' 
#' @importFrom dplyr %>% 
#' 

#' @rdname estimateAlpha
#' @export
setGeneric("estimateAlpha",signature = c("x"),
           function(x,
                    assay.type = "counts", 
                    index = c("coverage_diversity", "fisher_diversity",
                              "faith_diversity", "gini_simpson_diversity",
                              "inverse_simpson_diversity",
                              "log_modulo_skewness_diversity", "shannon_diversity",
                              "absolute_dominance", "dbp_dominance",
                              "core_abundance_dominance", "gini_dominance",
                              "dmn_dominance", "relative_dominance",
                              "simpson_lambda_dominance",
                              "camargo_evenness", "pielou_evenness",
                              "simpson_evenness", "evar_evenness",
                              "bulla_evenness",
                              "ace_richness", "chao1_richness", "hill_richness",
                              "observed_richness"),
                    name = index,
                    ...,
                    n.iter=10,
                    rarefaction.depth=max(colSums(assay(x, assay.type)), na.rm = TRUE))
               standardGeneric("estimateAlpha"))

#' @rdname estimateAlpha
#' @export
setMethod("estimateAlpha", signature = c(x = "SummarizedExperiment"),
          function(x,
                   assay.type = "counts",
                   index = c("coverage_diversity", "fisher_diversity",
                             "faith_diversity", "gini_simpson_diversity",
                             "inverse_simpson_diversity",
                             "log_modulo_skewness_diversity", "shannon_diversity",
                             "absolute_dominance", "dbp_dominance",
                             "core_abundance_dominance", "gini_dominance",
                             "dmn_dominance", "relative_dominance",
                             "simpson_lambda_dominance",
                             "camargo_evenness", "pielou_evenness",
                             "simpson_evenness", "evar_evenness",
                             "bulla_evenness",
                             "ace_richness", "chao1_richness", "hill_richness",
                             "observed_richness"),
                   name = index,
                   ...,
                   n.iter=10,
                   rarefaction.depth=max(colSums(assay(x, assay.type)), na.rm = TRUE)){
              # checks
              if(is.null(index) & any(!sapply(index, .is_non_empty_string))) {
                  stop("'index' should be a character vector.", call. = FALSE)
                  }
              if(!.is_an_integer(n.iter)) {
                  stop("'n.iter' must be an integer.", call. = FALSE)
                  }
              if(!(is.numeric(rarefaction.depth) & rarefaction.depth > 0)) {
                  stop("'rarefaction.depth' must be a non-zero positive double.",
                       call. = FALSE)
                  }
              # if multiple indices to be estimated, name should a vector of same length
              if(length(index)!=length(name)) {
                  stop("'index' and 'name' should be vectors of the same length.",
                       call. = FALSE)
                  }
              # Looping over the vector of indices to be estimated
              for (i in seq_along(index)) {
                  # Getting the corresponding alpha indices function by parsing the index
                  FUN <- NULL
                  if(any(grepl(index[i], .get_indices("diversity")))) {
                      # making name having the alpha type suffix or leave it as is if 
                      # user defined
                      name[i] <- .parse_name(index[i], name[i], "diversity")
                      # cleaning index from suffix to be used with the corresponding index
                      # function
                      index[i] <- gsub("_diversity", "", index[i])
                      FUN <- .estimate_diversity
                  } else if (any(grepl(index[i], .get_indices("dominance")))) {
                      name[i] <-  .parse_name(index[i], name[i], "dominance")
                      index[i] <- gsub("_dominance", "", index[i])
                      FUN <- .estimate_dominance
                  } else if (any(grepl(index[i], .get_indices("evenness")))) {
                      name[i] <- .parse_name(index[i], name[i], "evenness")
                      if (index[i]!="simpson_evenness") {
                          index[i] <- gsub("_evenness", "", index[i])
                          }
                      FUN <- .estimate_evenness
                  } else if (any(grepl(index[i], .get_indices("richness")))) {
                      name[i] <- .parse_name(index[i], name[i], "richness")
                      index[i] <- gsub("_richness", "", index[i])
                      FUN <- .estimate_richness
                  } else {
                  stop("'index' is coresponding to none of the alpha diversity indices.
                  'index' should be one of: coverage_diversity, fisher_diversity,
                  faith_diversity, gini_simpson_diversity,
                  inverse_simpson_diversity, log_modulo_skewness_diversity,
                  shannon_diversity, absolute_dominance, dbp_dominance,
                  core_abundance_dominance, gini_dominance,
                  dmn_dominance, relative_dominance, simpson_lambda_dominance,
                  camargo_evenness, pielou_evenness, simpson_evenness,
                  evar_evenness, bulla_evenness, ace_richness, chao1_richness,
                  hill_richness or observed_richness.",
                       call. = FALSE)
                  }
                  # Performing rarefaction if rarefaction.depth is specified to be less 
                  # the max of the total read counts
                  if (rarefaction.depth < max(colSums(assay(x, assay.type)), na.rm = TRUE)) {
                      x <- .alpha_rarefaction(x, n.iter = n.iter,
                                              args.sub = list(assay.type=assay.type,
                                                              min_size=rarefaction.depth,
                                                              verbose=FALSE),
                                              FUN=FUN,
                                              args.fun=list(index=index[i], assay.type="subsampled"),
                                              ...,
                                              name=name[i])
                  } else {
                      # Estimate index without rarefaction
                      # warning is supressed due to the deprication of the functions called.
                      suppressWarnings(x <- do.call(FUN, args = c(list(x, assay.type=assay.type,
                                                                       index=index[i],
                                                                       name=name[i]),
                                                                  list(...))))
                  }
              }
              return(x)
    }
)
## Helper functions

.get_indices <- function(
        measure) {
    switch(measure,
           "diversity" = c("coverage_diversity", "faith_diversity", 
                           "fisher_diversity", "gini_simpson_diversity",
                           "inverse_simpson_diversity",
                           "log_modulo_skewness_diversity", "shannon_diversity"),
           "dominance" = c("absolute_dominance",
                           "dbp_dominance", "core_abundance_dominance",
                           "gini_dominance", "dmn_dominance", "relative_dominance",
                           "simpson_lambda_dominance"),
           "evenness" = c("camargo_evenness", "pielou_evenness", "simpson_evenness",
                          "evar_evenness", "bulla_evenness"),
           "richness" = c("ace_richness", "chao1_richness", "hill_richness", 
                          "observed_richness"))
}

.alpha_rarefaction <- function(
        x,
        n.iter=1L,
        args.sub=list(assay.type="counts", min_size=min(colSums(assay(x, "counts")),
                                                        na.rm = TRUE),
                      verbose=FALSE),
        FUN=.estimate_diversity,
        args.fun=c(index="shannon",
                   assay.type="subsampled"),
        ...,
        name = args.fun$index) {
    # Calculating the mean of the subsampled alpha estimates ans storing them
    colData(x)[, name] <- lapply(seq(n.iter), function(j){
        # subsampling the counts from the original tse object
        x_sub <- do.call(subsampleCounts, args = c(list(x), args.sub))
        # calculating the diversity indices on the subsampled object
        # warnings are supressed due to the depricated warning of the alpha
        # measure functions
        suppressWarnings(x_sub <- do.call(FUN, args = c(list(x_sub),
                                                        args.fun,
                                                        list(...))))
        # Storing estimate results
        colData(x_sub)[, args.fun$index, drop=FALSE]
    }) %>% data.frame(.) %>% rowMeans(.) %>% data.frame(.)
    return(x)
}

.parse_name <- function(
        index, name, measure) {
    # parsing name string to use as a column name at colData when storing estimates
    if (name==index) {
        # check if suffix of the alpha indices if present at index
        # otherwise keeping suffix as a name if name not defined by user.
        if (measure %in% unlist(strsplit(index, "\\_"))) {
            name <- index
        } else {
            name <- paste0(index, "_", measure)
        }
    } else {
        # don't change name if defined by user
        return(name)
    }
}

.estimate_diversity <- function(
        x, assay.type = "counts",
        index = c("coverage", "fisher", "gini_simpson",  "inverse_simpson",
                  "log_modulo_skewness","shannon"),
        name = index, ...) {
    estimateDiversity(x, assay.type=assay.type, index=index, name=name, ...)
}

.estimate_dominance <- function(
        x,
        assay.type = "counts", 
        index = c("absolute", "dbp", "core_abundance", "gini", "dmn",  "relative",
                  "simpson_lambda"),
        ntaxa = 1,
        aggregate = TRUE,
        name = index,
        ...) {
    estimateDominance(x, assay.type=assay.type, index=index, ntaxa=ntaxa,
                      aggregate=aggregate, name=name, ...)
}

.estimate_evenness <- function(
        x, assay.type = "counts",
        index = c("camargo", "pielou", "simpson_evenness", "evar", "bulla"),
        name = index, ...) {
    estimateEvenness(x, assay.type = assay.type, index=index, name=name, ...)
}

.estimate_richness <- function(
        x,
        assay.type = "counts",
        index = c("ace", "chao1", "hill", "observed"),
        name = index,
        detection = 0,
        ...) {
    estimateRichness(x, assay.type = assay.type, index=index, name=name,
                     detection=detection, ...)
}