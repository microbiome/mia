#' Estimate alpha indices using rarefaction
#' 
#' The function estimates alpha diversity measures optionally using n rounds of rarefaction,
#' given the rarefaction depth, then stores results at \code{\link{colData}}.
#' 
#' @param x a \code{\link{SummarizedExperiment}} object.
#' 
#' @param assay.type the name of the assay used for
#'   calculation of the sample-wise estimates.
#'   
#' @param assay_name a single \code{character} value for specifying which
#'   assay to use for calculation.
#'   (Please use \code{assay.type} instead. At some point \code{assay_name}
#'   will be disabled.)
#'   
#' @param index a \code{character} vector, specifying the alpha diversity measures
#'   to be calculated
#'   
#' @param name a name for the column(s) of the colData the results should be
#'   stored in. By default this will use the original names of the calculated
#'   indices specifying the alpha diversity measures used.
#'   
#' @param ... optional arguments.
#' 
#' 
#' @param rarify logical scalar: Should the alpha diversity measures be estimated
#'   using rarefaction? (default: \code{FALSE})
#' 
#' @param seed a single \code{integer} value as the seed used for the nround
#'  rarefaction.
#' 
#' @param nrounds a single \code{integer} value for the number of rarefaction
#' rounds.
#' 
#' @param rarefaction_depth a \code{double} value as for the minimim size or 
#' rarefaction_depth. (default: \code{min(colSums(assay(x, "counts")), na.rm = TRUE)})
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
#' colData(tse)$shannon_diversity
#'
#'# Calculate observed richness with 10 rarefaction rounds
#' tse <- estimateAlpha(tse, assay.type = "counts", index = "observed_richness",
#' rarify=TRUE, nrounds=10)
#' 
#' # Shows the estimated observed richness
#' colData(tse)$observed_richness
#' 
#' @importFrom dplyr %>% 
#' 
#' @rdname estimateAlpha
#' @export
estimateAlpha <- function(x, assay.type = "counts", assay_name = NULL,
                          index = c("coverage_diversity", "fisher_diversity",
                                    "faith_diversity", "faith",
                                    "gini_simpson_diversity", "inverse_simpson_diversity",
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
                          rarify=FALSE,
                          seed = 123,
                          nrounds=10,
                          rarefaction_depth=min(colSums(assay(x, "counts")), na.rm = TRUE)){
    if(!.is_non_empty_string(index)) {
        stop("'index' should be a non empty string.",
             call. = FALSE)
    }
    if(!.is_a_bool(rarify)){
        stop("'rarify' must be TRUE or FALSE.", call. = FALSE)
    }
    if(!.is_an_integer(seed)) {
        stop("'seed' must be an interger.",
             call. = FALSE)
    }
    if(!.is_an_integer(nrounds)) {
        stop("'nrounds' must be an integer.",
             call. = FALSE)
    }
    if(!(is.double(rarefaction_depth) & rarefaction_depth > 0)) {
        stop("'rarefaction_depth' must be a non-zero positive double.",
             call. = FALSE)
    }
    FUN <- NULL
    if(index %in% .get_indices("diversity")) {
        name <- .parse_name(index, name, "diversity")
        index <- gsub("_diversity", "", index)
        FUN <- .estimate_diversity
    } else if(index %in% .get_indices("dominance")) {
        name <-  .parse_name(index, name, "dominance")
        index <- gsub("_dominance", "", index)
        FUN <- .estimate_dominance
    } else if (index %in% .get_indices("evenness")) {
        name <- .parse_name(index, name, "evenness")
        if (index!="simpson_evenness") {
            index <- gsub("_evenness", "", index)
        }
        FUN <- .estimate_evenness
    } else if (index %in% .get_indices("richness")) {
        name <- .parse_name(index, name, "richness")
        index <- gsub("_richness", "", index)
        FUN <- .estimate_richness
    } else {
        stop("'index' is coresponding to none of the alpha diversity measures.",
             call. = FALSE)
    }
    
    if (rarify) {
        .alpha_rarefaction(x, nrounds = nrounds, seed = seed,
                           args.sub = list(assay.type=assay.type,
                                           min_size=rarefaction_depth,
                                           verbose=FALSE),
                           FUN=FUN,
                           args.fun=list(index=index,
                                         assay.type="subsampled",
                                         ...),
                           name=name)
    } else {
        suppressWarnings(do.call(FUN, list(x, assay.type=assay.type, assay_name=assay_name,
                          index=index, name=name, ...)))
    }

}

## Helper functions

.get_indices <- function(measure) {
    switch(measure,
           "diversity" = c("coverage_diversity", "coverage",
                           "faith_diversity", "faith",
                           "fisher_diversity", "fisher",
                           "gini_simpson_diversity", "gini_simpson",
                           "inverse_simpson_diversity", "inverse_simpson",
                           "log_modulo_skewness_diversity", "log_modulo_skewness",
                           "shannon_diversity", "shannon"),
           "dominance" = c("absolute_dominance", "absolute",
                           "dbp_dominance", "dbp",
                           "core_abundance_dominance", "core_abundance",
                           "gini_dominance", "gini", 
                           "dmn_dominance", "dmn",
                           "relative_dominance", "relative",
                           "simpson_lambda_dominance", "simpson_lambda"),
           "evenness" = c("camargo_evenness", "camargo",
                          "pielou_evenness", "pielou",
                          "simpson_evenness",
                          "evar_evenness", "evar",
                          "bulla_evenness", "bulla"),
           "richness" = c("ace_richness", "ace",
                          "chao1_richness", "chao1",
                          "hill_richness", "hill",
                          "observed_richness", "observed"))
}

.alpha_rarefaction <- function(x,
                               nrounds=1L,
                               seed=123,
                               args.sub=list(assay.type="counts",
                                             min_size=min(colSums(assay(x, "counts")),
                                                          na.rm = TRUE),
                                             verbose=FALSE),
                               FUN=estimateDiversity,
                               args.fun=list(index="shannon",
                                             assay.type="subsampled",
                                             ...),
                               name = args.fun$index) {
    set.seed(seed)
    colData(x)[, name] <- lapply(seq(nrounds), function(i){
        x_sub <- do.call(subsampleCounts, append(list(x), args.sub))
        suppressWarnings(x_sub <- do.call(FUN, append(list(x_sub), args.fun)))
        colData(x_sub)[, args.fun$index, drop=FALSE]
    }) %>% as.data.frame() %>% rowMeans() %>% as.data.frame()
    return(x)
}

.parse_name <- function(index, name, measure) {
    # don't change name if defined by user
    if (name==index) {
        if (measure %in% unlist(strsplit(index, "\\_"))) {
            name = index
        } else {
            name = paste0(index, "_", measure)
        }
    } else {
        return(name)
    }
}

.estimate_diversity <- function(x, assay.type = "counts", assay_name = NULL,
                                index = c("coverage", "fisher", "gini_simpson", 
                                          "inverse_simpson", "log_modulo_skewness",
                                          "shannon"),
                                name = index, ...) {
    estimateDiversity(x, assay.type=assay.type, assay_name=assay_name,
                      index=index, name=name, ...)
}

.estimate_dominance <- function(x,
                                assay.type = assay_name, assay_name = "counts",
                                index = c("absolute", "dbp", "core_abundance",
                                          "gini", "dmn",  "relative",
                                          "simpson_lambda"),
                                ntaxa = 1,
                                aggregate = TRUE,
                                name = index,
                                ...) {
    estimateDominance(x, assay.type=assay.type, assay_name=assay_name,
                      index=index, ntaxa=ntaxa, aggregate=aggregate,
                      name=name, ...)
}

.estimate_evenness <- function(x, assay.type = assay_name, assay_name = "counts",
                               index = c("camargo", "pielou", "simpson_evenness",
                                         "evar", "bulla"),
                               name = index, ...) {
    estimateEvenness(x, assay.type = assay.type, assay_name = assay_name,
                     index=index, name=name, ...)
}

.estimate_richness <- function(x,
                               assay.type = assay_name, assay_name = "counts",
                               index = c("ace", "chao1", "hill", "observed"),
                               name = index,
                               detection = 0,
                               ...) {
    estimateRichness(x, assay.type = assay.type, assay_name = assay_name,
                     index=index, name=name, detection=detection, ...)
}