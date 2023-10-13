#' Estimate alpha indices using rarefaction
#' 
#' The function estimates alpha diversity indices using n rounds of rarefaction,
#' then stores results at \code{\link{colData}}.
#' 
#' @param x a \code{\link{SummarizedExperiment}} object.
#' 
#' @param nrounds a single \code{integer} value for the number of rarefaction
#' rounds.
#' 
#' @param seed a single \code{integer} value that creates the seeds for the 
#' nround rarefaction.
#' 
#' @param args.sub argument list passed to \code{\link[mia:subsampleCounts]{subsampleCounts}}
#' 
#' @param FUN the alpha diversity function to be used; e.g.
#'  \code{\link[mia:estimateDiversity]{estimateDiversity}},
#'  \code{\link[mia:estimateEvenness]{estimateEvenness}}, 
#'  \code{\link[mia:estimateRichness]{estimateRichness}}.
#' 
#' @param args.fun argument list passed to the alpha diversity function \code{FUN}
#' 
#' @param name The column name where to place results at \code{\link{colData}}. 
#' 
#' @return \code{x} with additional \code{\link{colData}} named after the index 
#' used.
#' 
#' @examples
#' 
#' data("GlobalPatterns")
#' tse <- GlobalPatterns
#' 
#' # Calculate the default Shannon index with 1 rarefaction round
#' tse <- estimateAlphaWithRarefaction(tse)
#' 
#' # Shows the estimated Shannon index
#' colData(tse)$shannon
#'
#'# Calculate the default observed richness with 10 rarefaction rounds
#' tse <- estimateAlphaWithRarefaction(tse, nrounds=10,
#'  FUN=mia::estimateRichness, args.fun=list(index="observed"))
#' 
#' # Shows the estimated observed richness
#' colData(tse)$richness
#' 
#' @importFrom dplyr %>% 
#' 
#' @rdname estimateAlphaWithRarefaction
#' @export
estimateAlphaWithRarefaction <- function(x,
                                         nrounds=1L,
                                         seed=123,
                                         args.sub=list(assay.type="counts",
                                                       min_size=min(colSums(assay(x, "counts")), na.rm = TRUE),
                                                       verbose=FALSE),
                                         FUN=mia::estimateDiversity,
                                         args.fun=list(index="shannon",
                                                       assay.type="subsampled"),
                                         name = args.fun$index){
    # checks
    if(!.is_an_integer(nrounds)) {
        stop("'nrounds' must be an interger.",
             call. = FALSE)
        }
    if(!.is_an_integer(seed)) {
        stop("'seed' must be an interger.",
             call. = FALSE)
    }
    if(!.is_non_empty_string(name)) {
        stop("'name' should be a non empty string.",
             call. = FALSE)
    }
    
    # Generating seeds for every round
    set.seed(seed)
    SEEDS <- sample.int(10000, size = nrounds)
    colData(x)[, name] <- lapply(seq(nrounds), function(i){
        x_sub <- do.call(subsampleCounts, append(list(x, seed = SEEDS[i]),
                                                 args.sub))
        x_sub <- do.call(FUN, append(list(x_sub), args.fun))
        colData(x_sub)[, name, drop=FALSE]
        }) %>% as.data.frame() %>% rowMeans() %>% as.data.frame()
    return(x)
}
