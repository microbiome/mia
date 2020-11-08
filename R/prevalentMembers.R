#' @title Prevalent Taxa
#' @description Determine members of the prevalent microbiota with a given abundance
#' and prevalence thresholds.
#' @param x \code{\link{MicrobiomeExperiment-class}} object
#' @param detection Detection threshold for absence/presence
#' (strictly greater by default; in [0,1]).
#' @param prevalence Prevalence threshold (in [0, 1]). The
#' required prevalence is strictly greater by default. To include the
#' limit, set include.lowest to TRUE. 
#' @param include.lowest Include the lower boundary of the detection and
#' prevalence cutoffs. FALSE by default.
#' @param ... Arguments to pass.
#' @return Vector of prevalent members
#' @details For a MicrobiomeExperiment object, lists taxa that are more prevalent with the
#' given detection threshold based on relative abundances. 
#' @examples
#' data(enterotype)
#' # Detection threshold 1 (strictly greater by default);
#' # Note that the data (enterotype) is here in absolute counts
#' # (and not compositional, relative abundances)
#' # Prevalence threshold 50 percent (strictly greater by default)
#' a <- prevalentMembers(enterotype, 1/100, 50/100)
#' @export
#' @references 
#' A Salonen et al. The adult intestinal core microbiota is determined by 
#' analysis depth and health status. Clinical Microbiology and Infection 
#' 18(S4):16 20, 2012. 
#' To cite the R package, see citation('mia') 
#' @author Leo Lahti 
#' @keywords utilities
prevalentMembers <- function(x, detection=1/100, prevalence=50/100,
    include.lowest=FALSE, ...) {

    if ((prevalence < 0) | (prevalence > 1)) {
        stop("The prevalence argument should be in [0, 1].")
    }

    if (include.lowest) {
        taxa <- names(which(prevalence(x, detection,
            include.lowest=include.lowest, count=FALSE) >= prevalence))
    } else {
        taxa <- names(which(prevalence(x, detection,
        include.lowest=include.lowest, count=FALSE) > prevalence))
    }
    
    taxa
    
}


