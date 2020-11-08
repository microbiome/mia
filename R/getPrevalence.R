#' @title Taxa Prevalence 
#' @description Simple prevalence measure.
#' @param x A vector, data matrix or \code{\link{MicrobiomeExperiment}} object
#' @param sort Sort the groups by prevalence
#' @param as_relative Logical. Indicate prevalence as fraction of samples
#' (in percentage between 0 to 1 default); or in absolute counts indicating
#' the number of samples where the OTU is detected (strictly) above the given
#' abundance threshold.
#' @inheritParams prevalentMembers
#' @details For vectors, calculates the fraction (as_relative=FALSE) or
#' number (as_relative=TRUE) of samples that exceed the
#' detection. For matrices, calculates this for each matrix
#' column. For MicrobiomeExperiment objects, calculates this for each OTU. The
#' relative prevalence (as_relative=FALSE) is simply the absolute
#' prevalence (as_relative=TRUE) divided by the number of samples.
#' @return For each OTU, the fraction of samples where a given OTU is
#' detected. The output is readily given as a percentage.
#' @references 
#' A Salonen et al. The adult intestinal core microbiota is determined by 
#' analysis depth and health status. Clinical Microbiology and Infection 
#' 18(S4):16 20, 2012. 
#' To cite the R package, see citation('mia') 
#' @author Leo Lahti
#' @keywords utilities
#' @export
#' @examples
#' data(enterotype)     
#' pr <- getPrevalence(enterotype, detection=0, sort=TRUE, as_relative=TRUE)
#' pr <- getPrevalence(enterotype, detection=0, sort=TRUE, as_relative=FALSE)
getPrevalence <- function(x, detection=0, sort=FALSE, as_relative=FALSE,
    include.lowest=FALSE) {
    
    if (is.null(detection)) {
        detection <- (-Inf)
    }
    
    if (is.null(x)) {
        warning("x is NULL - returning NULL")
        return(NULL)
    }

    # Add relative abundances if not yet available
    if (!"relabundance" %in% names(x@assays)) {
        x <- relAbundanceCounts(x)
    }

    # Convert to matrix
    x <- assays(x)$relabundance
    
    if (is.vector(x)) {
    
        if (include.lowest) {
            prev <- sum(x >= detection)
        } else {
            prev <- sum(x > detection)
        }
   
    } else if (is.matrix(x) || is.data.frame(x)) {
        
        if (include.lowest) {
            prev <- rowSums(x >= detection)
        } else {
            prev <- rowSums(x > detection)
        }
    }
    
    if (!as_relative) {
        prev <- prev/prevalence_nsamples(x)
    }
    
    if (sort) {
        prev <- rev(sort(prev))
    }
    
    # Return
    prev
    
}




# Internal auxiliary function
prevalence_nsamples <- function(x) {
    
    if (is.vector(x)) {
        n <- length(x)
    } else if (is.matrix(x) || is.data.frame(x)) {
        n <- ncol(x)
    }
    
    n
    
}
