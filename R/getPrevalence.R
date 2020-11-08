#' @title Taxa Prevalence
#' @name getPrevalence 
#' @description Simple prevalence measure.
#' @param x A vector, data matrix or \code{\link{MicrobiomeExperiment}} object
#' @inheritParams getPrevalentTaxa
#' @details Calculates the frequency (between 0 and 1) of samples that exceed the
#' detection threshold. For MicrobiomeExperiment objects, calculates this for each
#' taxa. The absolute population prevalence can be obtained by multiplying the
#' prevalence by the number of samples.
#' @return For each taxa, the frequency of samples where a given taxa is
#' detected at a given detection threshold. The output is provided as a
#' frequency (between 0 and 1).
#' @references 
#' A Salonen et al. The adult intestinal core microbiota is determined by 
#' analysis depth and health status. Clinical Microbiology and Infection 
#' 18(S4):16 20, 2012. 
#' To cite the R package, see citation('mia') 
#' @author Leo Lahti
#' @keywords utilities
#' @export
#' @examples
#' data(GlobalPatterns)     
#' pr <- getPrevalence(GlobalPatterns, detection=0, sort=TRUE, as_relative=TRUE)
#' pr <- getPrevalence(GlobalPatterns, detection=0, sort=TRUE, as_relative=FALSE)

#' @rdname getPrevalence
#' @export
setGeneric("getPrevalence", signature = "x",
           function(x, ...)
               standardGeneric("getPrevalence"))

#' @rdname getPrevalence
#' @export
setMethod("getPrevalence", signature = c(x = "ANY"),
    function(x, detection, include_lowest, sort, as_relative, ...){
    
        # do vector or matrix stuff here
        if (is.null(detection)) {
            detection <- (-Inf)
        }
    
        if (is.null(x)) {
            warning("x is NULL - returning NULL")
            return(NULL)
        } 

        if (include_lowest) {
            prev <- rowSums(x >= detection)
        } else {
            prev <- rowSums(x > detection)
        }

        # Always return prevalence as a relative frequency.
	# This helps to avoid confusion with detection limit, which
	# is applied on either relative abundances or absolute counts
	# with as_relative argument
        prev <- prev/ncol(x)
    
        if (sort) {
            prev <- rev(sort(prev))
        }
    
        # Return
        prev

    }
)

#' @rdname getPrevalence
#' @export
setMethod("getPrevalence", signature = c(x = "SummarizedExperiment"),
    function(x, detection=0, include_lowest=FALSE, sort=FALSE, as_relative=FALSE, ...){

        if (as_relative) { 

            # Add relative abundances if not yet available
            x <- relAbundanceCounts(x)

            # Convert to matrix	
            x <- relabundance(x)

        } else {

            x <- assay(x, "counts")

        }

        getPrevalence(x, detection, include_lowest, sort, as_relative, ...)
	
    }
)



