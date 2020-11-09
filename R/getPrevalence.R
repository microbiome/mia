#' @title Taxa Prevalence
#' 
#' @name getPrevalence
#' 
#' @description This function estimates population prevalence for microbial taxa in a \code{\link{SummarizedExperiment-class}} object.
#' 
#' @param x \code{\link{SummarizedExperiment-class}} object
#' 
#' @param detection Detection threshold for absence/presence
#' (strictly greater by default; in 0 to 1).
#' 
#' @param include_lowest Include the lower boundary of the detection and
#' prevalence cutoffs. FALSE by default.
#' 
#' @param sort Sort the groups by prevalence (default: FALSE)
#' 
#' @param as_relative Logical. Apply detection threshold on compositional (relative) abundances.
#' 
#' @param abund_values A single character value for selecting the
#'   \code{\link[SummarizedExperiment:SummarizedExperiment-class]{assay}}
#'   to use for calculating the relative abundance.
#' 
#' @param ... Arguments to pass.
#' 
#' @details Calculates the frequency (between 0 and 1) of samples that exceed the
#' detection threshold. For SummarizedExperiment objects, calculates this for each
#' taxa. The absolute population prevalence can be obtained by multiplying the
#' prevalence by the number of samples.
#' 
#' @return For each taxa, the frequency of samples where a given taxa is
#' detected at a given detection threshold. The output is provided as a
#' frequency (between 0 and 1).
#' 
#' @references 
#' A Salonen et al. The adult intestinal core microbiota is determined by 
#' analysis depth and health status. Clinical Microbiology and Infection 
#' 18(S4):16 20, 2012. 
#' To cite the R package, see citation('mia')
#' 
#' @author Leo Lahti
#' 
#' @keywords utilities
#' 
#' @export
#' 
#' @examples
#'
#' # Get prevalence estimates for microbial taxa (population frequencies)
#' data(GlobalPatterns)     
#' prevalence.frequency <- getPrevalence(GlobalPatterns, detection=0, sort=TRUE, as_relative=TRUE)
#' print(head(prevalence.frequency))
#'
#' # Get prevalence estimates for microbial taxa (population counts)
#' # - the getPrevalence function itself always returns population frequencies
#' # - to obtain population counts, multiply frequencies with the sample size:
#' prevalence.count <- prevalence.frequency * ncol(GlobalPatterns)


#' @rdname getPrevalence
#' @export
setGeneric("getPrevalence", signature = "x",
           function(x, ...)
               standardGeneric("getPrevalence"))

#' @rdname getPrevalence
#' @export
setMethod("getPrevalence", signature = c(x = "ANY"),
    function(x, detection, include_lowest, sort, as_relative, ...){

        if (!is.numeric(detection)) {
            stop("The detection argument in getPrevalence function should be numeric.")
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
    function(x, detection=0, include_lowest=FALSE, sort=FALSE, as_relative=TRUE, abund_values = "counts", ...){

        if (as_relative) { 
            # Add relative abundances if not yet available
            x <- relAbundanceCounts(x)
            abund_values <- "relabundance"
        } 

        # check assay
        .check_abund_values(abund_values, x)

        # retrieve abundance matrix
        x <- assay(x, abund_values)
    
        getPrevalence(x, detection=detection,
            include_lowest=include_lowest,
            sort=sort,
            as_relative=as_relative, ...)
    
    }
)



