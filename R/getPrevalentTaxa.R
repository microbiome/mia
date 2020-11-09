#' @title Prevalent Taxa
#' 
#' @name getPrevalentTaxa
#' 
#' @description This function retrieves prevalent taxa from a \code{\link{SummarizedExperiment-class}} object.
#' 
#' @inheritParams getPrevalence
#' 
#' @param prevalence Prevalence threshold (in 0 to 1). The
#' required prevalence is strictly greater by default. To include the
#' limit, set include_lowest to TRUE. 
#' 
#' @return Vector of prevalent taxa
#' 
#' @details Lists taxa that are more prevalent with the
#' given detection threshold based on relative abundances or absolute counts
#' (see as_relative argument). 
#' 
#' @examples
#' data(GlobalPatterns)
#' # Detection threshold 1 (strictly greater by default);
#' # Note that the data (GlobalPatterns) is here in absolute counts
#' # (and not compositional, relative abundances)
#' # Prevalence threshold 50 percent (strictly greater by default)
#' taxa <- getPrevalentTaxa(GlobalPatterns, detection=1/100, prevalence=50/100, as_relative=TRUE)
#' print(head(taxa))
#' 
#' @export
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

#' @rdname getPrevalentTaxa
#' @export
setGeneric("getPrevalentTaxa", signature = "x",
           function(x, ...)
               standardGeneric("getPrevalentTaxa"))

#' @rdname getPrevalentTaxa
#' @export
setMethod("getPrevalentTaxa", signature = c(x = "SummarizedExperiment"),
    function(x, detection=1/100, prevalence=50/100,
    include_lowest=FALSE, sort=FALSE, as_relative=TRUE, ...){

        pr <- getPrevalence(x, detection=detection, sort=sort,
            as_relative=as_relative, include_lowest=include_lowest)

        if (include_lowest) {
    
            taxa <- names(which(pr >= prevalence))
        
        } else {
    
            taxa <- names(which(pr > prevalence))
        
        }
    
        return(taxa)
    
    }
)



