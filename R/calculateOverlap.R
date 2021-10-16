#' Estimate overlap
#' 
#' This function calculates overlap for all sample-pairs
#' in a \code{\link[TreeSummarizedExperiment:TreeSummarizedExperiment-class]{TreeSummarizedExperiment}}
#' object.
#'
#' @param x a
#'   \code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#'   object containing a tree.
#'   
#' @param abund_values A single character value for selecting the
#'   \code{\link[SummarizedExperiment:SummarizedExperiment-class]{assay}}
#'   to calculate the overlap.
#'   
#' @param detection A single numeric value for selecting detection threshold for 
#'   absence/presence of features. Feature that has abundance under threshold in
#'   either of samples, will be discarded when evaluating overlap between samples. 
#'   
#' @param ... Optional arguments not used.
#'   
#' @return A sample-by-sample distance matrix.
#' 
#' @details This function calculates overlap between all the sample-pairs. Overlap
#'   reflects similarity between sample-pairs. 
#'   
#'   When overlap is calculated using relative abundances, the higher the value the 
#'   higher the similarity is, When using relative abundances, overlap value 1 means that 
#'   all the abundances of features are equal between two samples, and 0 means that 
#'   samples have completely different relative abundances. 
#'
#' @seealso
#'   \code{\link[mia:calculateJSD]{calculateJSD}}
#'   \code{\link[mia:calculateUniFrac]{calculateUniFrac}}
#' 
#' 
#' @name calculateOverlap
#' @export
#'
#' @author Leo Lahti and Tuomas Borman. Contact: \url{microbiome.github.io}
#' 
#' @examples
#' data(esophagus)
#' tse <- esophagus
#' tse <- transformSamples(tse, method = "relabundance")
#' overlap <- calculateOverlap(tse, abund_values = "relabundance")
#' overlap
#' 
NULL


#' @rdname calculateOverlap
#' @export
setGeneric("calculateOverlap", signature = c("x"),
           function(x, abund_values = "counts", detection = 0, ...)
             standardGeneric("calculateOverlap"))

#' @rdname calculateOverlap
#' @export
setMethod("calculateOverlap", signature = c(x = "SummarizedExperiment"),
    function(x, abund_values = "counts", detection = 0, ...){
        ############################# INPUT CHECK ##############################
        # Check abund_values
        .check_assay_present(abund_values, x)
        # Check detection
        if (!.is_numeric_string(detection)) {
          stop("'detection' must be a single numeric value or coercible to ",
               "one.",
               call. = FALSE)
        }
        detection <- as.numeric(detection)
        ########################### INPUT CHECK END ############################
        # Get assay
        assay <- assay(x, abund_values)
        # Create empty matrices samples * samples
        result <- matrix(NA, nrow = ncol(assay), ncol = ncol(assay))
        
        # Loop through all samples, but not last one
        for (i in seq_len(ncol(assay)-1) ) {
          # Loop through all samples that have larger index than i --> all 
          # combinations are went through only once
          for (j in seq(i+1, ncol(assay)) ) {
            # Get samples
            sample1 <- assay[ , i]
            sample2 <- assay[ , j]
            
            # Galculate overlap
            temp_result <- .calculate_overlap(sample1, 
                                              sample2, 
                                              detection)
            
            # Store values of two samples to coordinates i,j and j,i
            result[i, j] <- result[j, i] <- temp_result
          }
        }
        # Change sample names
        colnames(result) <- colnames(x)
        rownames(result) <- colnames(x)
        # Convert into distances
        result <- stats::as.dist(result)
        return(result)
    }
)

################################ HELP FUNCTIONS ################################

.calculate_overlap <- function (x, y, detection) {
  # Take those taxa that have abundance over threshold
  inds <- which(x > detection & y > detection)
  x <- x[inds]
  y <- y[inds]
  # Overlap is the average of the sums of the values in each sample
  overlap <- (sum(x) + sum(y))/2
  return(overlap)
}
