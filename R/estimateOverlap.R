#' Estimate overlap
#'
#' @param x a
#'   \code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#'   object containing a tree.
#'   
#' @param abund_values A single character value for selecting the
#'   \code{\link[SummarizedExperiment:SummarizedExperiment-class]{assay}}
#'   to calculate the overlap.
#'   
#' @param detection Detection threshold for absence/presence.
#'   
#' @return a list of matrices
#'
#' @name estimateOverlap
#' @export
#'
#' @author Leo Lahti and Tuomas Borman. Contact: \url{microbiome.github.io}
#' 
#' @examples
#' data(esophagus)
#' tse <- esophagus
#' overlap <- estimateOverlap(tse)
#' head(overlap, 5)
#' 
NULL


#' @rdname estimateOverlap
#' @export
setGeneric("estimateOverlap", signature = c("x"),
           function(x, abund_values = "counts", detection = 0, ...)
             standardGeneric("estimateOverlap"))

#' @rdname estimateOverlap
#' @export
setMethod("estimateOverlap", signature = c(x = "SummarizedExperiment"),
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
        # Create empty matrix samples * samples
        result <- matrix(NA, nrow = ncol(assay), ncol = ncol(assay))
        # Loop through all samples, but not last one
        for (i in seq_len(ncol(assay)-1) ) {
          # Loop through all samples that have larger index than i --> all 
          # combinations are went through only once
          for (j in seq(i+1, ncol(assay)) ) {
            # Get samples
            x <- assay[ , i]
            y <- assay[ , j]
            # Store overlap of two samples to coordinates i,j and j,i
            result[i, j] <- result[j, i] <- .calculate_overlap(x, y, detection)
          }
        }
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





