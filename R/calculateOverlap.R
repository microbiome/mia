#' Estimate overlap
#' 
#' This function calculates overlap for all sample-pairs
#' in a \code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#' object.
#'
#' @param x a
#'   \code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#'   object containing a tree.
#'   
#' @param assay_name A single character value for selecting the
#'   \code{\link[SummarizedExperiment:SummarizedExperiment-class]{assay}}
#'   to calculate the overlap.
#'   
#' @param abund_values a single \code{character} value for specifying which
#'   assay to use for calculation.
#'   (Please use \code{assay_name} instead. At some point \code{abund_values}
#'   will be disabled.)
#'   
#' @param detection A single numeric value for selecting detection threshold for 
#'   absence/presence of features. Feature that has abundance under threshold in
#'   either of samples, will be discarded when evaluating overlap between samples. 
#'   
#' @param ... Optional arguments not used.
#'   
#' @return calculateOverlap returns sample-by-sample distance matrix. 
#'   runOverlap returns \code{x} that includes overlap matrix in its 
#'   reducedDim. 
#' 
#' @details These function calculates overlap between all the sample-pairs. Overlap
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
#' overlap <- calculateOverlap(tse, assay_name = "relabundance")
#' overlap
#' 
#' # Store result to reducedDim
#' tse <- runOverlap(tse, assay_name = "relabundance", name = "overlap_between_samples")
#' head(reducedDims(tse)$overlap_between_samples)
#' 
NULL


#' @rdname calculateOverlap
#' @export
setGeneric("calculateOverlap", signature = c("x"),
           function(x, assay_name = abund_values, abund_values = "counts", 
                    detection = 0, ...)
             standardGeneric("calculateOverlap"))

#' @rdname calculateOverlap
#' @export
setMethod("calculateOverlap", signature = c(x = "SummarizedExperiment"),
    function(x, assay_name = abund_values, abund_values = "counts", 
             detection = 0, ...){
        ############################# INPUT CHECK ##############################
        # Check assay_name
        .check_assay_present(assay_name, x)
        # Check detection
        if (!.is_numeric_string(detection)) {
          stop("'detection' must be a single numeric value or coercible to ",
               "one.",
               call. = FALSE)
        }
        detection <- as.numeric(detection)
        ########################### INPUT CHECK END ############################
        # Get assay
        assay <- assay(x, assay_name)
        
        # All the sample pairs
        sample_pairs <- as.matrix(expand.grid(colnames(x), colnames(x)))
        
        # Loop through all sample pairs
        result <- apply(sample_pairs, 1, FUN = function(sample_pair){
          # Get samples
          sample1 <- assay[ , sample_pair[1]]
          sample2 <- assay[ , sample_pair[2]]
          # Calculate overlap
          temp_result <- .calculate_overlap(sample1, 
                                            sample2, 
                                            detection)
        })
        # Create a matrix from result vector and give name to rownames and colnames
        result <- matrix(result, ncol = ncol(assay))
        colnames(result) <- colnames(assay)
        rownames(result) <- colnames(assay)
        
        # Convert into distances
        result <- stats::as.dist(result)
        return(result)
    }
)

#' @rdname calculateOverlap
#' @export
setGeneric("runOverlap", signature = c("x"),
           function(x, ...)
               standardGeneric("runOverlap"))

#' @rdname calculateOverlap
#' 
#' @param name A single character value specifying the name of overlap matrix that
#' is stored in reducedDim(x).
#'   
#' @export
#' @importFrom SingleCellExperiment reducedDim<-
setMethod("runOverlap", signature = c(x = "SummarizedExperiment"),
    function(x, name = "overlap", ...){
        # Check name
        if(!.is_non_empty_string(name)){
            stop("'name' must be a non-empty single character value.",
                call. = FALSE)
        }
        # Calculate overlap
        mat <- calculateOverlap(x, ...)
        # Convert it into matrix so that nrow equals number of samples
        mat <- as.matrix(mat)
        # Store it to reducedDim
        reducedDim(x, name) <- mat
        return(x)
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
