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
#' @param assay.type A single character value for selecting the
#'   \code{\link[SummarizedExperiment:SummarizedExperiment-class]{assay}}
#'   to calculate the overlap.
#'   
#' @param assay_name a single \code{character} value for specifying which
#'   assay to use for calculation.
#'   (Please use \code{assay.type} instead. At some point \code{assay_name}
#'   will be disabled.)
#'   
#' @param detection A single numeric value for selecting detection threshold for 
#'   absence/presence of features. Feature that has abundance under threshold in
#'   either of samples, will be discarded when evaluating overlap between samples. 
#'   
#' @param ... Optional arguments not used.
#'   
#' @return getOverlap returns sample-by-sample distance matrix. 
#'   addOverlap returns \code{x} that includes overlap matrix in its 
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
#'   \code{\link[mia:getJSD]{getJSD}}
#'   \code{\link[mia:getUnifrac]{getUnifrac}}
#' 
#' 
#' @name getOverlap
#' @export
#'
#' @author Leo Lahti and Tuomas Borman. Contact: \url{microbiome.github.io}
#' 
#' @examples
#' data(esophagus)
#' tse <- esophagus
#' tse <- transformAssay(tse, method = "relabundance")
#' overlap <- getOverlap(tse, assay_name = "relabundance")
#' overlap
#' 
#' # Store result to reducedDim
#' tse <- addOverlap(tse, assay.type = "relabundance", name = "overlap_between_samples")
#' head(reducedDims(tse)$overlap_between_samples)
#' 
NULL


#' @rdname getOverlap
#' @export
setGeneric("getOverlap", signature = c("x"),
           function(x, assay.type = assay_name, assay_name = "counts", 
                    detection = 0, ...)
             standardGeneric("getOverlap"))

#' @rdname getOverlap
#' @export
setMethod("getOverlap", signature = c(x = "SummarizedExperiment"),
    function(x, assay.type = assay_name, assay_name = "counts", 
             detection = 0, ...){
        ############################# INPUT CHECK ##############################
        # Check assay.type
        .check_assay_present(assay.type, x)
        # Check detection
        if (!.is_numeric_string(detection)) {
          stop("'detection' must be a single numeric value or coercible to ",
               "one.",
               call. = FALSE)
        }
        detection <- as.numeric(detection)
        ########################### INPUT CHECK END ############################
        # Get assay
        res <- getOverlap()
        return(res)
    }
)

#' @rdname getOverlap
#' @export
setMethod("getOverlap", signature = c(x = "ANY"),
    function(x,detection = 0, ...){
        ############################# INPUT CHECK ##############################
        # Check detection
        if (!.is_numeric_string(detection)) {
          stop("'detection' must be a single numeric value or coercible to ",
               "one.",
               call. = FALSE)
        }
        detection <- as.numeric(detection)
        ########################### INPUT CHECK END ############################
        x <- t(x)
        # All the sample pairs
        sample_pairs <- as.matrix(expand.grid(colnames(x), colnames(x)))
        
        # Loop through all sample pairs
        result <- apply(sample_pairs, 1, FUN = function(sample_pair){
          # Get samples
          sample1 <- x[ , sample_pair[1]]
          sample2 <- x[ , sample_pair[2]]
          # Calculate overlap
          temp_result <- .calculate_overlap(sample1, 
                                            sample2, 
                                            detection)
        })
        # Create a matrix from result vector and give name to rownames and colnames
        result <- matrix(result, ncol = ncol(x))
        colnames(result) <- colnames(x)
        rownames(result) <- colnames(x)
        
        # Convert into distances
        result <- stats::as.dist(result)
        return(result)
    }
)

#' @rdname getOverlap
#' @export
setGeneric("addOverlap", signature = c("x"),
           function(x, ...)
               standardGeneric("addOverlap"))

#' @rdname getOverlap
#' 
#' @param name A single character value specifying the name of overlap matrix that
#' is stored in reducedDim(x).
#'   
#' @export
#' @importFrom SingleCellExperiment reducedDim<-
setMethod("addOverlap", signature = c(x = "SummarizedExperiment"),
    function(x, name = "overlap", ...){
        # Check name
        if(!.is_non_empty_string(name)){
            stop("'name' must be a non-empty single character value.",
                call. = FALSE)
        }
        # Calculate overlap
        mat <- getOverlap(x, ...)
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
