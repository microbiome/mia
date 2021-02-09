#' Calculate the Jensen-Shannon Divergence
#'
#' This function calculates the Jensen-Shannon Divergence (JSD) in a
#' \code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#' object.
#'
#' @param x a numeric matrix or a
#'   \code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}.
#'
#' @param exprs_values a single \code{character} value for specifying which
#'   assay to use for calculation.
#'
#' @param transposed Logical scalar, is x transposed with cells in rows?
#'
#' @param ... optional arguments not used.
#'
#' @return a sample-by-sample distance matrix, suitable for NMDS, etc.
#'
#' @seealso
#' \url{http://en.wikipedia.org/wiki/Jensen-Shannon_divergence}
#'
#' @references
#' Jensen-Shannon Divergence and Hilbert space embedding.
#' Bent Fuglede and Flemming Topsoe University of Copenhagen,
#' Department of Mathematics
#' \url{http://www.math.ku.dk/~topsoe/ISIT2004JSD.pdf}
#'
#' @name calculateJSD
#'
#' @author
#' Susan Holmes \email{susan@@stat.stanford.edu}.
#' Adapted for phyloseq2 by Paul J. McMurdie.
#' Adapted for mia by Felix G.M. Ernst
#'
#' @export
#'
#' @examples
#' data(enterotype)
#' calculateJSD(enterotype)
#' enterotype <- runMDS2(enterotype, FUN = calculateJSD, name = "JSD")
#' reducedDim(enterotype)
NULL

setGeneric("calculateJSD", signature = c("x"),
           function(x, ...)
             standardGeneric("calculateJSD"))

#' @rdname calculateJSD
#' @export
setMethod("calculateJSD", signature = c(x = "ANY"),
    function(x){
        calculateDistance(x, FUN = runJSD)
    }
)

#' @rdname calculateJSD
#'
#' @importFrom SummarizedExperiment assay
#'
#' @export
setMethod("calculateJSD", signature = c(x = "SummarizedExperiment"),
    function(x, exprs_values = "counts", transposed = FALSE, ...){
        mat <- assay(x, exprs_values)
        if(!transposed){
            mat <- t(mat)
        }
        calculateJSD(mat, ...)
    }
)

.JSD <- function(x, y){
    # Function to compute Shannon-Jensen Divergence
    # x and y are the frequencies for the same p categories
    # Assumes relative abundance transformation already happened (for efficiency)

    # Define the mean point
    m <- (x+y)/2
    # Define each samples component
    P1 <- x*log(x/m)
    P2 <- y*log(y/m)
    # In the case of zeroes entries log is undefined, JSD is defined as zero
    P1[!is.finite(P1)] <- 0
    P2[!is.finite(P2)] <- 0
    d <- (P1+P2)/2
    return(rowSums(d, na.rm = TRUE))
}

#' @rdname calculateJSD
#'
#' @importFrom utils combn
#' @importFrom stats as.dist
#'
#' @export
runJSD <- function(x){
    # input check
    if(is.null(rownames(x))){
        rownames(x) <- seq_len(nrow(x))
    }
    # Coerce to relative abundance by sample (row)
    x <- sweep(x, 1L, rowSums(x), "/")
    # create N x 2 matrix of all pairwise combinations of samples.
    spn <- utils::combn(rownames(x), 2, simplify = FALSE)
    #
    A <- vapply(spn,"[",character(1),1L)
    B <- vapply(spn,"[",character(1),2L)
    distlist <- .JSD(x[A,], x[B,])
    # reformat
    # initialize distmat with NAs
    distmat <- matrix(NA_real_, nrow(x), nrow(x))
    rownames(distmat) <- colnames(distmat) <- rownames(x)
    matIndices <- matrix(c(B, A), ncol = 2)
    distmat[matIndices] <- distlist
    #
    stats::as.dist(distmat)
}
