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
#' @param BPPARAM A
#'   \code{\link[BiocParallel:BiocParallelParam-class]{BiocParallelParam}}
#'   object specifying whether the JSD calculation should be parallelized.
#'
#' @param chunkSize a integer scalar, defining the size of data send
#'   to the individual worker. Only has an effect, if \code{BPPARAM} defines
#'   more than one worker. (default: \code{chunkSize = nrow(x)})
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
#' Adapted for phyloseq by Paul J. McMurdie.
#' Adapted for mia by Felix G.M. Ernst
#'
#' @export
#'
#' @examples
#' data(enterotype)
#' jsd <- calculateJSD(enterotype)
#' class(jsd)
#' head(jsd)
#' enterotype <- runMDS(enterotype, FUN = calculateJSD, name = "JSD")
#' head(reducedDim(enterotype))
#' head(attr(reducedDim(enterotype),"eig"))
#' attr(reducedDim(enterotype),"GOF")
NULL

setGeneric("calculateJSD", signature = c("x"),
           function(x, ...)
             standardGeneric("calculateJSD"))

#' @rdname calculateJSD
#' @export
setMethod("calculateJSD", signature = c(x = "ANY"),
    function(x, ...){
        calculateDistance(x, FUN = runJSD, ...)
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

# written by Susan Holmes \email{susan@@stat.stanford.edu}.
# Adapted for phyloseq by Paul J. McMurdie.
# Adapted for mia by Felix G.M. Ernst
#' @importFrom DelayedArray rowSums
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
#' @importFrom BiocParallel SerialParam register bplapply bpisup bpstart bpstop
#' @importFrom DelayedArray getAutoBPPARAM setAutoBPPARAM
#'
#' @export
runJSD <- function(x, BPPARAM = SerialParam(), chunkSize = nrow(x)){
    # input check
    if(is.null(rownames(x))){
        rownames(x) <- seq_len(nrow(x))
    }
    if(missing(chunkSize) || is.na(chunkSize) || is.null(chunkSize) ||
       !is.integer(chunkSize)){
        chunkSize <- nrow(x)
    } else if(length(chunkSize) != 1L) {
        chunkSize <- chunkSize[1L]
    }
    #
    old <- getAutoBPPARAM()
    setAutoBPPARAM(BPPARAM)
    on.exit(setAutoBPPARAM(old))
    if (!(bpisup(BPPARAM) || is(BPPARAM, "MulticoreParam"))) {
        bpstart(BPPARAM)
        on.exit(bpstop(BPPARAM), add = TRUE)
    }
    # Coerce to relative abundance by sample (row)
    x <- sweep(x, 1L, rowSums(x), "/")
    # create N x 2 matrix of all pairwise combinations of samples.
    spn <- utils::combn(rownames(x), 2, simplify = TRUE)
    #
    N <- ncol(spn)
    f <- ceiling(seq_len(N)/chunkSize)
    A <- split(spn[1L,], f)
    B <- split(spn[2L,], f)
    FUN <- function(X, a, b){
        .JSD(X[a,,drop=FALSE], X[b,,drop=FALSE])
    }
    distlist <- BiocParallel::bpmapply(FUN, A, B,
                                       MoreArgs = list(X = x),
                                       BPPARAM = BPPARAM,
                                       SIMPLIFY = FALSE)
    distlist <- do.call(c, unname(distlist))
    # reformat
    # initialize distmat with NAs
    distmat <- matrix(NA_real_, nrow(x), nrow(x))
    rownames(distmat) <- colnames(distmat) <- rownames(x)
    matIndices <- matrix(c(unlist(B), unlist(A)), ncol = 2)
    distmat[matIndices] <- distlist
    #
    stats::as.dist(distmat)
}
