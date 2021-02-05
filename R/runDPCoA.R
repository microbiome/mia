#' Calculation of Double Principal Correspondance analysis
#'
#' ToDo
#'
#' @name runDPCoA
#' @seealso
#'
#' @examples
NULL

setGeneric("calculateDPCoA", signature = c("x", "y"),
           function(x, y, ...)
               standardGeneric("calculateDPCoA"))

.calculate_dpcoa <- function(x, y, ncomponents = 2, ntop = 500,
                             subset_row = NULL, scale = FALSE,
                             transposed = FALSE, ...)
{
    .require_package("ade4")
    # input check
    y <- as.matrix(y)
    if(length(unique(dim(y))) != 1L){
        stop("'y' must be symmetric.", call. = FALSE)
    }
    #
    if(!transposed) {
        x <- .get_mat_for_reddim(x, subset_row = subset_row, ntop = ntop,
                                 scale = scale)
    }
    if(nrow(y) != ncol(x)){
        stop("x and y must have corresponding dimensions.", call. = FALSE)
    }
    y <- y[colnames(x),colnames(x)]
    y <- sqrt(y)
    y <- as.dist(y)
    #
    dpcoa <- ade4::dpcoa(data.frame(x), y, scannf = FALSE, nf = ncomponents)
    ans <- as.matrix(dpcoa$li)
    rownames(ans) <- rownames(x)
    colnames(ans) <- NULL
    attr(ans,"eig") <- dpcoa$eig
    tmp <- as.matrix(dpcoa$dls)
    rownames(tmp) <- colnames(x)
    colnames(tmp) <- NULL
    attr(ans,"sample_red") <- tmp
    attr(ans,"feature_weights") <- unname(dpcoa$dw)
    attr(ans,"sample_weights") <- unname(dpcoa$lw)
    ans
}

#' @export
#' @rdname runDPCoA
setMethod("calculateDPCoA", c("ANY","ANY"), .calculate_dpcoa)

#' @export
#' @importFrom ape cophenetic.phylo
#' @rdname runDPCoA
setMethod("calculateDPCoA", signature = c("TreeSummarizedExperiment","missing"),
    function(x, ..., exprs_values = "logcounts", dimred = NULL, n_dimred = NULL)
    {
        .require_package("ade4")
        mat <- assay(x, exprs_values)
        dist <- cophenetic.phylo(rowTree(x))
        calculateDPCoA(mat, dist, ...)
    }
)

#' @export
#' @rdname runDPCoA
#' @importFrom SingleCellExperiment reducedDim<-
runDPCoA <- function(x, ..., altexp = NULL, name = "DPCoA"){
    if (!is.null(altexp)) {
        y <- altExp(x, altexp)
    } else {
        y <- x
    }
    reducedDim(x, name) <- calculateDPCoA(y, ...)
    x
}
