#' Calculation of Double Principal Correspondance analysis
#' 
#' ToDo
#'
#' @name runDPCoA
#' @seealso
#' 
#' @examples
NULL

setGeneric("calculateDPCoA", signature = c("x"),
           function(x, ...)
               standardGeneric("calculateDPCoA"))

setGeneric("runDPCoA", signature = c("x"),
           function(x, ...)
               standardGeneric("runDPCoA"))


.calculate_dpcoa <- function(x, tree, ncomponents = 50, ntop = 500,
                             subset_row = NULL, scale = FALSE, transposed=FALSE,
                             BPPARAM = SerialParam())
{
    old <- getAutoBPPARAM()
    setAutoBPPARAM(BPPARAM)
    on.exit(setAutoBPPARAM(old))
    if (!(bpisup(BPPARAM) || is(BPPARAM, "MulticoreParam"))) {
        bpstart(BPPARAM)
        on.exit(bpstop(BPPARAM), add = TRUE)
    }
    
    if (!transposed) {
        out <- .get_mat_for_reddim(x, subset_row=subset_row, ntop=ntop, scale=scale, get.var=TRUE) 
        x <- out$x
        cv <- out$v
    } else {
        cv <- colVars(DelayedArray(x))
    }
    
    pca <- runPCA(x, rank=ncomponents, BSPARAM=BSPARAM, BPPARAM=BPPARAM)
    varExplained <- pca$sdev^2
    percentVar <- varExplained / sum(cv) * 100
    
    # Saving the results
    pcs <- pca$x
    rownames(pcs) <- rownames(x)
    attr(pcs, "varExplained") <- varExplained
    attr(pcs, "percentVar") <- percentVar
    rownames(pca$rotation) <- colnames(x)
    attr(pcs, "rotation") <- pca$rotation
    pcs
}

#' @export
#' @rdname runDPCoA
setMethod("calculateDPCoA", "ANY", .calculate_dpcoa)


#' @export
#' @rdname runDPCoA
setMethod("calculateDPCoA", "MicrobiomeExperiment",
    function(x, ..., exprs_values = "logcounts", dimred = NULL, n_dimred = NULL)
    {
        mat <- .get_mat_from_sce(x, exprs_values = exprs_values,
                                 dimred = dimred, n_dimred = n_dimred)
        .calculate_dpcoa(mat, transposed = !is.null(dimred), ...)
    }
)

#' @export
#' @rdname runDPCoA
#' @importFrom SingleCellExperiment reducedDim<-
setMethod("runDPCoA", "MicrobiomeExperiment",
    function(x, ..., altexp=NULL, name="PCA")
    {
        if (!is.null(altexp)) {
            y <- altExp(x, altexp)
        } else {
            y <- x
        }
        reducedDim(x, name) <- calculateDPCoA(y, ...)
        x
    }
)