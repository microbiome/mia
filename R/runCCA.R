#' Canonical Correspondance Analysis
#' 
#' ToDo
#'
#' @name runCCA
#' @seealso
#' 
#' @examples
NULL

setGeneric("calculateCCA", signature = c("x"),
           function(x, ...)
               standardGeneric("calculateCCA"))

setGeneric("runCCA", signature = c("x"),
           function(x, ...)
               standardGeneric("runCCA"))


.calculate_cca <- function(x, cov, ncomponents = 50, ntop = 500,
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
    
}

#' @export
#' @rdname runDPCoA
setMethod("calculateCCA", "ANY", .calculate_cca)


#' @export
#' @rdname runCCA
setMethod("calculateCCA", "MicrobiomeExperiment",
          function(x, ..., exprs_values = "relabundance", dimred = NULL,
                   n_dimred = NULL)
          {
              mat <- .get_mat_from_sce(x, exprs_values = exprs_values,
                                       dimred = dimred, n_dimred = n_dimred)
              cov <- data.frame()
              .calculate_cca(mat, cov, transposed = !is.null(dimred), ...)
          }
)

#' @export
#' @rdname runCCA
#' @importFrom SingleCellExperiment reducedDim<-
setMethod("runCCA", "MicrobiomeExperiment",
          function(x, ..., altexp=NULL, name="PCA")
          {
              if (!is.null(altexp)) {
                  y <- altExp(x, altexp)
              } else {
                  y <- x
              }
              reducedDim(x, name) <- calculateCCA(y, ...)
              x
          }
)