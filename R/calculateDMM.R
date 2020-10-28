#' Dirichlet-Multinomial Mixture Model Machine Learning for Microbiome Data
#' 
#' These function are accessors for functions implemented in the 
#' \code{\link[DirichletMultinomial:DirichletMultinomial-package]{DirichletMultinomial}} package
#' 
#' @param x a numeric matrix with samples as rows or a 
#'   \code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#'   object.
#'   
#' @param exprs_values a single \code{character} value for specifying which
#'   assay to use for calculation.
#'   
#' @param BPPARAM A
#'   \code{\link[BiocParallel:BiocParallelParam-class]{BiocParallelParam}}
#'   object specifying whether the UniFrac calculation should be parallelized.
#'   
#' @param transposed Logical scalar, is x transposed with samples in rows?
#'   
#' @param ... optional arguments not used.
#' 
#' @name calculateDMN
#' 
#' @examples 
#' fl <- system.file(package="DirichletMultinomial", "extdata", "Twins.csv")
#' counts <- as.matrix(read.csv(fl, row.names=1))
#' se <- SummarizedExperiment(assays = list(counts = counts))
#' dmn <- calculateDMN(se)
#' dmn[[1L]]
NULL

setGeneric("calculateDMN", signature = c("x"),
           function(x, ...)
               standardGeneric("calculateDMN"))

# setGeneric("calculateDMNgroup", signature = c("x"),
#            function(x, ...)
#                standardGeneric("calculateDMNgroup"))
# 
# setGeneric("calculateCvDMNgroup", signature = c("x"),
#            function(x, ...)
#                standardGeneric("calculateCvDMNgroup"))

#' @importFrom DirichletMultinomial dmn
.calculate_DMN <- function(x, k = 1, BPPARAM = SerialParam(), ...){
    if(!is.numeric(k) || 
       length(k) == 0 ||
       anyNA(k) || 
       any(k <= 0) ||
       any(k != as.integer(k))){
        stop("'k' must be an integer vector with positive values only.",
             call. = FALSE)
    }
    #
    old <- getAutoBPPARAM()
    setAutoBPPARAM(BPPARAM)
    on.exit(setAutoBPPARAM(old))
    if (!(bpisup(BPPARAM) || is(BPPARAM, "MulticoreParam"))) {
        bpstart(BPPARAM)
        on.exit(bpstop(BPPARAM), add = TRUE)
    }
    
    ans <- BiocParallel::bplapply(k, DirichletMultinomial::dmn, count = x, ...,
                                  BPPARAM = BPPARAM)
    ans
}

#' @rdname calculateDMN
#' @export
setMethod("calculateDMN", signature = c(x = "ANY"), .calculate_DMN)

#' @rdname calculateDMN
#' @export
setMethod("calculateDMN", signature = c(x = "SummarizedExperiment"),
    function(x, exprs_values = "counts", transposed = FALSE, ...){
        mat <- assay(x, exprs_values)
        if(!transposed){
            mat <- t(mat)
        }
        calculateDMN(mat, ...)
    }
)