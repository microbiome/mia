#' Estimate alpha diversity
#'
#' Several functions for calculation of alpha diversity indices available via 
#' wrapper functions. They are implemented via the \code{breakaway} package.
#' 
#' This includes the \sQuote{Shannon},  \sQuote{Simpson}, 
#'
#' @param x a \code{\link{MicrobiomeExperiment}} object
#'
#' @param abund_values the name of the assay used for calculation of the
#'   sample-wise estimates
#'
#' @param name a name for the column of the colData the results should be stored
#'   in.
#'
#' @param BPPARAM A
#'   \code{\link[BiocParallel:BiocParallelParam-class]{BiocParallelParam}}
#'   object specifying whether calculation of estimates should be parallelized.
#'
#' @param ... additional parameters passed to the function implemented by the
#'   \code{breakaway} package
#'
#' @return \code{x} with additional \code{\link{colData}} named
#'   \code{*name*}
#'
#' @seealso 
#' \itemize{
#'   \item{\code{\link[=estimateBreakway]{estimateBreakway}}}
#'   \item{\code{\link[breakaway:sample_inverse_simpson]{sample_inverse_simpson}}}
#'   \item{\code{\link[breakaway:sample_simpson]{sample_simpson}}}
#'   \item{\code{\link[breakaway:sample_richness]{sample_richness}}}
#'   \item{\code{\link[breakaway:sample_shannon]{sample_shannon}}}
#'   \item{\code{\link[breakaway:sample_shannon_e]{sample_shannon_e}}}
#' }
#'
#' @name estimateDiversity
#'
#' @examples
#' data(esophagus, package = "MicrobiomeExperiment")
#' esophagus <- estimateShannon(esophagus)
#' colData(esophagus)$Shannon
#'
#' esophagus <- estimateSimpson(esophagus)
#' colData(esophagus)$Simpson
NULL

#' @rdname estimateDiversity
#' @export
setGeneric("estimateShannon",signature = c("x"),
           function(x, abund_values = "counts", name = "Shannon", ...)
               standardGeneric("estimateShannon"))

#' @rdname estimateDiversity
#' @export
setGeneric("estimateSimpson",signature = c("x"),
           function(x, abund_values = "counts", name = "Simpson", ...)
               standardGeneric("estimateSimpson"))

#' @rdname estimateDiversity
#' @export
setGeneric("estimateInvSimpson",signature = c("x"),
           function(x, abund_values = "counts", name = "InvSimpson", ...)
               standardGeneric("estimateInvSimpson"))

#' @rdname estimateDiversity
#' @export
setGeneric("estimateRichness",signature = c("x"),
           function(x, abund_values = "counts", name = "Richness", ...)
               standardGeneric("estimateRichness"))

#' @rdname estimateDiversity
#' @export
setGeneric("estimateShannonE",signature = c("x"),
           function(x, abund_values = "counts", name = "ShannonE", ...)
               standardGeneric("estimateShannonE"))

#' @rdname estimateDiversity
#' @export
setMethod("estimateShannon", signature = c(x = "MicrobiomeExperiment"),
    function(x, abund_values = "counts", name = "Shannon",
             BPPARAM = SerialParam(), ...){
        .run_brkwy_dvrsty(x = x, abund_values = abund_values, name = name,
                          FUN = breakaway::sample_shannon, BPPARAM = BPPARAM, ...)
    }
)

#' @rdname estimateDiversity
#' @export
setMethod("estimateSimpson", signature = c(x = "MicrobiomeExperiment"),
    function(x, abund_values = "counts", name = "Simpson",
             BPPARAM = SerialParam(), ...){
        .run_brkwy_dvrsty(x = x, abund_values = abund_values, name = name,
                          FUN = breakaway::sample_simpson, BPPARAM = BPPARAM, ...)
    }
)

#' @rdname estimateDiversity
#' @export
setMethod("estimateInvSimpson", signature = c(x = "MicrobiomeExperiment"),
          function(x, abund_values = "counts", name = "InvSimpson",
                   BPPARAM = SerialParam(), ...){
              .run_brkwy_dvrsty(x = x, abund_values = abund_values, name = name,
                                FUN = breakaway::sample_inverse_simpson,
                                BPPARAM = BPPARAM, ...)
          }
)

#' @rdname estimateDiversity
#' @export
setMethod("estimateRichness", signature = c(x = "MicrobiomeExperiment"),
          function(x, abund_values = "counts", name = "Richness",
                   BPPARAM = SerialParam(), ...){
              .run_brkwy_dvrsty(x = x, abund_values = abund_values, name = name,
                                FUN = breakaway::sample_richness,
                                BPPARAM = BPPARAM, ...)
          }
)

#' @rdname estimateDiversity
#' @export
setMethod("estimateShannonE", signature = c(x = "MicrobiomeExperiment"),
          function(x, abund_values = "counts", name = "ShannonE",
                   BPPARAM = SerialParam(), ...){
              .run_brkwy_dvrsty(x = x, abund_values = abund_values, name = name,
                                FUN = breakaway::sample_shannon_e,
                                BPPARAM = BPPARAM, ...)
          }
)


#' @importFrom SummarizedExperiment assay assays
.run_brkwy_dvrsty <- function(x, abund_values = "counts", name, FUN,
                               BPPARAM = SerialParam(), ...){
    # input check
    .require_package("breakaway")
    # input checks
    if(!.is_non_empty_string(abund_values)){
        stop("'abund_values' must be a non empty single character value.",
             call. = FALSE)
    }
    if(!(abund_values %in% names(assays(x)))){
        stop("'abund_values' must reference a name of an assay in 'x'")
    }
    if(!.is_non_empty_string(name)){
        stop("'name' must be a single non-empty character value.",
             call. = FALSE)
    }
    
    brkwy_values <- .run_brkwy_dvrsty_on_assay(assay(x, abund_values), FUN,
                                               BPPARAM = BPPARAM, ...)
    .add_brkwy_dvrsty_values_to_colData(x, brkwy_values, name, abund_values)
}

#' @importFrom S4Vectors DataFrame
#' @importFrom IRanges NumericList
.run_brkwy_dvrsty_on_assay <- function(mat, FUN, BPPARAM, ...){
    val <- apply(mat, 2, breakaway::convert)
    
    #
    old <- getAutoBPPARAM()
    setAutoBPPARAM(BPPARAM)
    on.exit(setAutoBPPARAM(old))
    if (!(bpisup(BPPARAM) || is(BPPARAM, "MulticoreParam"))) {
        bpstart(BPPARAM)
        on.exit(bpstop(BPPARAM), add = TRUE)
    }
    #
    estimates <- BiocParallel::bplapply(val, FUN, ..., BPPARAM=BPPARAM)
    ans <- vapply(estimates, "[[", numeric(1),"estimate")
    ans
}

#' @importFrom SummarizedExperiment colData colData<-
.add_brkwy_dvrsty_values_to_colData <- function(x, brkwy_values, name,
                                                abund_values){
    colData <- colData(x)
    colData[[name]] <- brkwy_values
    colData(x) <- colData
    x
}
