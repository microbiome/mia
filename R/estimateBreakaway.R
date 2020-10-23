#' Functions for diversity estimation from \code{breakaway} package
#'
#' Several functions for calculation of diversity indices implemented by the 
#' \code{breakaway} package are available via wrapper functions.
#' 
#' Other diversity measures implemented by breakaway are also available,
#' such as \code{\link[=estimateDiversity]{Shannon}}.
#'
#' @param x a \code{\link[=MicrobiomeExperiment-class]{MicrobiomeExperiment}}
#'   object
#'
#' @param abund_values the name of the assay used for calculation of the
#'   sample-wise estimates
#'
#' @param name a name for the column of the \code{colData} the results should be
#'   stored in.
#'
#' @param BPPARAM A
#'   \code{\link[BiocParallel:BiocParallelParam-class]{BiocParallelParam}}
#'   object specifying whether calculation of estimates should be parallelized.
#'
#' @param ... additional parameters passed to the function implemented by the
#'   \code{breakaway} package
#'
#' @return \code{x} with additional \code{\link{colData}} named
#'   \code{*name*} containing the main diversity score and 
#'   \code{*name*_summary} containing a summary \code{DataFrame}
#'
#' @seealso
#' \itemize{
#'   \item{\code{\link[=estimateDiversity]{estimateDiversity}}}
#'   \item{\code{\link[breakaway:breakaway]{breakaway}}}
#'   \item{\code{\link[breakaway:kemp]{kemp}}}
#'   \item{\code{\link[breakaway:wlrm_transformed]{wlrm_transformed}}}
#'   \item{\code{\link[breakaway:wlrm_untransformed]{wlrm_untransformed}}}
#'   \item{\code{\link[breakaway:chao_bunge]{chao_bunge}}}
#'   \item{\code{\link[breakaway:poisson_model]{poisson_model}}}
#'   \item{\code{\link[breakaway:chao1]{chao1}}}
#' }
#'
#' @name estimateBreakway
#'
#' @examples
#' data(esophagus, package = "MicrobiomeExperiment")
#' esophagus <- estimateBreakaway(esophagus)
#' colData(esophagus)$breakaway
#' colData(esophagus)$breakaway_summary
#'
#' esophagus <- estimateWLRMuntransformed(esophagus)
#' colData(esophagus)$WLRMun
#' colData(esophagus)$WLRMun_summary
NULL

#' @rdname estimateBreakway
#' @export
setGeneric("estimateBreakaway",signature = c("x"),
           function(x, abund_values = "counts", name = "breakaway", ...)
             standardGeneric("estimateBreakaway"))

#' @rdname estimateBreakway
#' @export
setGeneric("estimateKemp",signature = c("x"),
           function(x, abund_values = "counts", name = "kemp", ...)
             standardGeneric("estimateKemp"))

#' @rdname estimateBreakway
#' @export
setGeneric("estimateWLRMuntransformed",signature = c("x"),
           function(x, abund_values = "counts", name = "WLRMun", ...)
             standardGeneric("estimateWLRMuntransformed"))

#' @rdname estimateBreakway
#' @export
setGeneric("estimateWLRMtransformed",signature = c("x"),
           function(x, abund_values = "counts", name = "wlrm", ...)
             standardGeneric("estimateWLRMtransformed"))

#' @rdname estimateBreakway
#' @export
setGeneric("estimatePoissonModel",signature = c("x"),
           function(x, abund_values = "counts", name = "poisson", ...)
             standardGeneric("estimatePoissonModel"))

#' @rdname estimateBreakway
#' @export
setGeneric("estimateChaoBunge",signature = c("x"),
           function(x, abund_values = "counts", name = "chao_bunge", ...)
             standardGeneric("estimateChaoBunge"))

#' @rdname estimateBreakway
#' @export
setGeneric("estimateChao1",signature = c("x"),
           function(x, abund_values = "counts", name = "chao1", ...)
             standardGeneric("estimateChao1"))

#' @rdname estimateBreakway
#' @export
setMethod("estimateBreakaway", signature = c(x = "MicrobiomeExperiment"),
    function(x, abund_values = "counts", name = "breakaway",
             BPPARAM = SerialParam(), ...){
        .run_breakaway_fun(x = x, abund_values = abund_values, name = name,
                           FUN = breakaway::breakaway, BPPARAM = BPPARAM, ...)
    }
)

#' @rdname estimateBreakway
#' @export
setMethod("estimateKemp", signature = c(x = "MicrobiomeExperiment"),
    function(x, abund_values = "counts", name = "kemp",
             BPPARAM = SerialParam(), ...){
        .run_breakaway_fun(x = x, abund_values = abund_values, name = name,
                           FUN = breakaway::kemp, BPPARAM = BPPARAM, ...)
    }
)

#' @rdname estimateBreakway
#' @export
setMethod("estimateWLRMuntransformed",
    signature = c(x = "MicrobiomeExperiment"),
    function(x, abund_values = "counts", name = "WLRMun",
             BPPARAM = SerialParam(), ...){
        .run_breakaway_fun(x = x, abund_values = abund_values, name = name,
                           FUN = breakaway::wlrm_untransformed,
                           BPPARAM = BPPARAM, ...)
    }
)

#' @rdname estimateBreakway
#' @export
setMethod("estimateWLRMtransformed",
    signature = c(x = "MicrobiomeExperiment"),
    function(x, abund_values = "counts", name = "WLRM",
             BPPARAM = SerialParam(), ...){
        .run_breakaway_fun(x = x, abund_values = abund_values, name = name,
                           FUN = breakaway::wlrm_transformed,
                           BPPARAM = BPPARAM, ...)
    }
)

#' @rdname estimateBreakway
#' @export
setMethod("estimatePoissonModel", signature = c(x = "MicrobiomeExperiment"),
    function(x, abund_values = "counts", name = "poisson",
             BPPARAM = SerialParam(), ...){
        .run_breakaway_fun(x = x, abund_values = abund_values, name = name,
                           FUN = breakaway::poisson_model, BPPARAM = BPPARAM,
                           ...)
    }
)

#' @rdname estimateBreakway
#' @export
setMethod("estimateChaoBunge", signature = c(x = "MicrobiomeExperiment"),
    function(x, abund_values = "counts", name = "chao_bunge",
             BPPARAM = SerialParam(), ...){
        .run_breakaway_fun(x = x, abund_values = abund_values, name = name,
                           FUN = breakaway::chao_bunge, BPPARAM = BPPARAM, ...)
    }
)

#' @rdname estimateBreakway
#' @export
setMethod("estimateChao1", signature = c(x = "MicrobiomeExperiment"),
    function(x, abund_values = "counts", name = "chao1",
             BPPARAM = SerialParam(), ...){
        .run_breakaway_fun(x = x, abund_values = abund_values, name = name,
                           FUN = breakaway::chao1, BPPARAM = BPPARAM, ...)
    }
)

#' @importFrom SummarizedExperiment assay assays
.run_breakaway_fun <- function(x, abund_values = "counts", name, FUN,
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
    #
    brkwy_values <- .run_brkwy_on_assay(assay(x, abund_values), FUN,
                                        BPPARAM = BPPARAM, ...)
    .add_brkwy_values_to_colData(x, brkwy_values, name, abund_values)
}

#' @importFrom S4Vectors DataFrame
#' @importFrom IRanges NumericList
.run_brkwy_on_assay <- function(mat, FUN, BPPARAM, ...){
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

    estimate <- vapply(estimates, "[[", numeric(1),"estimate")
    error <- vapply(estimates, "[[", numeric(1),"error")
    interval <- lapply(estimates, "[[", "interval")
    model <- vapply(estimates, "[[", character(1),"model")
    warnings <- lapply(estimates, "[[", "warnings")
    reasonable <- vapply(estimates, "[[", logical(1),"reasonable")

    ans <- DataFrame(estimate = estimate,
                     error = error,
                     interval = NumericList(interval),
                     model = model,
                     reasonable = reasonable)
    ans$warnings <- warnings
    ans
}

#' @importFrom SummarizedExperiment colData colData<-
.add_brkwy_values_to_colData <- function(x, brkwy_values, name, abund_values){
    colData <- colData(x)
    colData[[name]] <- brkwy_values$estimate
    colData[[paste0(name,"_summary")]] <- brkwy_values
    colData(x) <- colData
    x
}
