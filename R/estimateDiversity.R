#' Estimate alpha diversity
#'
#' Several functions for calculation of alpha diversity indices available via
#' wrapper functions. They are implemented via the \code{breakaway} package.
#'
#' This includes the \sQuote{Shannon}, \sQuote{Shannon-E}, \sQuote{Simpson},
#' \sQuote{inverse Simpson} and \sQuote{Richness} diversity measures.
#'
#' @param x a \code{\link{SummarizedExperiment}} object
#'
#' @param abund_values the name of the assay used for calculation of the
#'   sample-wise estimates
#'
#' @param index a diversity measurement
#'
#' @param name a name for the column of the colData the results should be stored
#'   in.
#'
#' @param BPPARAM A
#'   \code{\link[BiocParallel:BiocParallelParam-class]{BiocParallelParam}}
#'   object specifying whether calculation of estimates should be parallelized.
#'
#' @param ... additional parameters passed to \code{estimateDiversity} or the
#'   corresponding function implemented by the \code{breakaway} package.
#'
#' @return \code{x} with additional \code{\link{colData}} named
#'   \code{*name*}
#'
#' @seealso
#' \code{\link[scater:plotColData]{plotColData}}
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
#' data(GlobalPatterns)
#' se <- GlobalPatterns
#'
#' se <- estimateShannon(se)
#' colData(se)$shannon
#'
#' esophagus <- estimateSimpson(se)
#' colData(se)$simpson
#'
#' # calculating all the diversites
#' se <- estimateDiversity(se)
#'
#' # plotting the diversities
#' library(scater)
#' plotColData(se, "shannon")
#' # ... by sample type
#' plotColData(se, "shannon", "SampleType")
#' \donttest{
#' # combining different plots
#' plots <- lapply(c("shannon","simpson","richness"),
#'                 plotColData,
#'                 object = se,
#'                 x = "SampleType",
#'                 colour_by = "SampleType")
#' plots <- lapply(plots,"+", theme(axis.text.x = element_text(angle=45,hjust=1)))
#' ggpubr::ggarrange(plotlist = plots, nrow = 1, common.legend = TRUE, legend = "right")
#' }
NULL

#' @rdname estimateDiversity
#' @export
setGeneric("estimateDiversity",signature = c("x"),
           function(x, abund_values = "counts",
                    index = c("shannon","shannon_e","simpson","inv_simpson","richness"),
                    name = index, ...)
               standardGeneric("estimateDiversity"))

#' @rdname estimateDiversity
#' @export
setGeneric("estimateShannon",signature = c("x"),
           function(x, ...)
               standardGeneric("estimateShannon"))

#' @rdname estimateDiversity
#' @export
setGeneric("estimateSimpson",signature = c("x"),
           function(x, ...)
               standardGeneric("estimateSimpson"))

#' @rdname estimateDiversity
#' @export
setGeneric("estimateInvSimpson",signature = c("x"),
           function(x, ...)
               standardGeneric("estimateInvSimpson"))

#' @rdname estimateDiversity
#' @export
setGeneric("estimateRichness",signature = c("x"),
           function(x, ...)
               standardGeneric("estimateRichness"))

#' @rdname estimateDiversity
#' @export
setGeneric("estimateShannonE",signature = c("x"),
           function(x, ...)
               standardGeneric("estimateShannonE"))


#' @rdname estimateDiversity
#' @export
setMethod("estimateDiversity", signature = c(x = "SummarizedExperiment"),
    function(x, abund_values = "counts",
             index = c("shannon","shannon_e","simpson","inv_simpson","richness"),
             name = index, BPPARAM = SerialParam(), ...){
        # input check
        index<- match.arg(index,
                          c("shannon","shannon_e","simpson","inv_simpson",
                            "richness"),
                          several.ok = TRUE)
        if(!.is_non_empty_character(name) || length(name) != length(index)){
            stop("'name' must be a non-empty character value and have the ",
                 "same length then 'index'.",
                 call. = FALSE)
        }
        #
        FUN <- function(i){
            dvrsty_FUN <- switch(i,
                                 shannon = breakaway::sample_shannon,
                                 shannon_e = breakaway::sample_shannon_e,
                                 simpson = breakaway::sample_simpson,
                                 inv_simpson = breakaway::sample_inverse_simpson,
                                 richness = breakaway::sample_richness)
            .run_brkwy_dvrsty(x = x, abund_values = abund_values,
                              FUN = dvrsty_FUN, BPPARAM = BPPARAM, ...)
        }
        dvrsts <- lapply(index, FUN)
        .add_brkwy_dvrsty_values_to_colData(x, dvrsts, name)
    }
)

#' @rdname estimateDiversity
#' @export
setMethod("estimateShannon", signature = c(x = "SummarizedExperiment"),
    function(x, ...){
        estimateDiversity(x, index = "shannon", ...)
    }
)

#' @rdname estimateDiversity
#' @export
setMethod("estimateSimpson", signature = c(x = "SummarizedExperiment"),
    function(x, ...){
        estimateDiversity(x, index = "simpson", ...)
    }
)

#' @rdname estimateDiversity
#' @export
setMethod("estimateInvSimpson", signature = c(x = "SummarizedExperiment"),
    function(x, ...){
        estimateDiversity(x, index = "inv_simpson",  ...)
    }
)

#' @rdname estimateDiversity
#' @export
setMethod("estimateRichness", signature = c(x = "SummarizedExperiment"),
    function(x, ...){
        estimateDiversity(x, index = "richness", ...)
    }
)

#' @rdname estimateDiversity
#' @export
setMethod("estimateShannonE", signature = c(x = "SummarizedExperiment"),
    function(x, ...){
        estimateDiversity(x, index = "shannon_e", ...)
    }
)


#' @importFrom SummarizedExperiment assay assays
.run_brkwy_dvrsty <- function(x, abund_values = "counts", FUN,
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
    dvrsty <- .run_brkwy_dvrsty_on_assay(assay(x, abund_values), FUN,
                                         BPPARAM = BPPARAM, ...)
    dvrsty
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
#' @importFrom S4Vectors DataFrame
.add_brkwy_dvrsty_values_to_colData <- function(x, dvrsts, name){
    colData <- colData(x)
    colData[,name] <- DataFrame(dvrsts)
    colData(x) <- colData
    x
}
