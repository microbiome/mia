#' Estimate alpha diversity
#'
#' Several functions for calculation of alpha diversity indices available via
#' wrapper functions. They are implemented via the \code{vegan} package.
#'
#' This includes the \sQuote{Shannon}, \sQuote{Simpson} and
#' \sQuote{inverse Simpson} diversity measures.
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
#'   (Currently not used)
#'
#' @param ... additional parameters passed to \code{estimateDiversity}
#'
#' @return \code{x} with additional \code{\link{colData}} named
#'   \code{*name*}
#'
#' @seealso
#' \code{\link[scater:plotColData]{plotColData}}
#' \itemize{
#'   \item{\code{\link[vegan:diversity]{diversity}}}
#'   \item{\code{\link[vegan:specpool]{estimateR}}}
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
#' plots <- lapply(c("shannon","simpson"),
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
                    index = c("shannon","simpson","inv_simpson", "richness",
                              "chao1", "ACE"),
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
setMethod("estimateDiversity", signature = c(x = "SummarizedExperiment"),
    function(x, abund_values = "counts",
             index = c("shannon","simpson","inv_simpson", "richness", "chao1",
                       "ACE"),
             name = index, ..., BPPARAM = SerialParam()){
        # input check
        index<- match.arg(index, several.ok = TRUE)
        if(!.is_non_empty_character(name) || length(name) != length(index)){
            stop("'name' must be a non-empty character value and have the ",
                 "same length then 'index'.",
                 call. = FALSE)
        }
        .check_abund_values(abund_values, x)
        .require_package("vegan")
        #
        dvrsts <- lapply(index,
                         .run_dvrsty,
                         x = x,
                         mat = assay(x, abund_values),
                         BPPARAM = BPPARAM, ...)
        .add_dvrsty_values_to_colData(x, dvrsts, name)
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
        estimateDiversity(x, index = "richness",  ...)
    }
)

.get_shannon <- function(x, ...){
    vegan::diversity(t(x), index="shannon")
}

.get_simpson <- function(x, ...){
    vegan::diversity(t(x), index="simpson")
}

.get_inverse_simpson <- function(x, ...){
    vegan::diversity(t(x), index="invsimpson")
}

.get_observed <- function(x, ...){
    vegan::estimateR(t(x))["S.obs",]
}

.get_chao1 <- function(x, ...){
    ans <- t(vegan::estimateR(t(x))[c("S.chao1","se.chao1"),])
    colnames(ans) <- c("","se")
    ans
}

.get_ACE <- function(x, ...){
    ans <- t(vegan::estimateR(t(x))[c("S.ACE","se.ACE"),])
    colnames(ans) <- c("","se")
    ans
}

#' @importFrom SummarizedExperiment assay assays
.run_dvrsty <- function(x, i, mat, BPPARAM = SerialParam(), ...){
    dvrsty_FUN <- switch(i,
                         shannon = .get_shannon,
                         simpson = .get_simpson,
                         inv_simpson = .get_inverse_simpson,
                         richness = .get_observed,
                         chao1 = .get_chao1,
                         ACE = .get_ACE)
    dvrsty <- dvrsty_FUN(mat, BPPARAM = BPPARAM, ...)
    dvrsty
}

#' @importFrom SummarizedExperiment colData colData<-
#' @importFrom S4Vectors DataFrame
.add_dvrsty_values_to_colData <- function(x, dvrsts, name){
    dvrsts <- mapply(
        function(dvrsty, n){
            dvrsty <- DataFrame(dvrsty)
            colnames(dvrsty)[1L] <- n
            if(ncol(dvrsty) > 1L){
                i <- seq.int(2,ncol(dvrsty))
                colnames(dvrsty)[i] <- paste0(n,"_",colnames(dvrsty)[i])
            }
            dvrsty
        },
        dvrsts,
        name)
    colData(x) <- cbind(colData(x),DataFrame(dvrsts))
    x
}
