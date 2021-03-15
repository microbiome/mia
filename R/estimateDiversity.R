#' Estimate diversity measures
#'
#' Several functions for calculation of diversity indices available via
#' wrapper functions. They are implemented via the \code{vegan} package.
#'
#' These include the \sQuote{Shannon}, \sQuote{Gini-Simpson} and
#' \sQuote{inverse Simpson}, \sQuote{coverage}, and \sQuote{Fisher's alpha}
#' diversity measures.
#'
#' @param x a \code{\link{SummarizedExperiment}} object
#'
#' @param abund_values the name of the assay used for calculation of the
#'   sample-wise estimates
#'
#' @param index a \code{character} vector, specifying the diversity measures
#'   to be calculated.
#'
#' @param name a name for the column(s) of the colData the results should be
#'   stored in.
#'
#' @param BPPARAM A
#'   \code{\link[BiocParallel:BiocParallelParam-class]{BiocParallelParam}}
#'   object specifying whether calculation of estimates should be parallelized.
#'
#' @param ... additional parameters passed to \code{estimateDiversity}
#' \itemize{
#'   \item{\code{threshold} is a numeric value for selecting threshold for \code{coverage} index.
#'   By default, the threshold is 0.9.
#'   }
#' }
#'
#' @return \code{x} with additional \code{\link{colData}} named
#'   \code{*name*}
#'
#' @details
#'
#' Diversity is a joint quantity that combines elements or community richness and evenness.
#' Diversity increases, in general, when species richness or evenness increase.
#'
#' By default, this function returns all indices.
#'
#' The available diversity indices include the following:
#' \itemize{
#' \item{inverse_simpson }{Inverse Simpson diversity:
#' $1/lambda$ where $lambda=sum(p^2)$ and $p$ are relative abundances.
#' This corresponds to the diversity index
#' 'invsimpson' in the vegan::diversity. This should not be confused with the
#' closely related Gini-Simpson index}
#'
#' \item{gini_simpson }{Gini-Simpson diversity i.e. $1 - lambda$, where lambda is the
#' Simpson index, calculated as the sum of squared relative abundances.
#' This corresponds to the diversity index
#' 'simpson' in the vegan::diversity.
#' This is also called Gibbs–Martin, or Blau index in sociology,
#' psychology and management studies. The Gini-Simpson index (1-lambda) should not be
#' confused with Simpson's dominance (lambda), Gini index, or inverse Simpson index (1/lambda).}
#'
#' \item{shannon }{Shannon diversity ie entropy}
#'
#' \item{fisher }{Fisher's alpha; as implemented in the \pkg{vegan} package (Fisher et al. (1943)).}
#'
#' \item{coverage }{Number of species needed to cover 50\% of the ecosystem.
#' For other quantiles, apply the function coverage directly.}
#' }
#'
#' @references
#'
#' Beisel J-N. et al. (2003)
#' A Comparative Analysis of Diversity Index Sensitivity.
#' _Internal Rev. Hydrobiol._ 88(1):3-15.
#' \url{https://portais.ufg.br/up/202/o/2003-comparative_evennes_index.pdf}
#'
#' Bulla L. (1994)
#' An  index of diversity and its associated diversity measure.
#' _Oikos_ 70:167--171
#'
#' Fisher, R.A., Corbet, A.S. & Williams, C.B. (1943).
#' The relation between the number of species and the number of individuals in a
#' random sample of animal population.
#' _Journal of Animal Ecology_ *12*, 42-58.
#'
#' Magurran AE, McGill BJ, eds (2011)
#' Biological Diversity: Frontiers in Measurement and Assessment
#' (Oxford Univ Press, Oxford), Vol 12.
#'
#' Smith B and Wilson JB. (1996)
#' A Consumer's Guide to Diversity Indices.
#' _Oikos_ 76(1):70-82.
#'
#' @seealso
#' \code{\link[scater:plotColData]{plotColData}}
#' \itemize{
#'   \item{\code{\link[mia:estimateRichness]{estimateRichness}}}
#'   \item{\code{\link[mia:estimateEvenness]{estimateEvenness}}}
#'   \item{\code{\link[mia:estimateDominance]{estimateDominance}}}
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
#' se <- estimateSimpsonDiversity(se)
#' colData(se)$gini_simpson
#'
#' # Calculates all the diversity indices
#' se <- estimateDiversity(se)
#' # All the indices' names
#' indices <- c("shannon","gini_simpson","inverse_simpson", "coverage", "fisher")
#' colData(se)[,indices]
#'
#' # plotting the diversities
#' library(scater)
#' plotColData(se, "shannon")
#' # ... by sample type
#' plotColData(se, "shannon", "SampleType")
#' \donttest{
#' # combining different plots
#' plots <- lapply(c("shannon","gini_simpson"),
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
                    index = c("shannon","gini_simpson", "inverse_simpson",
                              "coverage", "fisher"),
                    name = index, ...)
               standardGeneric("estimateDiversity"))

#' @rdname estimateDiversity
#' @export
setGeneric("estimateShannon",signature = c("x"),
           function(x, ...)
               standardGeneric("estimateShannon"))

#' @rdname estimateDiversity
#' @export
setGeneric("estimateSimpsonDiversity",signature = c("x"),
           function(x, ...)
               standardGeneric("estimateSimpsonDiversity"))

#' @rdname estimateDiversity
#' @export
setGeneric("estimateInvSimpson",signature = c("x"),
           function(x, ...)
               standardGeneric("estimateInvSimpson"))

#' @rdname estimateDiversity
#' @export
setGeneric("estimateCoverage",signature = c("x"),
           function(x, ...)
               standardGeneric("estimateCoverage"))

#' @rdname estimateDiversity
#' @export
setGeneric("estimateFisher",signature = c("x"),
           function(x, ...)
               standardGeneric("estimateFisher"))

#' @rdname estimateDiversity
#' @export
setMethod("estimateDiversity", signature = c(x = "SummarizedExperiment"),
    function(x, abund_values = "counts",
             index = c("shannon","gini_simpson","inverse_simpson",
                       "coverage", "fisher"),
             name = index, ..., BPPARAM = SerialParam()){

        # input check
        index<- match.arg(index, several.ok = TRUE)
        if(!.is_non_empty_character(name) || length(name) != length(index)){
            stop("'name' must be a non-empty character value and have the ",
                 "same length than 'index'.",
                 call. = FALSE)
        }
        .check_abund_values(abund_values, x)
        .require_package("vegan")

        dvrsts <- BiocParallel::bplapply(index,
                                         .get_diversity_values,
                                         x = x,
                                         mat = assay(x, abund_values),
                                         BPPARAM = BPPARAM,
                                         ...)
        .add_values_to_colData(x, dvrsts, name)
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
setMethod("estimateSimpsonDiversity", signature = c(x = "SummarizedExperiment"),
    function(x, ...){
        estimateDiversity(x, index = "gini_simpson", ...)
    }
)

#' @rdname estimateDiversity
#' @export
setMethod("estimateInvSimpson", signature = c(x = "SummarizedExperiment"),
    function(x, ...){
        estimateDiversity(x, index = "inverse_simpson",  ...)
    }
)


#' @rdname estimateDiversity
#' @export
setMethod("estimateCoverage", signature = c(x = "SummarizedExperiment"),
    function(x, ...){
        estimateDiversity(x, index = "coverage", ...)
    }
)


#' @rdname estimateDiversity
#' @export
setMethod("estimateFisher", signature = c(x = "SummarizedExperiment"),
    function(x, ...){
        estimateFisher(x, index = "fisher", ...)
    }
)

.calc_shannon <- function(mat, ...){
    vegan::diversity(t(mat), index="shannon")
}

# NOTE: vegan::diversity(x, index = "simpson")
# gives Simpson diversity, also called Gini-Simpson
# index: 1-lambda, where lambda is the Simpson index
# (lambda). This may cause confusion if your familiarity
# with diversity indices is limited.
# Moreover, Simpson's lambda is simply the
# squared sum of relative abundances so we can
# just use that for clarity and simplicity.
#.get_simpson <- function(x, ...){
.simpson_lambda <- function(mat, ...){

    # The function gives 1-lambda
    # Therefore we must take
    # lambda = 1 - vegan::diversity(t(x), index="simpson")
    # sum((x/sum(x))^2)

    # Convert table to relative values
    rel <- .calc_rel_abund(mat)

    # Squared sum of relative abundances
    colSums2(rel^2)
}

.calc_gini_simpson <- function(mat, ...){
    1 - .simpson_lambda(mat, ...)
}

.calc_inverse_simpson <- function(mat, ...){
    1 / .simpson_lambda(mat, ...)
}

.calc_coverage <- function(mat, threshold = 0.9, ...){

    # Threshold must be a numeric value between 0-1
    if( !( is.numeric(threshold) && (threshold >= 0 || threshold <= 1) ) ){
        stop("'threshold' must be a numeric value between 0-1.",
             call. = FALSE)
    }

    # Convert table to relative values
    rel <- .calc_rel_abund(mat)

    # Number of groups needed to have threshold (e.g. 50 %) of the ecosystem occupied
    coverage <- apply(rel, 2, function(x) {
        min(which(cumsum(rev(sort(x/sum(x)))) >= threshold))
    })
    names(coverage) <- colnames(rel)
    coverage
}

.calc_fisher <- function(mat, ...){
    vegan::fisher.alpha(t(mat))
}

#' @importFrom SummarizedExperiment assay assays
.get_diversity_values <- function(x, i, mat, ...){
    dvrsty_FUN <- switch(i,
                        shannon = .calc_shannon,
                        gini_simpson = .calc_gini_simpson,
                        inverse_simpson = .calc_inverse_simpson,
                        coverage = .calc_coverage,
                        fisher = .calc_fisher
                        )

    dvrsty <- dvrsty_FUN(mat, ...)
    dvrsty
}
