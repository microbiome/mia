#' Estimate diversity measures
#'
#' Several functions for calculation of diversity indices available via
#' wrapper functions. Some of them are implemented via the \code{vegan} package.
#'
#' The available indices include the \sQuote{Shannon}, \sQuote{Gini-Simpson},
#' \sQuote{Inverse Simpson}, \sQuote{Coverage}, and \sQuote{Fisher alpha}
#' diversity indices. See details for more information and references.
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
#' @param ... optional arguments:
#' \itemize{
#'   \item{threshold}{ a numeric value in the unit interval,
#'   determining the threshold for coverage index. By default, the threshold is 0.9.}
#' }
#'
#' @return \code{x} with additional \code{\link{colData}} named \code{*name*}
#'
#' @details
#'
#' Diversity is a joint quantity that combines elements or community richness and evenness.
#' Diversity increases, in general, when species richness or evenness increase.
#'
#' By default, this function returns all indices.
#'
#' \itemize{
#' \item{'inverse_simpson' }{Inverse Simpson diversity:
#' \eqn{1/lambda} where \eqn{lambda=sum(p^2)} and p refers to relative abundances.
#' This corresponds to the diversity index
#' 'invsimpson' in vegan::diversity. Don't confuse this with the
#' closely related Gini-Simpson index}
#'
#' \item{'gini_simpson' }{Gini-Simpson diversity i.e. \eqn{1 - lambda}, where \eqn{lambda} is the
#' Simpson index, calculated as the sum of squared relative abundances.
#' This corresponds to the diversity index
#' 'simpson' in \code{\link[vegan:diversity]{vegan::diversity}}.
#' This is also called Gibbs–Martin, or Blau index in sociology,
#' psychology and management studies. The Gini-Simpson index (1-lambda) should not be
#' confused with Simpson's dominance (lambda), Gini index, or inverse Simpson index (1/lambda).}
#'
#' \item{'shannon' }{Shannon diversity (entropy).}
#'
#' \item{'fisher' }{Fisher's alpha; as implemented in
#' \code{\link[vegan:fisher.alpha]{vegan::fisher.alpha}}. (Fisher et al. (1943)).}
#'
#' \item{'coverage' }{Number of species needed to cover a given fraction of the ecosystem (50\% by default).
#' Tune this with the threshold argument.}
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
#' # All index names as known by the function
#' index <- c("shannon","gini_simpson","inverse_simpson", "coverage", "fisher")
#'
#' # Corresponding polished names
#' name <- c("Shannon","GiniSimpson","InverseSimpson", "Coverage", "Fisher")
#'
#' # Calculate diversities
#' se <- estimateDiversity(se, index = index)
#'
#' # The colData contains the indices with their code names by default
#' colData(se)[, index]
#'
#' # Removing indices
#' colData(se)[, index] <- NULL
#'
#' # It is recommended to specify also the final names used in the output.
#' se <- estimateDiversity(se,
#'   index = c("shannon", "gini_simpson", "inverse_simpson", "coverage", "fisher"),
#'    name = c("Shannon", "GiniSimpson",  "InverseSimpson",  "Coverage", "Fisher"))
#'
#' # The colData contains the indices by their new names provided by the user
#' colData(se)[, name]
#'
#' # Compare the indices visually
#' pairs(colData(se)[, name])
#'
#' # Plotting the diversities - use the selected names
#' library(scater)
#' plotColData(se, "Shannon")
#' # ... by sample type
#' plotColData(se, "Shannon", "SampleType")
#' \dontrun{
#' # combining different plots
#' library(patchwork)
#' plot_index <- c("Shannon","GiniSimpson")
#' plots <- lapply(plot_index,
#'                 plotColData,
#'                 object = se,
#'                 x = "SampleType",
#'                 colour_by = "SampleType")
#' plots <- lapply(plots,"+", theme(axis.text.x = element_text(angle=45,hjust=1)))
#' names(plots) <- plot_index
#' plots$Shannon + plots$GiniSimpson + plot_layout(guides = "collect")
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
        .check_assay_present(abund_values, x)
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
.get_diversity_values <- function(index, x, mat, ...){

    FUN <- switch(index,
                        shannon = .calc_shannon,
                        gini_simpson = .calc_gini_simpson,
                        inverse_simpson = .calc_inverse_simpson,
                        coverage = .calc_coverage,
                        fisher = .calc_fisher
                        )

    FUN(x = x, mat = mat, ...)

}
