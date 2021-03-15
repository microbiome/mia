#' Estimate richness measures
#'
#' Several functions for calculation of community richness indices available via
#' wrapper functions. They are implemented via the \code{vegan} package.
#'
#' These include the \sQuote{Observed}, \sQuote{Chao1} and
#' \sQuote{ACE} and \sQuote{Hill} richness measures.
#' See details for more information and references.
#'
#' @param x a \code{\link{SummarizedExperiment}} object
#'
#' @param abund_values the name of the assay used for calculation of the
#'   sample-wise estimates
#'
#' @param index a \code{character} vector, specifying the richness measures
#'   to be calculated.
#'
#' @param name a name for the column(s) of the colData the results should be
#'   stored in.
#'
#' @param detection a numeric value for selecting detection threshold
#' for the abundances. The default detection threshold is 0.
#'
#' @param BPPARAM A
#'   \code{\link[BiocParallel:BiocParallelParam-class]{BiocParallelParam}}
#'   object specifying whether calculation of estimates should be parallelized.
#'
#' @param ... additional parameters passed to \code{estimateRichness}
#'
#' @return \code{x} with additional \code{\link{colData}} named
#'   \code{*name*}
#'
#' @details
#'
#' The richness is calculated per sample. This is a standard index in community
#' ecology, and it provides an estimate of the number of unique species in the
#' community. This is often not directly observed for the whole community but only
#' for a limited sample from the community. This has led to alternative richness
#* indices that provide different ways to estimate the species richness.
#'
#' Richness index differs from the concept of species diversity or evenness in
#' that it ignores species abundance, and focuses on the binary presence/absence
#' values that indicate simply whether the species was detected.
#'
#' The following richness indices are provided (not case-sensitive):
#'
#' \itemize{
#'   \item{Observed}{ The _observed richness_ gives the number of species that
#'   is detected above a given \code{detection} threshold in the observed sample
#'   (default 0). This is conceptually the simplest richness index. The corresponding
#'   index in the \pkg{vegan} package is "richness".}
#'
#'   \item{Chao1}{ This is a nonparametric estimator of species richness. It
#'   assumes that rare species carry information about the (unknown) number
#'   of unobserved species. We use here the bias-corrected version
#'   (O'Hara 2005, Chiu et al. 2014) implemented in
#'   \code{\link[vegan:specpool]{estimateR}}. This index implicitly
#'   assumes that every taxa has equal probability of being observed. Note
#'   that it gives a lower bound to species richness. The bias-corrected
#'   For an exact formulation, see \code{\link[vegan:specpool]{estimateR}}.
#'   This estimator uses only the singleton and doubleton counts, and
#'   hence it gives more weight to the low abundance species.
#'   Note that this index comes with an additional column with standard
#'   error information.}
#'
#'   \item{ACE}{ Abundance-based coverage estimator (ACE) is another nonparametric richness
#'   index that uses sample coverage, defined based on the sum of the probabilities
#'   of the observed species. This method divides the species into abundant (more than 10
#'   reads or observations) and rare groups
#'   in a sample and tends to underestimate the real number of species. The ACE index
#'   ignores the abundance information for the abundant species,
#'   based on the assumption that the abundant species are observed regardless of their
#'   exact abundance. We use here the bias-corrected version
#'   (O'Hara 2005, Chiu et al. 2014) implemented in
#'   \code{\link[vegan:specpool]{estimateR}}.
#'   For an exact formulation, see \code{\link[vegan:specpool]{estimateR}}.
#'   Note that this index comes with an additional column with standard
#'   error information.}
#'
#'   \item{Hill}{ Effective species richness aka Hill index (see e.g. Chao et al. 2016).
#'   Currently only the case ${}^1D$ is implemented. This corresponds to the exponent
#'   of Shannon diversity. Intuitively, the effective richness indicates the number of
#'   species whose even distribution would lead to tha same diversity than the observed
#'   community, where the species abundances are unevenly distributed.}
#' }
#'
#'
#' @references
#'
#' Chao A. (1984)
#' Non-parametric estimation of the number of classes in a population.
#' _Scand J Stat._ 11:265–270.
#'
#' Chao A, Chun-Huo C, Jost L (2016).
#' Phylogenetic Diversity Measures and Their Decomposition:
#' A Framework Based on Hill Numbers. Biodiversity Conservation and Phylogenetic Systematics,
#' Springer International Publishing, pp. 141–172, doi:10.1007/978-3-319-22461-9_8.
#'
#' Chiu, C.H., Wang, Y.T., Walther, B.A. & Chao, A. (2014).
#' Improved nonparametric lower bound of species richness via a modified
#' Good-Turing frequency formula.
#' _Biometrics_ 70, 671-682.
#'
#' O'Hara, R.B. (2005).
#' Species richness estimators: how many species can dance on the head of a pin?
#' _J. Anim. Ecol._ 74, 375-386.
#'
#' @seealso
#' \code{\link[scater:plotColData]{plotColData}}
#' \itemize{
#'   \item{\code{\link[vegan:specpool]{estimateR}}}
#' }
#'
#' @name estimateRichness
#'
#' @export
#'
#' @author Leo Lahti. Contact: \url{microbiome.github.io}
#'
#' @examples
#' data(esophagus)
#'
#' # Calculates observed richness index
#' esophagus <- estimateRichness(esophagus, index="observed")
#' # Shows all indices
#' colData(esophagus)
#'
#' # Calculate richness excluding singletons (detection limit 1)
#' esophagus <- estimateRichness(esophagus, index="observed", detection = 1)
#'
#' # Indices must be written correctly (e.g. ACE, not ace), otherwise an error
#' # gets thrown
#' \dontrun{esophagus <- estimateRichness(esophagus, index="ace")}
#'
#' # Calculates Chao1 and ACE indices
#' esophagus <- estimateRichness(esophagus, index=c("Chao1", "ACE"))
#' # Shows all indices
#' colData(esophagus)
#' # Shows ACE index
#' colData(esophagus)$ACE
#'
#' # Deletes ACE index
#' colData(esophagus)$ACE <- NULL
#' # Shows all indices, ACE is deleted
#' colData(esophagus)
#' # Deletes all indices
#' colData(esophagus) <- NULL
#'
#' # Names of columns can be chosen, but the length of arguments must match.
#' esophagus <- estimateRichness(esophagus,
#'                                index = c("ACE", "Chao1"),
#'                                name = c("index1", "index2"))
#' # Shows all indices
#' colData(esophagus)
#'
#' # Deletes all indices
#' colData(esophagus) <- NULL
#'
#' # Calculates all indices
#' esophagus <- estimateRichness(esophagus)
#' # Shows all indices
#' colData(esophagus)
NULL

#' @rdname estimateRichness
#' @export
setGeneric("estimateRichness",signature = c("x"),
        function(x, abund_values = "counts",
                 index = c("observed", "Chao1", "ACE", "Hill"),
                 name = index,
                 detection = 0,
                 ...,
                 BPPARAM = SerialParam())
            standardGeneric("estimateRichness"))

#' @rdname estimateRichness
#' @export
setMethod("estimateRichness", signature = c(x = "SummarizedExperiment"),
    function(x,
            abund_values = "counts",
            index = c("observed", "Chao1", "ACE", "Hill"),
            name = index,
            detection = 0,
            ...,
            BPPARAM = SerialParam()){

        # Input check
        # Check abund_values
        .check_abund_values(abund_values, x)
        # Check indices
        index <- match.arg(index, several.ok = TRUE)
        if(!.is_non_empty_character(name) || length(name) != length(index)){
            stop("'name' must be a non-empty character value and have the ",
                "same length than 'index'.",
                call. = FALSE)
        }
        # Calculates richness indices
        richness <- BiocParallel::bplapply(index,
                                            FUN = .get_richness_values,
                                            assay = assay(x, abund_values),
                                            detection = detection,
                                            BPPARAM = BPPARAM)
        # Add richness indices to colData
        .add_values_to_colData(x, richness, name)
    }
)


.calc_observed <- function(mat, detection, ...){
    # vegan::estimateR(t(mat))["S.obs",]
    colSums(mat > detection)
}

.calc_chao1 <- function(mat, ...){
    # Required to work with DelayedArray
    if(is(mat, "DelayedArray")) {
      mat <- matrix(mat, nrow = nrow(mat))
    }

    ans <- t(vegan::estimateR(t(mat))[c("S.chao1","se.chao1"),])
    colnames(ans) <- c("","se")
    ans
}

.calc_ACE <- function(mat, ...){
    # Required to work with DelayedArray
    if(is(mat, "DelayedArray")) {
      mat <- matrix(mat, nrow = nrow(mat))
    }

    ans <- t(vegan::estimateR(t(mat))[c("S.ACE","se.ACE"),])
    colnames(ans) <- c("","se")
    ans
}

.calc_hill <- function(mat, ...){
    # Exponent of Shannon diversity
    exp(vegan::diversity(t(mat), index="shannon"))
}

.get_richness_values <- function(index, assay, detection) {

    FUN <- switch(index,
                observed = .calc_observed,
                Chao1 = .calc_chao1,
                ACE = .calc_ACE,
                Hill = .calc_hill
        )
    do.call(FUN,
            list(mat = assay,
            detection = detection
        ))
}
