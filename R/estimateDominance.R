#' Estimate Dominance
#'
#' This function calculates the community dominance index.
#' This includes the \sQuote{DBP}, \sQuote{DMN}, \sQuote{Absolute},
#' \sQuote{Relative}, \sQuote{Core Abundance}, \sQuote{Gini}
#'
#' @param x a
#'   \code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#'   object
#'
#' @param abund_values A single character value for selecting the
#'   \code{\link[SummarizedExperiment:SummarizedExperiment-class]{assay}}
#'   to use for prevalence calculation.
#'
#' @param index Specifies the indices which are calculated.
#'
#' @param ntaxa Optional. The rank of the dominant taxa to consider. Disregarded
#'   for the \code{index} \dQuote{gini}, \dQuote{simpson},
#'   \dQuote{core_abundance}, \dQuote{DBP} and \dQuote{DMN}.
#'
#' @param aggregate (Optional, default = TRUE) Aggregate the top members or not.
#'   If aggregate=TRUE, then the sum of relative abundances is returned.
#'   Otherwise the relative abundance is returned for the single taxa with
#'   the indicated rank. Disregarded
#'   for the \code{index} \dQuote{gini}, \dQuote{simpson},
#'   \dQuote{core_abundance}, \dQuote{DMN}.
#'
#' @param name A name for the column of the colData where the calculated
#' Dominance indices should be stored in.
#'
#' @param ... additional arguments
#'
#' @param BPPARAM A
#'   \code{\link[BiocParallel:BiocParallelParam-class]{BiocParallelParam}}
#'   object specifying whether calculation of estimates should be parallelized.
#'   (Currently not used)
#'
#' @details
#' \code{estimateDominance} calculates the following community dominance indices.
#'
#' \itemize{
#' \item{'DBP'}{ Berger-Parker index is calculated similarly than relative index.
#' DBP is the relative abundance of the most abundant species of the sample.
#' Index gives values in interval 0 to 1, where bigger value represent greater dominance.
#' (See e.g. Berger & Parker 1970.)
#'
#' \deqn{DBP = \frac{N_1}{N_{tot}}}{%
#' DBP = N_1/N_tot}
#' where $N_1$ is the absolute abundance of the most dominant species and
#' $N_{tot}$ is the sum of absolute abundances of all species.}
#'
#' \item{'DMN'}{ McNaughton’s index is the sum of relative abundances of the two most
#' abundant species of the sample. Index gives values in interval 0 to 1,
#' where bigger value represent greater dominance.
#'
#' \deqn{DMN = \frac{N_1 + N_2}{N_{tot}}}{%
#' DMN = (N_1 + N_2)/N_tot}
#' where $N_1$ and $N_2$ are the absolute
#' abundances of the two most dominant species and $N_{tot}$ is the sum of absolute
#' abundances of all species.}
#'
#' \item{'absolute'}{ Absolute index equals to the absolute abundance of the most
#' dominant species of the sample. Index gives positive integer values.}
#'
#' \item{'relative'}{ Relative index equals to the relative abundance of the most
#' dominant species of the sample. Index gives values in interval 0 to 1, where bigger
#' value represent greater dominance.
#'
#' \deqn{relative = \frac{N_1}{N_{tot}}}{%
#' relative = N_1/N_tot}
#' where $N_1$ is the absolute abundance of
#' the most dominant species and $N_{tot}$ is the sum of absolute abundances of all species.}
#'
#' \item{'simpson'}{ Simpson's index, or Simpson's dominance index, is calculated
#' by raising all relative abundances of species to the power of 2, and then summing them together.
#' Index gives values in interval 0 to 1. Value equals the probability that two randomly
#' chosen individuals belongs to the same species. The higher the probability, the greater
#' the dominance is. (See e.g. Simpson 1949.)
#'
#' \deqn{simpson = \sum(p^2)}{%
#' simpson = \sum(p^2)}
#' where $p$ is relative abundances.}
#'
#' \item{'core_abundance'}{ Core abundance index is related to core species.
#' Core species are species that are most abundant in all samples, i.e., in whole data set.
#  Core species are defined as those species that have prevalence over 50\%.
#  It means that in order to belong to core species, species must be prevalent in 50\% of samples.
#' Core species are used to calculate the core abundance index.
#' Core abundance index is sum of relative abundances of core species in the sample.
#' Index gives values in interval 0 to 1, where bigger value represent greater dominance.
#'
#' \deqn{core_abundance = \frac{N_{core}}{N_{tot}}}{%
#' core_abundance = N_core/N_tot}
#' where $N_{core}$ is the sum of absolute abundance of the core species and
#' $N_{tot}$ is the sum of absolute abundances of all species.}
#'
#' \item{'gini'}{ Gini index is probably best-known from socio-economic contexts.
#' In economics, it is used to measure, for example, how unevenly income
#' is distributed among population. Here, Gini index is used similarly, but income
#' is replaced with abundance. If there is small group of species that represent
#' large portion of total abundance of microbes, the inequality is large and
#' Gini index closer to 1. If all species has equally large abundances, the equality
#' is perfect and Gini index equals 0. (See e.g. Gini 1921.)}
#' }
#'
#' @references
#' Berger WH & Parker FL (1970) Diversity of Planktonic Foraminifera in Deep-Sea Sediments.
#' Science 168(3937): 1345-1347. doi: 10.1126/science.168.3937.1345
#'
#' Gini C (1921) Measurement of Inequality of Incomes.
#' The Economic Journal 31(121): 124-126. doi: 10.2307/2223319
#'
#' Simpson EH (1949) Measurement of Diversity.
#' Nature 163(688). doi: 10.1038/163688a0
#'
#' @return \code{x} with additional \code{\link{colData}} named
#'   \code{*name*}
#'
#' @seealso
#' \itemize{
#'   \item{\code{\link[mia:estimateDiversity]{estimateDiversity}}}
#' }
#'
#' @name estimateDominance
#' @export
#'
#' @author Leo Lahti and Tuomas Borman. Contact: \url{microbiome.github.io}
#'
#' @examples
#' data(esophagus)
#'
#' # Calculates simpson dominance index
#' esophagus <- estimateDominance(esophagus, index="simpson")
#' # Shows all indices
#' colData(esophagus)
#'
#' # Indices must be written correctly (e.g. DBP, not dbp), otherwise an error
#' # gets thrown
#' \dontrun{esophagus <- estimateDominance(esophagus, index="dbp")}
#' # Calculates DBP and Core Abundance indices
#' esophagus <- estimateDominance(esophagus, index=c("DBP", "core_abundance"))
#' # Shows all indices
#' colData(esophagus)
#' # Shows DBP index
#' colData(esophagus)$DBP
#'
#' # Deletes DBP index
#' colData(esophagus)$DBP <- NULL
#' # Shows all indices, DBP is deleted
#' colData(esophagus)
#' # Deletes all indices
#' colData(esophagus) <- NULL
#'
#' # Names of columns can be chosen, but the length of arguments must match.
#' esophagus <- estimateDominance(esophagus,
#'                                index = c("DBP", "core_abundance"),
#'                                name = c("index1", "index2"))
#' # Shows all indices
#' colData(esophagus)
#' # If they do not match, gets an error.
#' \dontrun{
#' esophagus <- estimateDominance(esophagus,
#'                                index="simpson",
#'                                name = c("index3", "index4"))
#' }
#' # Shows all indices
#' colData(esophagus)
#' # Deletes all indices
#' colData(esophagus) <- NULL
#'
#' # Calculates all indices
#' esophagus <- estimateDominance(esophagus)
#' # Shows all indices
#' colData(esophagus)
NULL

#' @rdname estimateDominance
#' @export
setGeneric("estimateDominance",signature = c("x"),
           function(x,
                    abund_values = "counts",
                    index = c("DBP", "DMN", "absolute", "relative", "simpson", "core_abundance", "gini"),
                    ntaxa = 1,
                    aggregate = TRUE,
                    name = index,
                    ...,
                    BPPARAM = SerialParam())
               standardGeneric("estimateDominance"))


#' @rdname estimateDominance
#' @export
setMethod("estimateDominance", signature = c(x = "SummarizedExperiment"),
    function(x,
             abund_values = "counts",
             index = c("DBP", "DMN", "absolute", "relative", "simpson", "core_abundance", "gini"),
             ntaxa = 1,
             aggregate = TRUE,
             name = index,
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
        # Check ntaxa
        if(!(ntaxa>0 && ntaxa<3)){
            stop("'ntaxa' must be a numerical value 1 or 2.", call. = FALSE)
        }
        # Check aggregate
        if(!.is_a_bool(aggregate)){
            stop("'aggregate' must be TRUE or FALSE.", call. = FALSE)
        }

        # Calculates dominance indices
        dominances <- BiocParallel::bplapply(index,
                                             FUN = .get_dominances_values,
                                             assay = assay(x,abund_values),
                                             ntaxa = ntaxa,
                                             aggregate = aggregate,
                                             BPPARAM = BPPARAM)
        # Add dominance indices to colData
        .add_dominances_values_to_colData(x, dominances, name)
    }
)




#---------------------------Help functions----------------------------------------------------------------

# x: Species count vector
.simpson_dominance <- function(x, zeroes=TRUE) {
    x <- as.vector(x)

    if (!zeroes) {
        x[x > 0]
    }

    # Relative abundances
    p <- x/sum(x)

    # Simpson index (has interpretation as dominance)
    lambda <- sum(p^2)

    # More advanced Simpson dominance (Simpson 1949) However let us not use
    # this as it is not in [0,1] and it is very highly correlated with the
    # simpler variant lambda Species richness (number of species)
    # S <- length(x) sum(p * (p - 1)) / (S * (S - 1))

    lambda

}

.calc_simpson_dominance <- function(mat, ...){
    apply(mat, 2L, .simpson_dominance)
}

.gini_dominance <- function(x, w=rep(1, length(x))) {
    # See also reldist::gini for an independent implementation
    x <- as.vector(x)
    o <- order(x)
    x <- x[o]
    w <- w[o]/sum(w)
    p <- cumsum(w)
    nu <- cumsum(w * x)
    n <- length(nu)
    nu <- nu/nu[[n]]
    sum(nu[-1] * p[-n]) - sum(nu[-n] * p[-1])
}

.calc_gini_dominance <- function(mat, ...){
    apply(mat, 2L, .gini_dominance)
}

.calc_core_dominance <- function(mat, ...){
    getPrevalentAbundance(mat, detection = 0, as_relative = TRUE)
}

.calc_dominance <- function(mat, ntaxa, aggregate, index){
    if (index == "absolute") {
        # ntaxa=1 by default but can be tuned
        as_relative <- FALSE
    } else if (index == "relative") {
        # ntaxa=1 by default but can be tuned
        as_relative <- TRUE
    } else if (index == "DBP") {
        # Berger-Parker: if selected fix the following values
        ntaxa <- 1
        as_relative <- TRUE
    } else if (index == "DMN") {
        # McNaughton's dominance: if selected fix the following values
        ntaxa <- 2
        aggregate <- TRUE
        as_relative <- TRUE
    }

    if (as_relative) {
        # Calculates the relative abundance per sample
        mat <- .calc_rel_abund(mat)
    }

    # Aggregate or not
    if (!aggregate) {
        idx <- apply(mat, 2L,
                     function(mc) {
                         order(as.vector(mc), decreasing = TRUE)[[ntaxa]]
                     })
    } else {
        idx <- apply(mat, 2L,
                     function(mc) {
                         order(as.vector(mc), decreasing = TRUE)[seq_len(ntaxa)]
                     })
        idx <- split(as.vector(idx),
                     unlist(lapply(seq_len(length(idx) / ntaxa),rep.int,ntaxa)))
    }

    ans <- lapply(mapply(function(i,j,x){x[i,j]},
                         i = idx,
                         j = seq_len(ncol(mat)),
                         MoreArgs = list(x = mat),
                         SIMPLIFY = FALSE),
                  sum)
    ans <- unlist(ans)

    # Adds sample names to the table
    names(ans) <- colnames(mat)
    ans
}

.get_dominances_values <- function(index, assay, ntaxa = 1, aggregate = TRUE) {
    FUN <- switch(index,
                  simpson = .calc_simpson_dominance,
                  core_abundance = .calc_core_dominance,
                  gini = .calc_gini_dominance,
                  .calc_dominance)
    do.call(FUN,
            list(mat = assay,
                 ntaxa = ntaxa,
                 aggregate = aggregate,
                 index = index))
}

#' @importFrom SummarizedExperiment colData colData<-
#' @importFrom S4Vectors DataFrame
.add_dominances_values_to_colData <- function(x, dominances, name){
    dominances <- DataFrame(dominances)
    colnames(dominances) <- name
    colData(x)[,name] <- dominances
    x
}
