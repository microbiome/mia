#' Estimate Dominance
#'
#' This function calculates the community dominance index.
#'
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
#' \code{estimateDominance} calculates the community dominance indices.
#' \itemize{
#' \item{'DBP' }{Berger-Parker index is calculated similarly than relative index.
#' DBP is the relative abundance of the most abundant species of the sample.
#' ($(N_max/N_tot)$, where N_max is the absolute abundance of the most dominant species and
#' N_tot is the sum of absolute abundances of all species.) Index gives values in interval [0,1],
#' where bigger value represent greater dominance. (Maiti 2012 p. 29.)}
#' \item{'DMN' }{McNaughton’s index is the sum of relative abundances of the two most
#' abundant species of the sample. ($(N_{2 maxs}/N_tot)$, where N_{2 maxs} is the absolute
#' abundance of the two most dominant species and N_tot is the sum of absolute
#' abundances of all species.) Index gives values in interval [0,1],
#' where bigger value represent greater dominance.}
#' \item{'absolute' }{Absolute index equals to the absolute abundance of the most
#' dominant species of the sample. Index gives positive integer values.}
#' \item{'relative' }{Relative index equals to the relative abundance of the most
#' dominant species of the sample. ($(N_max/N_tot)$, where N_max is the absolute abundance of
#' the most dominant species and N_tot is the sum of absolute abundances of all species.)
#' Index gives values in interval [0,1], where bigger value represent greater dominance.}
#' \item{'simpson' }{Simpson's index, or Simpson's dominance index, is calculated
#' by raising all relative abundances of species to the power of 2, and then summing them together.
#' (($sum(p^2)$), where p is relative abundances) Index gives values in interval [0,1].
#' Value equals the probability that two randomly chosen individuals belongs to the same species.
#' The higher the probability, the greater the dominance is. (Thukral et al. 2019.)}
#' \item{'core_abundance' }{Core abundance index is related to core species.
#' Core species are species that are most abundant in all samples, i.e., in whole dataset.
#' Core species are defined as those species that have prevalence over 50\%.
#' It means that in order to belong to core species, species must be prevalent in 50 % of samples.
#' Core species are used to calculate the core abundance index.
#' Core abundance index is sum of relative abundances of core species in the sample.
#' ($(N_core/N_tot)$, where N_core is the sum of absolute abundance of the core species and
#' N_tot is the sum of absolute abundances of all species.) Index gives values in interval [0,1],
#' where bigger value represent greater dominance.}
#' \item{'gini' }{Gini index is probably best-known from economical relations.
#' In economics, it is used to measure, for example, how unevenly wealth
#' is distributed among population. Here, Gini index is used similarly, but wealth
#' is replaced with abundance. If there are small group of species that represent
#' large portion of total abundance of microbes in the sample, the inequality is large and
#' Gini index closer to 1. If all species has equally large abundances, the equality
#' is perfect and Gini index equals 0. (Sitthiyot & Holasut 2020.)}
#' }
#'
#' @references
#' Maiti SK (2012) Ecorestoration of the coalmine degraded lands.
#'
#' Sitthiyot T & Holasut K (2020) A simple method for measuring inequality.
#' Palgrave Commun 6(112) doi: 10.1057/s41599-020-0484-6.
#'
#' Thukral AK, Bhardwaj R, Kumar V & Sharma A (2019)
#' New indices regarding the dominance and diversity of communities,
#' derived from sample variance and standard deviation.
#' Heliyon 5(10). doi: 10.1016/j.heliyon.2019.e02606.
#'
#'
#' @return \code{x} with additional \code{\link{colData}} named
#'   \code{*name*}
#'
#' @seealso
#' \itemize{
#'   \item{\code{\link[=estimateDiversity]{estimateDiversity}}}
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
#' #Calculates simpson dominance index
#' esophagus <- estimateDominance(esophagus, index="simpson")
#' #Shows all indices
#' colData(esophagus)
#'
#' #Indices must be written correctly (e.g. DBP, not dbp), otherwise an error
#' # gets thrown
#' \dontrun{esophagus <- estimateDominance(esophagus, index="dbp")}
#' #Calculates DBP and Core Abundance indices
#' esophagus <- estimateDominance(esophagus, index=c("DBP", "core_abundance"))
#' #Shows all indices
#' colData(esophagus)
#' #Shows DBP index
#' colData(esophagus)$DBP
#'
#' #Deletes DBP index
#' colData(esophagus)$DBP <- NULL
#' #Shows all indices, DBP is deleted
#' colData(esophagus)
#' #Deletes all indices
#' colData(esophagus) <- NULL
#'
#' #Names of columns can be chosen, but the length of arguments must match.
#' esophagus <- estimateDominance(esophagus,
#'                                index = c("DBP", "core_abundance"),
#'                                name = c("index1", "index2"))
#' #Shows all indices
#' colData(esophagus)
#' #If they do not match, gets an error.
#' \dontrun{
#' esophagus <- estimateDominance(esophagus,
#'                                index="simpson",
#'                                name = c("index3", "index4"))
#' }
#' #Shows all indices
#' colData(esophagus)
#' #Deletes all indices
#' colData(esophagus) <- NULL
#'
#' #Calculates all indices
#' esophagus <- estimateDominance(esophagus)
#' #Shows all indices
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
        #Input check
        #Check abund_values
        .check_abund_values(abund_values, x)
        #Check indices
        index <- match.arg(index, several.ok = TRUE)
        if(!.is_non_empty_character(name) || length(name) != length(index)){
            stop("'name' must be a non-empty character value and have the ",
                 "same length than 'index'.",
                 call. = FALSE)
        }
        #Check ntaxa
        if(!(ntaxa>0 && ntaxa<3)){
            stop("'ntaxa' must be a numerical value 1 or 2.", call. = FALSE)
        }
        #Check aggregate
        if(!.is_a_bool(aggregate)){
            stop("'aggregate' must be TRUE or FALSE.", call. = FALSE)
        }
        #
        #Calculates dominance indices
        dominances <- BiocParallel::bplapply(index,
                                             FUN = .get_dominances_values,
                                             assay = assay(x,abund_values),
                                             ntaxa = ntaxa,
                                             aggregate = aggregate,
                                             BPPARAM = BPPARAM)
        #Add dominance indices to colData
        .add_dominances_values_to_colData(x, dominances, name)
    }
)




#---------------------------Help functions----------------------------------------------------------------

# x: Species count vector
.simpson_dominance <- function(x, zeroes=TRUE) {

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

.get_simpson_dominance <- function(x, ...){
    apply(x, 2L, .simpson_dominance)
}

.gini_dominance <- function(x, w=rep(1, length(x))) {
    # See also reldist::gini for an independent implementation
    o <- order(x)
    x <- x[o]
    w <- w[o]/sum(w)
    p <- cumsum(w)
    nu <- cumsum(w * x)
    n <- length(nu)
    nu <- nu/nu[[n]]
    sum(nu[-1] * p[-n]) - sum(nu[-n] * p[-1])
}

.get_gini_dominance <- function(x, ...){
    apply(x, 2L, .gini_dominance)
}

.get_core_dominance <- function(x, ...){
    getPrevalentAbundance(x, detection = 0, as_relative = TRUE)
}

.get_dominance <- function(x, ntaxa, aggregate, index){
    if (index == "absolute") {
        #ntaxa=1 by default but can be tuned
        as_relative <- FALSE
    } else if (index == "relative") {
        #ntaxa=1 by default but can be tuned
        as_relative <- TRUE
    } else if (index == "DBP") {
        #Berger-Parker: if selected fix the following values
        ntaxa <- 1
        as_relative <- TRUE
    } else if (index == "DMN") {
        #McNaughton's dominance: if selected fix the following values
        ntaxa <- 2
        aggregate <- TRUE
        as_relative <- TRUE
    }

    if (as_relative) {
        #Calculates the relative abundance per sample
        x <- apply(x, 2L,
                   function(x) {
                       x/sum(x, na.rm=TRUE)
                   })
    }

    #Aggregate or not
    if (!aggregate) {
        ans <- apply(x, 2L,
                     function(x) {
                         sort(x, decreasing = TRUE)[[ntaxa]]
                     })
    } else {
        ans <- apply(x, 2L,
                     function(x) {
                         sum(sort(x, decreasing = TRUE)[seq_len(ntaxa)])
                     })
    }

    #add sample names to the table
    names(ans) <- colnames(x)
    ans
}

.get_dominances_values <- function(index, assay, ntaxa = 1, aggregate = TRUE) {
    FUN <- switch(index,
                  simpson = .get_simpson_dominance,
                  core_abundance = .get_core_dominance,
                  gini = .get_gini_dominance,
                  .get_dominance)
    do.call(FUN,
            list(x = assay,
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
