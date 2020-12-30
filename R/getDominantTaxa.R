#' Get Dominant Taxa
#'
#' This function return the most dominant taxa.
#'
#' @param x
#' A \code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#' object
#'
#' @param rank A single character defining a taxonomic rank. Must be a value of
#'   \code{taxonomicRanks()} function.
#'
#' @details
#' \code{getDominantTaxa} extracts the most abundant \dQuote{FeatureID}s
#' in a \code{\link[=SummarizedExperiment-class]{SummarizedExperiment}} object.
#'
#' @return \code{x} with additional \code{\link{colData}} named
#'   \code{*name*}
#'
#' @name getDominantTaxa
#' @export
#'
#' @author Leo Lahti and Tuomas Borman. Contact: \url{microbiome.github.io}
#'
#' @examples
#' data(GlobalPatterns)
#'
#' #Finds the dominant taxa. Taxa are returned as a vector.
#' getDominantTaxa(GlobalPatterns)
#'
#' #If taxonomic information is available, it is possible to find the most dominant
#' #group from specific taxonomic level, here family level.
#' getDominantTaxa(GlobalPatterns, rank="Family")
NULL

TAXONOMY_RANKS <- c("domain","kingdom","phylum","class","order","family",
                    "genus","species")

#' @rdname getDominantTaxa
#' @export
setGeneric("getDominantTaxa",signature = c("x"),
           function(x, rank = taxonomyRanks(x)[1], name = "Dominant Taxa")
               standardGeneric("getDominantTaxa"))


#' @rdname getDominantTaxa
#' @export
setMethod("getDominantTaxa", signature = c(x = "SummarizedExperiment"),
          function(x, rank = taxonomyRanks(x)[1], name = "Dominant Taxa"){

              #Input check
              if(!.is_non_empty_string(rank)){
                  stop("'rank' must be an non empty single character value.",
                       call. = FALSE)
              }

              #If "rank" is not NULL, species are aggregated according to the taxonomic rank
              #that is specified by user.
              if (!is.null(rank)) {
                  x <- agglomerateByRank(x, rank = rank)
              }

              #assays()$counts returns counts of the taxa in the samples
              #apply() function finds the indices of taxa's that has the highest
              #amount of counts.
              #names() returns the names of taxa that are the most abundant.
              taxas <- names(x)[apply(assays(x)$counts, 2, which.max)]

              .add_dominant_taxas_to_colData(x, taxas, name)
          }

)

#--------------------------Help functions-----------------------------------------
#' @importFrom SummarizedExperiment colData colData<-
#' @importFrom S4Vectors DataFrame
.add_dominant_taxas_to_colData <- function(x, dominances, name){
    dominances <- DataFrame(dominances)
    colnames(dominances) <- name
    colData(x)[,name] <- dominances
    x
}
