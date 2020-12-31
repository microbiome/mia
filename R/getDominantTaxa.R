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
#' x <- GlobalPatterns
#'
#' #Finds the dominant taxa.
#' x <- getDominantTaxa(GlobalPatterns)
#' #Information is stored to colData
#' colData(x)
#'
#' #If taxonomic information is available, it is possible to find the most dominant
#' #group from specific taxonomic level, here family level. The name of column can be specified.
#' x <- getDominantTaxa(GlobalPatterns, rank="Family", name="Dominant Family")
#' colData(x)
#'
NULL

#' @rdname getDominantTaxa
#' @export
setGeneric("getDominantTaxa",signature = c("x"),
           function(x, rank = NULL, name = "Dominant Taxa")
               standardGeneric("getDominantTaxa"))


#' @rdname getDominantTaxa
#' @export
setMethod("getDominantTaxa", signature = c(x = "SummarizedExperiment"),
          function(x, rank = NULL, name = "Dominant Taxa"){

              #Input check
              if(!is.null(rank)){
                  if(!.is_a_string(rank)){
                      stop("'rank' must be an single character value.",
                           call. = FALSE)
                  }
                  TAXONOMY_RANKS <- c("domain","kingdom","phylum","class","order","family",
                                      "genus","species")
                  if(!rank %in% TAXONOMY_RANKS){
                      stop("'rank' must be a taxonomy rank.",
                           call. = FALSE)
                  }

              }

              if(!.is_non_empty_string(name)){
                  stop("'name' must be a non-empty single character value.",
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
