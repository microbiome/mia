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
#' @param group With group, it is possible group observations in overview.
#'
#' @param name A name for the column of the colData where the dominant taxa
#' will be stored in.
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
#' x <- getDominantTaxa(x)
#' #Information is stored to colData
#' colData(x)
#'
#' #If taxonomic information is available, it is possible to find the most dominant
#' #group from specific taxonomic level, here family level. The name of column can be specified.
#' x <- getDominantTaxa(x, rank="Family", name="Dominant Family")
#' colData(x)
#'
#' x <- microbiomeDataSets::dietswap()
#' #With group, it is possbile to group observations based on specific group
#' x <- getDominantTaxa(x, group = "nationality")
#'
NULL

#' @rdname getDominantTaxa
#' @export
setGeneric("getDominantTaxa",signature = c("x"),
           function(x, rank = NULL, group = NULL, name = "dominant_taxa")
               standardGeneric("getDominantTaxa"))


#' @rdname getDominantTaxa
#' @export
setMethod("getDominantTaxa", signature = c(x = "SummarizedExperiment"),
          function(x, rank = NULL, group = NULL, name = "dominant_taxa"){

              #Input check
              #rank check
              if(!is.null(rank)){
                  if(!.is_a_string(rank)){
                      stop("'rank' must be an single character value.",
                           call. = FALSE)
                  }
                  .check_taxonomic_rank(rank, x)
              }

              #group check
              if(!is.null(group)){
                  if(isFALSE(any(group %in% colnames(colData(x))))){
                      stop("'group' variable must be in colnames(colData(x))")
                  }
              }

              #name check
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

              #Add taxas to colData
              x <- .add_dominant_taxas_to_colData(x, taxas, name)

              #Creates a tibble df that contains dominant taxa and number of times that they present in samples
              #and relative portion of samples where they present.
              if (is.null(group)) {
                  overview <- S4Vectors::as.data.frame(colData(y)) %>%
                      dplyr::group_by(dominant_taxa) %>%
                      dplyr::tally() %>%
                      mutate(
                          rel.freq = round(100 * n / sum(n), 1),
                          rel.freq.pct = paste0(round(100 * n / sum(n), 0), "%")
                      ) %>%
                      dplyr::arrange(desc(n))
              } else {
                  group <- sym(group)
                  # dominant_taxa <-"dominant_taxa"
                  overview <- S4Vectors::as.data.frame(colData(y)) %>%
                      dplyr::group_by(!!group, dominant_taxa) %>%
                      dplyr::tally() %>%
                      mutate(
                          rel.freq = round(100 * n / sum(n), 1),
                          rel.freq.pct = paste0(round(100 * n / sum(n), 0), "%")
                      ) %>%
                      dplyr::arrange(desc(n))
              }

              #Prints overview
              print(overview)

              return(x)

          }

)

################################HELP FUNCTIONS#####################################
#' @importFrom SummarizedExperiment colData colData<-
#' @importFrom S4Vectors DataFrame
.add_dominant_taxas_to_colData <- function(x, dominances, name){
    dominances <- DataFrame(dominances)
    colnames(dominances) <- name
    colData(x)[,name] <- dominances
    x
}
