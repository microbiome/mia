#' Get Dominant Taxa
#'
#' This function returns the most dominant taxa.
#'
#' @param x
#' A \code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#' object
#'
#' @param abund_values A single character value for selecting the
#'   \code{\link[SummarizedExperiment:SummarizedExperiment-class]{assay}}
#'   to use for prevalence calculation.
#'
#' @param rank A single character defining a taxonomic rank. Must be a value
#' of the output of \code{taxonomicRanks()}.
#'
#' @param group With group, it is possible to group the observations. Must
#' be a one of the column names of colData.
#'
#' @param name A name for the column of the colData where the dominant taxa
#' will be stored in. Must be a string.
#'
#' @details
#' \code{getDominantTaxa} extracts the most abundant taxa
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
#' # Finds the dominant taxa.
#' x <- getDominantTaxa(x)
#' # Information is stored to colData
#' colData(x)
#'
#' # If taxonomic information is available, it is possible to find the most dominant
#' #group from specific taxonomic level, here family level. The name of column can be specified.
#' x <- getDominantTaxa(x, rank="Family", name="dominant_family")
#' colData(x)
#'
#' x <- microbiomeDataSets::dietswap()
#' # With group, it is possible to group observations based on group specified by user
#' x <- getDominantTaxa(x, group = "nationality")
#'
NULL

#' @rdname getDominantTaxa
#' @export
setGeneric("getDominantTaxa",signature = c("x"),
           function(x,
                    abund_values = "counts",
                    rank = NULL,
                    group = NULL,
                    name = "dominant_taxa")
               standardGeneric("getDominantTaxa"))


#' @rdname getDominantTaxa
#' @export
setMethod("getDominantTaxa", signature = c(x = "SummarizedExperiment"),
          function(x,
                   abund_values = "counts",
                   rank = NULL,
                   group = NULL,
                   name = "dominant_taxa"){

              # Input check
              # Check abund_values
              .check_abund_values(abund_values, x)

              # rank check
              if(!is.null(rank)){
                  if(!.is_a_string(rank)){
                      stop("'rank' must be an single character value.",
                           call. = FALSE)
                  }
                  .check_taxonomic_rank(rank, x)
              }

              # group check
              if(!is.null(group)){
                  if(isFALSE(any(group %in% colnames(colData(x))))){
                      stop("'group' variable must be in colnames(colData(x))")
                  }
              }

              # name check
              if(!.is_non_empty_string(name)){
                  stop("'name' must be a non-empty single character value.",
                       call. = FALSE)
              }

              # If "rank" is not NULL, species are aggregated according to the taxonomic rank
              # that is specified by user.
              if (!is.null(rank)) {

                  # Selects the level
                  col <- which( taxonomyRanks(x) %in% rank )

                  # Function from taxonomy.R. Divides taxas to groups where they belong.
                  tax_factors <- .get_tax_groups(x, col = col, onRankOnly = FALSE)

                  # Merges abundances within the groups
                  mat <- mergeRows(x, f = tax_factors)
              }

              # assays()$abund_values returns abundance of the taxa in the samples
              # apply() function finds the indices of taxa's that has the highest
              # abundance.
              # names() returns the names of taxa that are the most abundant.
              taxas <- names(mat)[apply(assay(mat, abund_values), 2, which.max)]

              # Adds taxa to colData
              mat <- .add_dominant_taxas_to_colData(mat, taxas, name)

              # Gets an overview
              overview <- .get_overview(mat, group, name)

              # Prints the overview
              message(overview)

              return(mat)

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

.get_overview <- function(x, group, name){

    # Creates a tibble df that contains dominant taxa and number of times that they present in samples
    # and relative portion of samples where they present.
    if (is.null(group)) {
        name <- sym(name)

        overview <- S4Vectors::as.data.frame(colData(x)) %>%
            dplyr::group_by(!!name) %>%
            dplyr::tally() %>%
            mutate(
                rel.freq = round(100 * n / sum(n), 1),
                rel.freq.pct = paste0(round(100 * n / sum(n), 0), "%")
            ) %>%
            dplyr::arrange(desc(n))
    } else {
        group <- sym(group)
        name <- sym(name)

        overview <- S4Vectors::as.data.frame(colData(x)) %>%
            dplyr::group_by(!!group, !!name) %>%
            dplyr::tally() %>%
            mutate(
                rel.freq = round(100 * n / sum(n), 1),
                rel.freq.pct = paste0(round(100 * n / sum(n), 0), "%")
            ) %>%
            dplyr::arrange(desc(n))
    }

    return(overview)
}

