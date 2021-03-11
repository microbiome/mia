#' Get Rare Taxa
#'
#' These functions return the rarest taxa.
#'
#' @param x a \code{\link{SummarizedExperiment}} object
#'
#' @param rank a single character defining a taxonomic rank. Must be a value
#' of the output of \code{taxonomicRanks()}. (default: \code{taxonomyRanks(x)[1L]})
#'
#' @param detection detection threshold for absence/presence. Either an
#'   absolutes value compared directly to the values of \code{x} or a relative
#'   value between 0 and 1, if \code{as_relative = TRUE}.
#'
#' @param prevalence prevalence threshold (in 0 to 1). The
#'   required prevalence is strictly less by default. To include the
#'   limit, set \code{include_lowest} to \code{TRUE}.
#'
#' @param as_relative logical scalar: Should the detection threshold be applied
#'   on compositional (relative) abundances? (default: \code{TRUE})
#'
#' @param include_highest logical scalar: Should the upper boundary of the
#'   detection and prevalence cutoffs be included? (default: \code{FALSE})
#'
#' @param ... additional parameters passed to \code{getRareTaxa}
#'
#' @details
#'   These functions return taxa, whose prevalence is under threshold.
#'
#' @return
#'   \code{getRareTaxa} returns a vector that contains the names of the rarest taxa.
#'   \code{getRare} returns \code{x} that includes only the rarest taxa.
#'
#' @seealso
#'   \code{\link[=getPrevalence]{getPrevalentTaxa}}
#'
#' @name getRareTaxa
#' @export
#'
#' @author Leo Lahti and Tuomas Borman. Contact: \url{microbiome.github.io}
#'
#' @examples
#' data(GlobalPatterns)
#'
#' # Detection threshold 20 (strictly lower by default);
#' # Note that the data (GlobalPatterns) is here in absolute counts
#' # (and not compositional, relative abundances)
#' # Prevalence threshold 50 percent (strictly lower by default)
#' taxa <- getRareTaxa(GlobalPatterns,
#'                          rank = "Phylum",
#'                          detection = 20,
#'                          prevalence = 50/100)
#' head(taxa)
#'
#' # Here data is compositional (as_relative = TRUE),
#' # and the upper limit is included (include_highest = TRUE)
#' taxa <- getRareTaxa(GlobalPatterns,
#'                          rank = "Family",
#'                          detection = 0,05,
#'                          prevalence = 0.1,
#'                          as_relative = TRUE,
#'                          include_highest = TRUE)
#' head(taxa)
#'
#' # Here SE object is returned. It includes only the rarest taxa.
#' taxa <- getRare(GlobalPatterns,
#'                          rank = "Order",
#'                          prevalence = 0.05)
#' taxa
#'
NULL

#' @rdname getRareTaxa
#' @export
setGeneric("getRareTaxa",signature = c("x"),
           function(x,
                    rank = taxonomyRanks(x)[1L],
                    detection = 0.1,
                    prevalence = 0.5,
                    as_relative = TRUE,
                    include_highest = FALSE,
                    ...)
               standardGeneric("getRareTaxa"))

#' @rdname getRareTaxa
#' @export
setGeneric("getRare",signature = c("x"),
           function(x, rank = taxonomyRanks(x)[1L], ...)
               standardGeneric("getRare"))



#' @rdname getRareTaxa
#' @export
setMethod("getRareTaxa", signature = c(x = "SummarizedExperiment"),
          function(x,
                   rank = taxonomyRanks(x)[1L],
                   detection = 0,
                   prevalence = 0.5,
                   as_relative = TRUE,
                   include_highest = FALSE,
                   ...){

              ######## Input check ########
              # rank must be one of the taxonomic ranks
              .check_taxonomic_rank(rank, x)

              # detection must be a number
              if (!.is_numeric_string(detection)) {
                  stop("'detection' must be a single numeric value or coercible to ",
                       "one.",
                       call. = FALSE)
              }
              # If value is a character, converts it to numeric value
              detection <- as.numeric(detection)

              # prevalence must be a number between 0-1
              if( !( .is_numeric_string(prevalence) && (prevalence >= 0 || prevalence <= 1) ) ){
                  stop("'prevalence' must be a single numeric value or coercible to ",
                       "one. Value must be between 0-1.",
                       call. = FALSE)
              }
              # If value is a character, converts it to numeric value
              prevalence <- as.numeric(prevalence)

              # as_relative must be a boolean value
              if( !.is_a_bool(as_relative) ){
                  stop("'as_relative' must be TRUE or FALSE.", call. = FALSE)
              }

              # include_highest must be a boolean value
              if( !.is_a_bool(include_highest) ){
                  stop("'include_highest' must be TRUE or FALSE.", call. = FALSE)
              }
              ###### Input check end ######

              # In getPrevalentTaxa, 'include_lowest'  means that the limit is included
              # to the most prevalent taxa, which means that it is excluded from the rarest
              # taxa. To include limit to the rarest taxa, 'include_lowest' must be the opposite.
              include_lowest <- !include_highest

              # Agglomerates object by rank
              x <- agglomerateByRank(x, rank = rank)

              # Stores rownames
              all_taxa_names <- rownames(x)

              # Fetches names of core taxa
              prevalent_taxa_names <- getPrevalentTaxa(x,
                                                       rank = rank,
                                                       detection = detection,
                                                       prevalence = prevalence,
                                                       as_relative = as_relative,
                                                       include_lowest = include_lowest)

              # Takes complement by removing core taxa from all the taxa
              rare_taxa_names <- setdiff(all_taxa_names, prevalent_taxa_names)

              return(rare_taxa_names)

          }
)

#' @rdname getRareTaxa
#' @export
setMethod("getRare", signature = c(x = "SummarizedExperiment"),
          function(x, rank = taxonomyRanks(x)[1L], ...){

              ####### Input check #########
              # rank must be one of the taxonomic ranks
              .check_taxonomic_rank(rank, x)
              ###### Input check end ######

              # Agglomerates object by rank
              x <- agglomerateByRank(x, rank = rank)

              # Gets the names of rare taxa
              rare_taxa_names <- getRareTaxa(x, rank = rank, ...)

              # Subsets the data
              rare_x <- x[rare_taxa_names]

              return(rare_x)

          }
)
