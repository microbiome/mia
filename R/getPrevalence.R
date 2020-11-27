#' Prevalence calculation
#'
#' These functions calculate the population prevalence for taxonomic ranks in a
#' \code{\link{SummarizedExperiment-class}} object.
#'
#' @param x a
#'   \code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#'   object
#'
#' @param detection Detection threshold for absence/presence. Either an
#'   absolutes value compared directly to the values of \code{x} or a relative
#'   value between 0 and 1, if \code{as_relative = TRUE}.
#'
#' @param include_lowest logical scalar: Should the lower boundary of the
#'   detection and prevalence cutoffs be included? (default: \code{FALSE})
#'
#' @param sort logical scalar: Should the result be sorted by prevalence?
#'   (default: \code{FALSE})
#'
#' @param as_relative logical scalar: Should the detection threshold be applied
#'   on compositional (relative) abundances? (default: \code{TRUE})
#'
#' @param abund_values A single character value for selecting the
#'   \code{\link[SummarizedExperiment:SummarizedExperiment-class]{assay}}
#'   to use for prevalence calculation.
#'
#' @param rank,... additional arguments
#' \itemize{
#'   \item{If \code{!is.null(rank)} arguments are passed on to
#'   \code{\link[=agglomerate-methods]{agglomerateByRank}}. See
#'   \code{\link[=agglomerate-methods]{?agglomerateByRank}} for more details.
#'   }
#' }
#'
#' @details
#' \code{getPrevalence} calculates the relative frequency of samples that exceed
#' the detection threshold. For \code{SummarizedExperiment} objects, the
#' prevalence is calculated for the selected taxonomic rank, otherwise for the
#' rows. The absolute population prevalence can be obtained by multiplying the
#' prevalence by the number of samples (\code{ncol(x)}). If \code{as_relative =
#' TRUE} the relative frequency (between 0 and 1) is used to check against the
#' \code{detection} threshold.
#'
#' @return
#' a named \code{numeric} vector. The names are either the row names of \code{x}
#' or the names after agglomeration. For \code{getPrevalentTaxa} only the names
#' exceeding the threshold set by \code{prevalence} are returned.
#'
#' @seealso
#' \code{\link[=agglomerate-methods]{agglomerateByRank}},
#' \code{\link[=getTopTaxa]{getTopTaxa}}
#'
#'
#' @name getPrevalence
#' @export
#'
#' @references
#' A Salonen et al. The adult intestinal core microbiota is determined by
#' analysis depth and health status. Clinical Microbiology and Infection
#' 18(S4):16 20, 2012.
#' To cite the R package, see citation('mia')
#'
#' @author Leo Lahti
#'
#' @export
#'
#' @examples
#' data(GlobalPatterns)
#' # Get prevalence estimates for individual ASV/OTU
#' prevalence.frequency <- getPrevalence(GlobalPatterns,
#'                                       detection = 0,
#'                                       sort = TRUE,
#'                                       as_relative = TRUE)
#' head(prevalence.frequency)
#'
#' # Get prevalence estimates for phylums
#' # - the getPrevalence function itself always returns population frequencies
#' # - to obtain population counts, multiply frequencies with the sample size,
#' #   which answers the question "In how many samples is this phylum detectable"
#' prevalence.frequency <- getPrevalence(GlobalPatterns,
#'                                       rank = "Phylum",
#'                                       detection = 0,
#'                                       sort = TRUE,
#'                                       as_relative = TRUE)
#' head(prevalence.frequency)
#' prevalence.count <- prevalence.frequency * ncol(GlobalPatterns)
#' head(prevalence.count)
#'
#' # Detection threshold 1 (strictly greater by default);
#' # Note that the data (GlobalPatterns) is here in absolute counts
#' # (and not compositional, relative abundances)
#' # Prevalence threshold 50 percent (strictly greater by default)
#' taxa <- getPrevalentTaxa(GlobalPatterns,
#'                          rank = "Phylum",
#'                          detection = 1/100,
#'                          prevalence = 50/100,
#'                          as_relative = TRUE)
#' head(taxa)
NULL

#' @rdname getPrevalence
#' @export
setGeneric("getPrevalence", signature = "x",
           function(x, ...)
               standardGeneric("getPrevalence"))

#' @rdname getPrevalence
#' @export
setMethod("getPrevalence", signature = c(x = "ANY"),
    function(x, detection = 0, include_lowest = FALSE, sort = FALSE, ...){
        # input check
        if (!.is_numeric_string(detection)) {
            stop("'detection' must be ainslge numeric value or coercibel to ",
                 "one.",
                 call. = FALSE)
        }
        detection <- as.numeric(detection)
        if(!.is_a_bool(include_lowest)){
            stop("'include_lowest' must be TRUE or FALSE.", call. = FALSE)
        }
        if(!.is_a_bool(sort)){
            stop("'sort' must be TRUE or FALSE.", call. = FALSE)
        }
        #
        if (include_lowest) {
            prev <- x >= detection
        } else {
            prev <- x > detection
        }
        prev <- rowSums(prev)
        # Always return prevalence as a relative frequency.
        # This helps to avoid confusion with detection limit
        prev <- prev / ncol(x)
        if (sort) {
            prev <- rev(sort(prev))
        }
        prev
    }
)

#' @rdname getPrevalence
#' @export
setMethod("getPrevalence", signature = c(x = "SummarizedExperiment"),
    function(x, abund_values = "counts", as_relative = TRUE,
             rank = NULL, ...){
        if(!.is_a_bool(as_relative)){
            stop("'as_relative' must be TRUE or FALSE.", call. = FALSE)
        }

        # check assay
        .check_abund_values(abund_values, x)
        if(!is.null(rank)){
            x <- agglomerateByRank(x, rank = rank, ...)
        }
        mat <- assay(x, abund_values)
        if (as_relative) {
            mat <- .calc_rel_abund(mat)
        }
        # retrieve abundance matrix
        getPrevalence(mat, ...)
    }
)

#' @rdname getPrevalence
#'
#' @param prevalence Prevalence threshold (in 0 to 1). The
#' required prevalence is strictly greater by default. To include the
#' limit, set include_lowest to TRUE.
#'
#' @details
#' \code{getPrevalentTaxa} returns taxa that are more prevalent with the
#' given detection threshold for the selected taxonomic rank.
#'
#' @export
setGeneric("getPrevalentTaxa", signature = "x",
           function(x, ...)
               standardGeneric("getPrevalentTaxa"))

#' @rdname getPrevalence
#' @export
setMethod("getPrevalentTaxa", signature = c(x = "SummarizedExperiment"),
    function(x, prevalence = 50/100, rank = taxonomyRanks(x)[1L],
             include_lowest = FALSE, ...){

        # input check
        if (!.is_numeric_string(prevalence)) {
            stop("'prevalence' must be a single numeric value or coercibel to ",
                 "one.",
                 call. = FALSE)
        }

        prevalence <- as.numeric(prevalence)
        if(!.is_a_bool(include_lowest)){
            stop("'include_lowest' must be TRUE or FALSE.", call. = FALSE)
        }

        pr <- getPrevalence(x, rank = rank, ...)

        if (include_lowest) {
          taxa <- pr >= prevalence
        } else {
          taxa <- pr > prevalence
        }

        taxa <- names(which(taxa))
        taxa
    }
)
