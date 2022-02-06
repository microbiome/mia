#' Get dominant taxa
#'
#' These functions return information about the most dominant taxa in a
#' \code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#' object.
#'
#' @param x A
#'   \code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#'   object.
#'
#' @param abund_values A single character value for selecting the
#'   \code{\link[SummarizedExperiment:SummarizedExperiment-class]{assay}}
#'   to use for identifying dominant taxa.
#'
#' @param rank A single character defining a taxonomic rank. Must be a value of
#'   the output of \code{taxonomyRanks()}.
#'
#' @param name A name for the column of the \code{colData} where the dominant
#'   taxa will be stored in when using \code{addPerSampleDominantTaxa}.
#'
#' @param ... Additional arguments passed on to \code{agglomerateByRank()} when
#' \code{rank} is specified.
#'
#' @details
#' \code{addPerSampleDominantTaxa} extracts the most abundant taxa in a
#' \code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#' object, and stores the information in the \code{colData}. It is a wrapper for
#' \code{perSampleDominantTaxa}.
#'
#' With \code{rank} parameter, it is possible to agglomerate taxa based on taxonomic
#' ranks. E.g. if 'Genus' rank is used, all abundances of same Genus are added
#' together, and those families are returned. See \code{agglomerateByRank()} for
#' additional arguments to deal with missing values or special characters.
#'
#' @return \code{perSampleDominantTaxa} returns a named character vector \code{x}
#' while \code{addPerSampleDominantTaxa} returns
#' \code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#' with additional column in \code{\link{colData}} named \code{*name*}.
#'
#' @name perSampleDominantTaxa
#' @export
#'
#' @author Leo Lahti, Tuomas Borman and Sudarshan A. Shetty.
#'
#' @examples
#' data(GlobalPatterns)
#' x <- GlobalPatterns
#'
#' # Finds the dominant taxa.
#' sim.dom <- perSampleDominantTaxa(x, rank="Genus")
#'
#' # Add information to colData
#' x <- addPerSampleDominantTaxa(x, rank = "Genus", name="dominant_genera")
#' colData(x)
NULL

#' @rdname perSampleDominantTaxa
#' @export
setGeneric("perSampleDominantTaxa",signature = c("x"),
           function(x, abund_values = "counts", rank = NULL, ...)
               standardGeneric("perSampleDominantTaxa"))

#' @rdname perSampleDominantTaxa
#' @importFrom IRanges relist
#' @export
setMethod("perSampleDominantTaxa", signature = c(x = "SummarizedExperiment"),
    function(x, abund_values = "counts", rank = NULL, ...){
        # Input check
        # Check abund_values
        .check_assay_present(abund_values, x)
        # rank check
        if(!is.null(rank)){
            if(!.is_a_string(rank)){
                stop("'rank' must be an single character value.",
                     call. = FALSE)
            }
            .check_taxonomic_rank(rank, x)
        }
        # If "rank" is not NULL, species are aggregated according to the
        # taxonomic rank that is specified by user.
        if (!is.null(rank)) {
            x <- agglomerateByRank(x, rank, ...)
            mat <- assay(x, abund_values)
        } # Otherwise, if "rank" is NULL, abundances are stored without ranking
        else {
            mat <- assay(x, abund_values)
        }
        # apply() function finds the indices of taxa's that has the highest
        # abundance.
        # rownames() returns the names of taxa that are the most abundant.
        idx <- as.list(apply(t(mat) == colMaxs(mat),1L,which))
        # Get rownames based on indices
        taxa <- rownames(mat)[unlist(idx)]
        # If multiple dominant taxa were found, names contain taxa in addition to 
        # sample name. Names are converted so that they include only sample names.
        names(taxa) <- rep( names(idx), times = lengths(idx) )
        return(taxa)
    }
)


#' @rdname perSampleDominantTaxa
#' @export
setGeneric("addPerSampleDominantTaxa", signature = c("x"),
           function(x, name = "dominant_taxa", ...)
               standardGeneric("addPerSampleDominantTaxa"))

#' @rdname perSampleDominantTaxa
#' @export
setMethod("addPerSampleDominantTaxa", signature = c(x = "SummarizedExperiment"),
    function(x, name = "dominant_taxa", ...){
        # name check
        if(!.is_non_empty_string(name)){
            stop("'name' must be a non-empty single character value.",
                 call. = FALSE)
        }
        dom.taxa <- perSampleDominantTaxa(x, ...)
        # If individual sample contains multiple dominant taxa (they have equal counts)
        if( length(dom.taxa) > nrow(colData(x)) ){
            # Store order
            order <- unique(names(dom.taxa))
            # there are multiple dominant taxa in one sample (counts are equal), length
            # of dominant is greater than rows in colData. --> create a list that contain
            # dominant taxa, and is as long as there are rows in colData
            dom.taxa <- split(dom.taxa, rep(names(dom.taxa), lengths(dom.taxa)) )
            # Order the data
            dom.taxa <- dom.taxa[order]
        }
        colData(x)[[name]] <- dom.taxa
        return(x)
    }
)
