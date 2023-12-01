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
#' @param assay.type A single character value for selecting the
#'   \code{\link[SummarizedExperiment:SummarizedExperiment-class]{assay}}
#'   to use for identifying dominant taxa.
#'
#' @param assay_name a single \code{character} value for specifying which
#'   assay to use for calculation.
#'   (Please use \code{assay.type} instead. At some point \code{assay_name}
#'   will be disabled.)
#'   
#' @param rank A single character defining a taxonomic rank. Must be a value of
#'   the output of \code{taxonomyRanks()}.
#'
#' @param name A name for the column of the \code{colData} where the dominant
#'   taxa will be stored in when using \code{addPerSampleDominantFeatures}.
#'   
#' @param complete 
#' 
#' @param n
#' 
#' @param complete
#'
#' @param ... Additional arguments passed on to \code{agglomerateByRank()} when
#' \code{rank} is specified.
#'
#' @details
#' \code{addPerSampleDominantFeatures} extracts the most abundant taxa in a
#' \code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#' object, and stores the information in the \code{colData}. It is a wrapper for
#' \code{perSampleDominantFeatures}.
#'
#' With \code{rank} parameter, it is possible to agglomerate taxa based on taxonomic
#' ranks. E.g. if 'Genus' rank is used, all abundances of same Genus are added
#' together, and those families are returned. See \code{agglomerateByRank()} for
#' additional arguments to deal with missing values or special characters.
#'
#' @return \code{perSampleDominantFeatures} returns a named character vector \code{x}
#' while \code{addPerSampleDominantFeatures} returns
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
#' sim.dom <- perSampleDominantFeatures(x, rank="Genus")
#'
#' # Add information to colData
#' x <- addPerSampleDominantFeatures(x, rank = "Genus", name="dominant_genera")
#' colData(x)
NULL

#' @rdname perSampleDominantTaxa
#' @aliases perSampleDominantTaxa
#' @export
setGeneric("perSampleDominantFeatures",signature = c("x"),
           function(x, assay.type = assay_name, assay_name = "counts", 
                    rank = NULL, other_name = "Other", n = NULL, complete = FALSE, ...)
               standardGeneric("perSampleDominantFeatures"))

#' @rdname perSampleDominantTaxa
#' @aliases perSampleDominantTaxa
#' @importFrom IRanges relist microbiome
#' @export
setMethod("perSampleDominantFeatures", signature = c(x = "SummarizedExperiment"),
    function(x, assay.type = assay_name, assay_name = "counts", 
             rank = NULL, other_name = "Other", n = NULL, complete = TRUE, ...){
        # Input check
        # Check assay.type
        .check_assay_present(assay.type, x)
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
            mat <- assay(x, assay.type)
        } # Otherwise, if "rank" is NULL, abundances are stored without ranking
        else {
            mat <- assay(x, assay.type)
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
        
        # If individual sample contains multiple dominant taxa (they have equal counts) and if 
        # complete is FALSE, the an arbitrarily chosen dominant taxa is returned
        if( length(taxa) > nrow(colData(x))){
            # Store order
            order <- unique(names(taxa))
            # there are multiple dominant taxa in one sample (counts are equal), length
            # of dominant is greater than rows in colData.
            taxa <- split(taxa, rep(names(taxa), lengths(taxa)) )
            # Order the data
            taxa <- taxa[order]
            
            # If complete is set FALSE, and there are multiple dominant taxa, one of them is arbitrarily chosen
            if( complete == FALSE) {
                # Randomly choose one element from each sublist
                taxa <- lapply(taxa, function(item) {
                                if (length(item) > 1) {
                                # If it's a list with more than one element, choose one randomly
                                    return(sample(item, 1))
                                } else {
                                    return(item)
                                }
                            })
                taxa <- unname(sapply(taxa, function (x) {unlist(x)}))
                warning("Multiple dominant taxa were found for some samples. Run perSampleDominantFeatures(x, complete = TRUE) for details.")
            }
        }     
        
        # Name "Other" the features that are not included in n the most abundant in the data
        if(!is.null(n)) {
            flat_taxa <- unlist(taxa, recursive = TRUE)
            top <- microbiome::top(flat_taxa, n=n)
            # Group the rest into the "Other" category
            
            # complete parameter defines the output format
            if (complete) { # return a list
                taxa <- lapply(taxa, function(item) {
                        replace <- !(item %in% names(top))
                        item[replace] <- other_name
                        return(item)
                })
            } else { # return a vector
                taxa[!taxa %in% names(top)] <- other_name
                taxa <- unname(sapply(taxa, function (x) {unlist(x)}))
            }
        }
        return(taxa)
    }
)

#' @rdname perSampleDominantTaxa
#' @aliases perSampleDominantFeatures
#' @export
setGeneric("perSampleDominantTaxa", signature = c("x"),
        function(x, ...) 
            standardGeneric("perSampleDominantTaxa"))

#' @rdname perSampleDominantTaxa
#' @aliases perSampleDominantFeatures
#' @export
setMethod("perSampleDominantTaxa", signature = c(x = "SummarizedExperiment"),
        function(x, ...){
            .Deprecated(old ="perSampleDominantTaxa", new = "perSampleDominantFeatures", msg = "The 'perSampleDominantTaxa' function is deprecated. Use 'perSampleDominantFeatures' instead.")
            perSampleDominantFeatures(x, ...)
        }
)

#' @rdname perSampleDominantTaxa
#' @aliases addPerSampleDominantTaxa
#' @export
setGeneric("addPerSampleDominantFeatures", signature = c("x"),
           function(x, name = "dominant_taxa", other_name = "Other", n = NULL, ...)
               standardGeneric("addPerSampleDominantFeatures"))

#' @rdname perSampleDominantTaxa
#' @aliases addPerSampleDominantTaxa
#' @export
setMethod("addPerSampleDominantFeatures", signature = c(x = "SummarizedExperiment"),
    function(x, name = "dominant_taxa", other_name = "Other", n = NULL, complete = FALSE, ...) {
        # name check
        if(!.is_non_empty_string(name)){
            stop("'name' must be a non-empty single character value.",
                 call. = FALSE)
        }
        # other_name check
        if(!.is_non_empty_string(other_name)){
            stop("'other_name' must be a non-empty single character value.",
                 call. = FALSE)
        }
        if(complete) {
            warning("Multiple dominant taxa are not allowed for samples. Run perSampleDominantFeatures(x, complete = TRUE) for details. ")
            complete <- FALSE
        }
        # complete must be FALSE to fit the dominant feature vector in colData
        dom.taxa <- perSampleDominantFeatures(x, other_name = other_name, n = n, complete = complete, ...)
        
        colData(x)[[name]] <- dom.taxa
        return(x)
    }
)

#' @rdname perSampleDominantTaxa
#' @aliases addPerSampleDominantFeatures
#' @export
setGeneric("addPerSampleDominantTaxa", signature = c("x"),
        function(x, ...) 
            standardGeneric("addPerSampleDominantTaxa"))

#' @rdname perSampleDominantTaxa
#' @aliases addPerSampleDominantFeatures
#' @export
setMethod("addPerSampleDominantTaxa", signature = c(x = "SummarizedExperiment"),
        function(x, ...){
            .Deprecated(old ="addPerSampleDominantTaxa", new = "addPerSampleDominantFeatures", msg = "The 'addPerSampleDominantTaxa' function is deprecated. Use 'addPerSampleDominantFeatures' instead.")
            addPerSampleDominantFeatures(x, ...)
        }
)
