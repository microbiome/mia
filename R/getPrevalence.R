#' Calculation prevalence information for features across samples
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
#' @param other_label A single \code{character} valued used as the label for the
#'   summary of non-prevalent taxa. (default: \code{other_label = "Other"})
#'
#' @param rank,... additional arguments
#' \itemize{
#'   \item{If \code{!is.null(rank)} arguments are passed on to
#'     \code{\link[=agglomerate-methods]{agglomerateByRank}}. See
#'     \code{\link[=agglomerate-methods]{?agglomerateByRank}} for more details.
#'   }
#'   \item{for \code{getPrevalentAbundance} additional parameters passed to
#'     \code{getPrevalentTaxa}}
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
#' The core abundance index from \code{getPrevalentAbundance} gives the relative
#' proportion of the core species (in between 0 and 1). The core taxa are
#' defined as those that exceed the given population prevalence threshold at the
#' given detection level as set for \code{getPrevalentTaxa}.
#' 
#' \code{subsetPrevalentTaxa} and \code{subsetRareTaxa} return a subset of \code{x}. 
#' The subset includes the most prevalent or rare taxa that are calculated with 
#' \code{getPrevalentTaxa} or \code{getRareTaxa} respectively.
#'
#' @return
#' \code{subsetPrevalentTaxa} and \code{subsetRareTaxa} return subset of \code{x}.
#' 
#' All other functions return a named \code{numeric} vector. For \code{getPrevalence} 
#' the names are either the row names of \code{x} or the names after agglomeration.
#'
#' For \code{getPrevalentAbundance} the names correspond to the column name
#' names of \code{x} and include the joint abundance of prevalent taxa.
#'
#' For \code{getPrevalentTaxa} only the names exceeding the threshold set by
#' \code{prevalence} are returned.
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
#' @author
#' Leo Lahti
#' For \code{getPrevalentAbundance}: Leo Lahti and Tuomas Borman.
#' Contact: \url{microbiome.github.io}
#'
#' @export
#'
#' @examples
#' data(GlobalPatterns)
#' tse <- GlobalPatterns
#' # Get prevalence estimates for individual ASV/OTU
#' prevalence.frequency <- getPrevalence(tse,
#'                                       detection = 0,
#'                                       sort = TRUE,
#'                                       as_relative = TRUE)
#' head(prevalence.frequency)
#'
#' # Get prevalence estimates for phylums
#' # - the getPrevalence function itself always returns population frequencies
#' prevalence.frequency <- getPrevalence(tse,
#'                                       rank = "Phylum",
#'                                       detection = 0,
#'                                       sort = TRUE,
#'                                       as_relative = TRUE)
#' head(prevalence.frequency)
#'
#' # - to obtain population counts, multiply frequencies with the sample size,
#' # which answers the question "In how many samples is this phylum detectable"
#' prevalence.count <- prevalence.frequency * ncol(tse)
#' head(prevalence.count)
#'
#' # Detection threshold 1 (strictly greater by default);
#' # Note that the data (GlobalPatterns) is here in absolute counts
#' # (and not compositional, relative abundances)
#' # Prevalence threshold 50 percent (strictly greater by default)
#' prevalent <- getPrevalentTaxa(tse,
#'                               rank = "Phylum",
#'                               detection = 10,
#'                               prevalence = 50/100,
#'                               as_relative = FALSE)
#' head(prevalent)
#' 
#' # Gets a subset of object that includes prevalent taxa
#' tse_prevalent <- subsetByPrevalentTaxa(tse,
#'                                        rank = "Family",
#'                                        detection = 0.001,
#'                                        prevalence = 0.55,
#'                                        as_relative = TRUE)
#' 
#' # Stores the subset to the original object as an alternative experiment
#' altExp(tse, "prevalent") <- tse_prevalent    
#' altExp(tse, "prevalent")                                 
#'
#' # getRareTaxa returns the inverse
#' rare <- getRareTaxa(tse,
#'                     rank = "Phylum",
#'                     detection = 1/100,
#'                     prevalence = 50/100,
#'                     as_relative = TRUE)
#' head(rare)
#' 
#' # Gets a subset of object that includes rare taxa
#' tse_rare <- subsetByRareTaxa(tse,
#'                              rank = "Class",
#'                              detection = 0.001,
#'                              prevalence = 0.001,
#'                              as_relative = TRUE)
#'                              
#' # Stores the subset to the original object as an alternative experiment
#' altExp(tse, "rare") <- tse_rare
#' altExp(tse, "rare")      
#' 
#' # Names of both experiments, prevalent and rare, can be found from slot altExpNames
#' tse
#'                          
#' data(esophagus)
#' getPrevalentAbundance(esophagus, abund_values = "counts")
#'
#' # data can be aggregated based on prevalent taxonomic results
#' agglomerateByPrevalence(tse,
#'                         rank = "Phylum",
#'                         detection = 1/100,
#'                         prevalence = 50/100,
#'                         as_relative = TRUE)
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
            stop("'detection' must be a single numeric value or coercible to ",
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
        # # If prev doesn't include names / if rownames(x) was NULL
        # if( is.null(names(prev)) ){
        #     # Adds indices as names
        #     names(prev) <- factor(1:length(prev))
        # }
        # Always return prevalence as a relative frequency.
        # This helps to avoid confusion with detection limit
        prev <- prev / ncol(x)
        if (sort) {
            prev <- rev(sort(prev))
        }
        prev
    }
)

.agg_for_prevalence <- function(x, rank, relabel = FALSE, make_unique = TRUE,
                                ...){
    if(!is.null(rank)){
        .check_taxonomic_rank(rank, x)
        args <- c(list(x = x, rank = rank), list(...))
        if(is.null(args[["na.rm"]])){
            args[["na.rm"]] <- TRUE
        }
        argNames <- c("x","rank","onRankOnly","na.rm","empty.fields",
                      "archetype","mergeTree","average","BPPARAM")
        args <- args[names(args) %in% argNames]
        x <- do.call(agglomerateByRank, args)
        if(relabel){
            rownames(x) <- getTaxonomyLabels(x, make_unique = make_unique)
        }
    }
    x
}

#' @rdname getPrevalence
#' @export
setMethod("getPrevalence", signature = c(x = "SummarizedExperiment"),
    function(x, abund_values = "counts", as_relative = TRUE,
             rank = NULL, ...){

        if(!.is_a_bool(as_relative)){
            stop("'as_relative' must be TRUE or FALSE.", call. = FALSE)
        }

        # check assay
        .check_assay_present(abund_values, x)
        x <- .agg_for_prevalence(x, rank = rank, ...)
        mat <- assay(x, abund_values)
        if (as_relative) {
            mat <- .calc_rel_abund(mat)
        }
        getPrevalence(mat, ...)
    }
)
############################# getPrevalentTaxa #################################
#' @rdname getPrevalence
#'
#' @param prevalence Prevalence threshold (in 0 to 1). The
#'   required prevalence is strictly greater by default. To include the
#'   limit, set \code{include_lowest} to \code{TRUE}.
#'
#' @details
#' \code{getPrevalentTaxa} returns taxa that are more prevalent with the
#' given detection threshold for the selected taxonomic rank.
#'
#' @export
setGeneric("getPrevalentTaxa", signature = "x",
           function(x, ...)
               standardGeneric("getPrevalentTaxa"))

.get_prevalent_taxa <- function(x, rank, prevalence = 50/100,
                                include_lowest = FALSE, ...){
    # input check
    if (!.is_numeric_string(prevalence)) {
        stop("'prevalence' must be a single numeric value or coercible to ",
             "one.",
             call. = FALSE)
    }

    prevalence <- as.numeric(prevalence)
    if(!.is_a_bool(include_lowest)){
        stop("'include_lowest' must be TRUE or FALSE.", call. = FALSE)
    }

    if(!is(x,"SummarizedExperiment")){
        pr <- getPrevalence(x, ...)
    } else {
        pr <- getPrevalence(x, rank = rank, ...)
    }

    if (include_lowest) {
        taxa <- pr >= prevalence
    } else {
        taxa <- pr > prevalence
    }
    
    # Gets indices of most prevalent taxa
    indices  <- which(taxa)
    # If vector contains names
    if( !is.null(names(indices)) ){
        # Gets the names
        taxa <- names(indices)
    } else {
        # Otherwise indices are returned
        taxa <- indices
    }
    taxa
}

#' @rdname getPrevalence
#' @export
setMethod("getPrevalentTaxa", signature = c(x = "ANY"),
    function(x, prevalence = 50/100, include_lowest = FALSE, ...){
        .get_prevalent_taxa(x, rank = NULL, prevalence = prevalence,
                            include_lowest = include_lowest, ...)
    }
)

#' @rdname getPrevalence
#' @export
setMethod("getPrevalentTaxa", signature = c(x = "SummarizedExperiment"),
    function(x, prevalence = 50/100, rank = NULL,
             include_lowest = FALSE, ...){
        .get_prevalent_taxa(x, rank = rank, prevalence = prevalence,
                            include_lowest = include_lowest, ...)
    }
)

############################# getRareTaxa ######################################

#' @rdname getPrevalence
#'
#' @details
#' \code{getRareTaxa} returns complement of \code{getPrevalentTaxa}.
#'
#' @export
setGeneric("getRareTaxa", signature = "x",
           function(x, rank = NULL, ...)
               standardGeneric("getRareTaxa"))

#' @rdname getPrevalence
#' @export
setMethod("getRareTaxa", signature = c(x = "SummarizedExperiment"),
    function(x, rank = NULL, ...){
        # Gets the prevalent taxa
        prev_taxa <- getPrevalentTaxa(x, rank = rank, ...)
        # If agglomeration was done
        if( !is.null(rank) ){
            # Gets names from specified taxonomic level
            taxa <- rowData(x)[[rank]]
        } else {
            # If rownames is not NULL
            if( !is.null(rownames(x)) ){
                # Gets rownames 
                taxa <- rownames(x)
            } else {
                # Gets indices of taxa
                taxa <- 1:nrow(tse)
            }
        }
        unique(taxa[!is.na(taxa) & !(taxa %in% prev_taxa)])
    }
)

############################# subsetByPrevalentTaxa ############################

#' @rdname getPrevalence
#' @export
setGeneric("subsetByPrevalentTaxa", signature = "x",
           function(x, ...)
               standardGeneric("subsetByPrevalentTaxa"))

#' @rdname getPrevalence
#' @export
setMethod("subsetByPrevalentTaxa", signature = c(x = "SummarizedExperiment"),
    function(x, rank = NULL, ...){
        # If rank is not NULL
        if( !is.null(rank) ){
            # Checks rank
            .check_taxonomic_rank(rank, x)
            # Agglomerates the object by rank. Taxa that do not have information
            # at specific rank are excluded.
            x <- agglomerateByRank(x, rank = rank, na.rm = TRUE)
            # Changes rank to NULL, so that agglomeration is not done in next steps
            rank <- NULL
        }
        # Gets the prevalent taxa
        prevalent_taxa <- getPrevalentTaxa(x, ...)
        # Subsets the object based on prevalent taxa
        x <- x[prevalent_taxa, ]
    }
)

############################# subsetByRareTaxa #################################

#' @rdname getPrevalence
#' @export
setGeneric("subsetByRareTaxa", signature = "x",
           function(x, ...)
               standardGeneric("subsetByRareTaxa"))

#' @rdname getPrevalence
#' @export
setMethod("subsetByRareTaxa", signature = c(x = "SummarizedExperiment"),
    function(x, rank = NULL, ...){
        # If rank is not NULL
        if( !is.null(rank) ){
            # Checks rank
            .check_taxonomic_rank(rank, x)
            # Agglomerates the object by rank. Taxa that do not have information
            # at specific rank are excluded.
            x <- agglomerateByRank(x, rank = rank, na.rm = TRUE)
            # Changes rank to NULL, so that agglomeration is not done in next steps
            rank <- NULL
        }
        # Gets the rare taxa
        rare_taxa <- getRareTaxa(x, ...)
        # Subsets the object based on rare taxa
        x <- x[rare_taxa, ]
    }
)

############################# getPrevalentAbundance ############################

#' @rdname getPrevalence
#' @export
setGeneric("getPrevalentAbundance", signature = "x",
           function(x, abund_values = "relabundance", ...)
               standardGeneric("getPrevalentAbundance"))

#' @rdname getPrevalence
#' @export
setMethod("getPrevalentAbundance", signature = c(x = "ANY"),
    function(x, ...){
        x <- .calc_rel_abund(x)
        cm <- getPrevalentTaxa(x, ...)
        if (length(cm) == 0) {
            stop("With the given abundance and prevalence thresholds, no taxa ",
                 "were found. Try to change detection and prevalence ",
                 "parameters.",
                 call. = FALSE)
        }
        colSums(x[cm, ,drop=FALSE])
    }
)

#' @rdname getPrevalence
#' @export
setMethod("getPrevalentAbundance", signature = c(x = "SummarizedExperiment"),
    function(x, abund_values = "counts", ...){
        # check assay
        .check_assay_present(abund_values, x)
        #
        getPrevalentAbundance(assay(x,abund_values))
    }
)


############################# agglomerateByPrevalence ##########################

#' @rdname getPrevalence
#' @export
setGeneric("agglomerateByPrevalence", signature = "x",
           function(x, ...)
               standardGeneric("agglomerateByPrevalence"))


#' @rdname getPrevalence
#' @export
setMethod("agglomerateByPrevalence", signature = c(x = "SummarizedExperiment"),
    function(x, rank = taxonomyRanks(x)[1L], other_label = "Other", ...){
        # input check
        if(!.is_a_string(other_label)){
            stop("'other_label' must be a single character value.",
                 call. = FALSE)
        }
        #
        x <- .agg_for_prevalence(x, rank, ...)
        pr <- getPrevalentTaxa(x, rank = NULL, ...)
        f <- rownames(x) %in% pr
        if(any(!f)){
            other_x <- mergeRows(x[!f,], factor(rep(1L,sum(!f))))
            rowData(other_x)[,colnames(rowData(other_x))] <- NA
            # set the other label
            rownames(other_x) <- other_label
            if(!is.null(rank)){
                rowData(other_x)[,rank] <- other_label
            }
            # temporary fix until TSE supports rbind
            class <- c("SingleCellExperiment","RangedSummarizedExperiment",
                       "SummarizedExperiment")
            class_x <- class(x)
            if(!(class_x %in% class)){
                class <- "SingleCellExperiment"
            } else {
                class <- class[class_x == class]
            }
            x <- rbind(as(x[f,],class),
                       as(other_x,class))
            x <- as(x,class_x)
        }
        x
    }
)

