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
#'   \item{for \code{getPrevalentTaxa}, \code{getRareaxa}, 
#'     \code{subsetByPrevalentTaxa} and \code{subsetByRareaxa} additional 
#'     parameters passed to \code{getPrevalence}}
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
#' All other functions return a named vectors:
#' \itemize{
#'   \item{\code{getPrevalence} returns a \code{numeric} vector with the 
#'     names being set to either the row names of \code{x} or the names after 
#'     agglomeration.}
#'   \item{\code{getPrevalentAbundance} returns a \code{numeric} vector with
#'     the names corresponding to the column name of \code{x} and include the 
#'     joint abundance of prevalent taxa.}
#'   \item{\code{getPrevalentTaxa} and \code{getRareTaxa} return a 
#'     \code{character} vector with only the names exceeding the threshold set
#'     by \code{prevalence}, if the \code{rownames} of \code{x} is set. 
#'     Otherwise an \code{integer} vector is returned matching the rows in
#'     \code{x}.}
#' }
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
#' altExp(tse, "prevalent") <- subsetByPrevalentTaxa(tse,
#'                                        rank = "Family",
#'                                        detection = 0.001,
#'                                        prevalence = 0.55,
#'                                        as_relative = TRUE)
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
#' altExp(tse, "rare") <- subsetByRareTaxa(tse,
#'                              rank = "Class",
#'                              detection = 0.001,
#'                              prevalence = 0.001,
#'                              as_relative = TRUE)
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
        # input check
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

.norm_rownames <- function(x){
    if(is.null(rownames(x))){
        rownames(x) <- seq_len(nrow(x))
    } else if(anyDuplicated(rownames(x))) {
        rownames(x) <- make.unique(rownames(x))
    }
    x
}

.get_prevalent_indices <- function(x, prevalence = 50/100,
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
    # rownames must bet set and unique, because if sort = TRUE, the order is 
    # not preserved
    x <- .norm_rownames(x)
    pr <- getPrevalence(x, rank = NULL, ...)
    
    # get logical vector which row does exceed threshold
    if (include_lowest) {
        f <- pr >= prevalence
    } else {
        f <- pr > prevalence
    }
    # get it back into order of x
    m <- match(rownames(x),names(f))
    taxa <- f[m]
    # Gets indices of most prevalent taxa
    indices <- which(taxa)
    # revert the order based on f
    m <- match(names(f),names(indices))
    m <- m[!is.na(m)]
    indices <- indices[m]
    # 
    indices
}

.get_prevalent_taxa <- function(x, rank = NULL, ...){
    if(is(x,"SummarizedExperiment")){
        x <- .agg_for_prevalence(x, rank = rank, ...)
    }
    indices <- .get_prevalent_indices(x, ...)
    # If named input return named output
    if( !is.null(rownames(x)) ){
        # Gets the names
        taxa <- rownames(x)[indices]
    } else {
        # Otherwise indices are returned
        taxa <- unname(indices)
    }
    unique(taxa)
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
    function(x, rank = NULL, prevalence = 50/100, 
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
           function(x, ...)
               standardGeneric("getRareTaxa"))

.get_rare_indices <- function(x, ...){
    indices <- .get_prevalent_indices(x = x, ...)
    # reverse the indices
    indices_x <- seq_len(nrow(x))
    f <- !(indices_x %in% indices)
    indices_new <- indices_x[f]
    indices_new
}
    
.get_rare_taxa <- function(x, rank = NULL, ...){
    if(is(x,"SummarizedExperiment")){
        x <- .agg_for_prevalence(x, rank = rank, ...)
    }
    indices <- .get_rare_indices(x, ...)
    #
    if( !is.null(rownames(x)) ){
        # Gets the names
        taxa <- rownames(x)[indices]
    } else {
        # Otherwise indices are returned
        taxa <- indices
    }
    unique(taxa)
}

#' @rdname getPrevalence
#' @export
setMethod("getRareTaxa", signature = c(x = "ANY"),
    function(x, prevalence = 50/100, include_lowest = FALSE, ...){
        .get_rare_taxa(x, rank = NULL, prevalence = prevalence,
                       include_lowest = include_lowest, ...)
    }
)

#' @rdname getPrevalence
#' @export
setMethod("getRareTaxa", signature = c(x = "SummarizedExperiment"),
    function(x, rank = NULL, prevalence = 50/100, 
             include_lowest = FALSE, ...){
        .get_rare_taxa(x, rank = rank, prevalence = prevalence,
                       include_lowest = include_lowest, ...)
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
        x <- .agg_for_prevalence(x, rank = rank, ...)
        prevalent_indices <- .get_prevalent_indices(x, ...)
        x[prevalent_indices, ]
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
        x <- .agg_for_prevalence(x, rank = rank, ...)
        rare_indices <- .get_rare_indices(x, ...)
        x[rare_indices, ]
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

