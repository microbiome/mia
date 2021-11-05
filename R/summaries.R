#' Summarizing microbiome data
#'
#' To query a \code{SummarizedExperiment} for interesting features, several
#' functions are available.
#'
#' @param x A
#'  \code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}} object.
#'
#' @param top Numeric value, how many top taxa to return. Default return top
#'   five taxa.
#'
#' @param method Specify the method to determine top taxa. Either sum, mean,
#'   median or prevalence. Default is 'mean'.
#'
#' @param abund_values A \code{character} value to select an
#'   \code{\link[SummarizedExperiment:SummarizedExperiment-class]{assayNames}}
#'
#' @param ... Additional arguments:
#'    \itemize{
#'        \item{\code{sort}}{A single boolean value for selecting 
#'        whether to sort taxa in alphabetical order or not. Enabled in functions
#'        \code{getUniqueTaxa}, and \code{getTopTaxa}.
#'        (By default: \code{sort = FALSE})}
#'        \item{\code{na.rm}}{A single boolean value for selecting 
#'        whether to remove missing values or not. Enabled in functions
#'        \code{getUniqueTaxa}, and \code{getTopTaxa}.
#'        (By default: \code{sort = FALSE})}
#'    }
#'    
#' @details
#' The \code{getTopTaxa} extracts the most \code{top} abundant \dQuote{FeatureID}s
#' in a \code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#' object.
#'
#' The \code{getUniqueTaxa} is a basic function to access different taxa at a
#' particular taxonomic rank.
#'
#' @return
#' The \code{getTopTaxa} returns a vector of the most \code{top} abundant
#' \dQuote{FeatureID}s
#'
#' @seealso
#' \code{\link[=getPrevalence]{getPrevalentTaxa}}
#'
#' @name summaries
#'
#' @author
#' Leo Lahti, Tuomas Borman and Sudarshan A. Shetty
#'
#' @examples
#' data(GlobalPatterns)
#' top_taxa <- getTopTaxa(GlobalPatterns,
#'                        method = "mean",
#'                        top = 5,
#'                        abund_values = "counts")
#' top_taxa
#'
#' # Gets the overview of dominant taxa
#' dominant_taxa <- countDominantTaxa(GlobalPatterns,
#'                                    rank = "Genus")
#' dominant_taxa
#'
#' # With group, it is possible to group observations based on specified groups
#' # Gets the overview of dominant taxa
#' dominant_taxa <- countDominantTaxa(GlobalPatterns,
#'                                    rank = "Genus",
#'                                    group = "SampleType",
#'                                    na.rm= TRUE)
#'
#' dominant_taxa
#'
#' # Get an overview of sample and taxa counts
#' summary(GlobalPatterns)
#'
#' # Get unique taxa at a particular taxonomic rank
#' # sort = TRUE means that output is sorted in alphabetical order. 
#' # sort can also be used in function getTopTaxa
#' getUniqueTaxa(GlobalPatterns, "Phylum", sort = TRUE)
#'
NULL

#' @rdname summaries
#' @export
setGeneric("getTopTaxa", signature = "x",
           function(x, top= 5L, method = c("mean","sum","median"),
                    abund_values = "counts", ...)
               standardGeneric("getTopTaxa"))

.check_max_taxa <- function(x, top, abund_values){
    if(!is.numeric(top) || as.integer(top) != top){
        stop("'top' must be integer value", call. = FALSE)
    }
    if(top > nrow(assay(x,abund_values))){
        stop("'top' must be <= nrow(x)", call. = FALSE)
    }
}

#' @rdname summaries
#'
#' @importFrom DelayedMatrixStats rowSums2 rowMeans2 rowMedians
#' @importFrom utils head
#'
#' @export
setMethod("getTopTaxa", signature = c(x = "SummarizedExperiment"),
    function(x, top = 5L, method = c("mean","sum","median","prevalence"),
             abund_values = "counts", ...){
        # input check
        method <- match.arg(method, c("mean","sum","median","prevalence"))
        # check max taxa
        .check_max_taxa(x, top, abund_values)
        # check assay
        .check_assay_present(abund_values, x)
        #
        if(method == "prevalence"){
            taxs <- getPrevalence(assay(x, abund_values), sort = TRUE,
                                  include_lowest = TRUE)
        } else {
            taxs <- switch(method,
                           mean = rowMeans2(assay(x, abund_values)),
                           sum = rowSums2(assay(x, abund_values)),
                           median = rowMedians(assay(x, abund_values)))
            names(taxs) <- rownames(assay(x))
            taxs <- sort(taxs,decreasing = TRUE)
        }
        names <- head(names(taxs), n = top)
        # Remove NAs and sort if specified
        names <- .remove_NAs_and_sort(names, ... )
        return(names)
    }
)

#' @rdname summaries
#'
#' @param rank A single character defining a taxonomic rank. Must be a value of
#' the output of \code{taxonomicRanks()}.
#'
#' @return
#' The \code{getUniqueTaxa} returns a vector of unique taxa present at a
#' particular rank
#' 
#' @export
setGeneric("getUniqueTaxa",
           signature = c("x"),
           function(x, ...)
               standardGeneric("getUniqueTaxa"))

#' @rdname summaries
#' @export
setMethod("getUniqueTaxa", signature = c(x = "SummarizedExperiment"),
    function(x, rank = NULL, ...){
        .check_taxonomic_rank(rank, x)
        names <- unique(rowData(x)[,rank])
        # Remove NAs and sort if specified
        names <- .remove_NAs_and_sort(names, ... )
        return(names)
    }
)


#' @rdname summaries
#'
#' @param group With group, it is possible to group the observations in an
#'   overview. Must be one of the column names of \code{colData}.
#'
#' @param ... Additional arguments passed on to \code{agglomerateByRank()} when
#'   \code{rank} is specified for \code{countDominantTaxa}.
#'
#' @details
#' \code{countDominantTaxa} returns information about most dominant
#' taxa in a tibble. Information includes their absolute and relative
#' abundances in whole data set.
#'
#'
#' @return
#' The \code{countDominantTaxa} returns an overview in a tibble. It contains dominant taxa
#' in a column named \code{*name*} and its abundance in the data set.
#' 
#' @export
setGeneric("countDominantTaxa",signature = c("x"),
           function(x, group = NULL,  ...)
               standardGeneric("countDominantTaxa"))

#' @rdname summaries
#' @export
setMethod("countDominantTaxa", signature = c(x = "SummarizedExperiment"),
    function(x, group = NULL, ...){
        # Input check
        # group check
        if(!is.null(group)){
            if(isFALSE(any(group %in% colnames(colData(x))))){
                stop("'group' variable must be in colnames(colData(x))",
                     call. = FALSE)
            }
        }
        # Adds dominant taxas to colData
        dominant_taxa <- perSampleDominantTaxa(x, ...)
        data <- colData(x)
        data$dominant_taxa <- dominant_taxa
        # Gets an overview
        .tally_col_data(data, group, name = "dominant_taxa")
    }
)

################################ HELP FUNCTIONS ################################

#' @importFrom S4Vectors as.data.frame
#' @importFrom dplyr n desc tally group_by arrange mutate
.tally_col_data <- function(data, group, name){
    # Creates a tibble that contains number of times that a column of "name"
    # is present in samples and relative portion of samples where they
    # present.
    if (is.null(group)) {
        name <- sym(name)
        data <- as.data.frame(data) %>%
            group_by(!!name)
    } else {
        group <- sym(group)
        name <- sym(name)
        data <- as.data.frame(data) %>%
            group_by(!!group, !!name)
    }
    tallied_data <- data %>%
        tally() %>%
        mutate(
            rel.freq = round(100 * n / sum(n), 1),
            rel.freq.pct = paste0(round(100 * n / sum(n), 0), "%")
        ) %>%
        arrange(desc(n))
    return(tallied_data)
}

#' @rdname summaries
#'
#' @param object A
#'  \code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}} object.
#'
#' @param abund_values a \code{character} value to select an
#'   \code{\link[SummarizedExperiment:SummarizedExperiment-class]{assayNames}}
#'   By default it expects count data.
#'
#' @details
#' The \code{summary} will return a summary of counts for all samples and
#' features in
#' \code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#' object.
#'
#' @return
#' The \code{summary} returns a list with two \code{tibble}s
#'
#' @seealso
#' \code{\link[scuttle:perCellQCMetrics]{perCellQCMetrics}},
#' \code{\link[scuttle:perFeatureQCMetrics]{perFeatureQCMetrics}},
#' \code{\link[scuttle:addPerCellQC]{addPerCellQC}},
#' \code{\link[scuttle:addPerFeatureQC]{addPerFeatureQC}},
#' \code{\link[scuttle:quickPerCellQC]{quickPerCellQC}}
#'
#' @export
setMethod("summary", signature = c(object = "SummarizedExperiment"),
    function(object, abund_values = "counts"){
        # check if counts
        .check_fraction_or_negative_values(object, abund_values)
        sample.summary <- .get_summary_col_data(object, abund_values)
        feature.summary <- .get_summary_row_data(object, abund_values)
        return(list("samples" = sample.summary, "features" = feature.summary))
    }
)

################################ HELP FUNCTIONS summary ####################

#' @importFrom DelayedMatrixStats colSums2
#' @importFrom stats sd median
#' @importFrom tibble tibble
.get_summary_col_data <- function(x, abund_values){
    # should check and extract assay
    assay.x <- .get_assay(x, abund_values)
    summary_col_data <- tibble(total_counts = sum(colSums2(assay.x)),
                               min_counts = min(colSums2(assay.x)),
                               max_counts = max(colSums2(assay.x)),
                               median_counts = median(colSums2(assay.x)),
                               mean_counts = mean(colSums2(assay.x)),
                               stdev_counts = sd(colSums2(assay.x)))
    return(summary_col_data)
}

#' @importFrom DelayedMatrixStats colSums2
#' @importFrom tibble tibble
.get_summary_row_data <- function(x, abund_values){
    # Should check and extract assay
    # Internal from splitByRanks
    assay.x <- .get_assay(x, abund_values)
    summary_row_data <- tibble(total = nrow(assay.x),
                               singletons = .get_singletons(assay.x),
                               per_sample_avg = mean(colSums2(assay.x != 0)))
    return(summary_row_data)
}

# Get singletons from assay matrix
#' @importFrom DelayedMatrixStats rowSums2
#' @importFrom tibble tibble
.get_singletons<-function(x){
    length(rowSums2(x)[rowSums2(x) == 1])
}

# Check if values in assay are fractions or or negative values
#' @importFrom DelayedMatrixStats colSums2
.check_fraction_or_negative_values <- function(x, abund_values){
    assay.x <- .get_assay(x, abund_values)
    if(any(colSums2(assay.x) < 1) | any(colSums2(assay.x) < 0)){
        stop("There are samples that sum to 1 or less counts. ",
             "Try to supply raw counts",
             call. = FALSE)
    }
}

# Remove NAs and order in alphabetical order
.remove_NAs_and_sort <- function(names, sort = FALSE, na.rm = FALSE, ...){
    # Check sort
    if( !.is_a_bool(sort) ){
        stop("'sort' must be a boolean value.", call. = FALSE)
    }
    # Check na.rm
    if( !.is_a_bool(na.rm) ){
        stop("'na.rm' must be a boolean value.", call. = FALSE)
    }
    # Remove NAs if specified
    if( na.rm ){
        names <- names[ !is.na(names) ]
    }
    # Sort in alphabetical order if sort is TRUE
    if( sort ){
        names <- sort(names, na.last = TRUE)
    }
    return(names)
}
