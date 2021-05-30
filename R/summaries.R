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
#' @param abund_values a \code{character} value to select an
#'   \code{\link[SummarizedExperiment:SummarizedExperiment-class]{assayNames}}
#'
#' @details
#' The \code{getTopTaxa} extracts the most \code{top} abundant \dQuote{FeatureID}s
#' in a \code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#' object.
#'
#' @return
#' For \code{getTopTaxa}: A vector of the most \code{top} abundant
#' \dQuote{FeatureID}s
#'
#' @seealso
#' \code{\link[=getPrevalence]{getPrevalentTaxa}}
#'
#' @name summaries-basic
#'
#' @author
#' Leo Lahti, Tuomas Borman and Sudarshan A. Shetty
#'
#' @examples
#' data("GlobalPatterns")
#' top_taxa <- getTopTaxa(GlobalPatterns,
#'                        method = "mean",
#'                        top = 5,
#'                        abund_values = "counts")
#' top_taxa
#'
#' # Gets the overview of dominant taxa
#' dominant_taxa <- summarizeDominantTaxa(GlobalPatterns, rank = "Family", name = "dominant_family")
#' dominant_taxa
#'
#' # With group, it is possible to group observations based on specified groups
#' # Gets the overview of dominant taxa
#' dominant_taxa <- summarizeDominantTaxa(GlobalPatterns, group = "SampleType")
#' dominant_taxa
#'
#' #
#' summarizeSE(GlobalPatterns)
#'
NULL

#' @rdname summaries-basic
#'
#' @export
setGeneric("getTopTaxa", signature = "x",
           function(x, top= 5L, method = c("mean","sum","median"),
                    abund_values = "counts")
               standardGeneric("getTopTaxa"))

.check_max_taxa <- function(x, top, abund_values){
    if(!is.numeric(top) || as.integer(top) != top){
        top("'top' must be integer value", call. = FALSE)
    }
    if(top > nrow(assay(x,abund_values))){
        stop("'top' must be <= nrow(x)", call. = FALSE)
    }
}

#' @rdname summaries-basic
#'
#' @importFrom DelayedMatrixStats rowSums2 rowMeans2 rowMedians
#' @importFrom utils head
#'
#' @export
setMethod("getTopTaxa", signature = c(x = "SummarizedExperiment"),
    function(x, top = 5L, method = c("mean","sum","median","prevalence"),
             abund_values = "counts"){
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
        head(names(taxs), n = top)
    }
)

#' @rdname summaries-basic
#'
#' @param rank A single character defining a taxonomic rank. Must be a value of
#'   the output of \code{taxonomicRanks()}.
#'
#' @param group With group, it is possible to group the observations in an
#'   overview. Must be a one of the column names of \code{colData}.
#'
#' @param name A name for the column of tibble table that includes taxa.
#'
#' @details
#' \code{summarizeDominantTaxa} returns information about most dominant
#' taxa in a tibble. Information includes their absolute and relative abundances in whole
#' data set.
#' \itemize{
#'         \item{\code{rank}}{ with \code{rank} parameter, it is possible to
#'         agglomerate taxa based on taxonomic ranks. E.g. if 'family' rank is
#'         used, all abundances of same family is added together, and those
#'         families are returned determining the threshold for coverage index.
#'         By default, \code{threshold} is 0.9.}
#'         \item{\code{group}}{ with \code{group} parameter, it is possible to
#'         group observations of returned overview based on samples' features.
#'         E.g., if samples contain information about patients' health status,
#'         it is possible to group observations, e.g. to 'healthy' and 'sick',
#'         and get the most dominant taxa of different health status.}
#' }
#'
#'
#' @return
#' \code{summarizeDominantTaxa} returns an overview in a tibble. It contains dominant taxa
#' in a column named \code{*name*} and its abundance in the data set.
#'
#' @export
setGeneric("summarizeDominantTaxa",signature = c("x"),
           function(x,
                    abund_values = "counts",
                    rank = NULL,
                    group = NULL,
                    name = "dominant_taxa")
               standardGeneric("summarizeDominantTaxa"))

#' @rdname summaries-basic
#' @export
setMethod("summarizeDominantTaxa", signature = c(x = "SummarizedExperiment"),
    function(x,
             abund_values = "counts",
             rank = NULL,
             group = NULL,
             name = "dominant_taxa"){
        # Input check
        # group check
        if(!is.null(group)){
            if(isFALSE(any(group %in% colnames(colData(x))))){
                stop("'group' variable must be in colnames(colData(x))")
                }
            }
        # Adds dominant taxas to colData
        tmp <- dominantTaxa(x, abund_values, rank, name)
        # Gets an overview
        overview <- .get_overview(tmp, group, name)
        overview
  }
)

################################ HELP FUNCTIONS ################################

#' @importFrom S4Vectors as.data.frame
#' @importFrom dplyr n desc tally group_by arrange mutate
.get_overview <- function(x, group, name){
    # Creates a tibble df that contains dominant taxa and number of times that
    # they present in samples and relative portion of samples where they
    # present.
    if (is.null(group)) {
        name <- sym(name)
        overview <- as.data.frame(colData(x)) %>%
            group_by(!!name) %>%
            tally() %>%
            mutate(
                rel.freq = round(100 * n / sum(n), 1),
                rel.freq.pct = paste0(round(100 * n / sum(n), 0), "%")
            ) %>%
            arrange(desc(n))
    } else {
        group <- sym(group)
        name <- sym(name)

        overview <- as.data.frame(colData(x)) %>%
            group_by(!!group, !!name) %>%
            tally() %>%
            mutate(
                rel.freq = round(100 * n / sum(n), 1),
                rel.freq.pct = paste0(round(100 * n / sum(n), 0), "%")
            ) %>%
            arrange(desc(n))
    }
    return(overview)
}

#' @rdname summaries-basic
#'
#' @param x A
#'  \code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}} object.
#'
#' @param abund_values a \code{character} value to select an
#'   \code{\link[SummarizedExperiment:SummarizedExperiment-class]{assayNames}}
#'   By default it expects count data.
#'
#' @param ... additional arguments not used.
#'
#' @details
#' The \code{summarizeSE} will return a summary of counts for all samples and
#' features in
#' \code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#' object.
#'
#' @return
#' For \code{summarizeSE}: A list with two \code{tibble}s
#'
#' @seealso
#' \code{\link[scuttle:perCellQCMetrics]{perCellQCMetrics}},
#' \code{\link[scuttle:perFeatureQCMetrics]{perFeatureQCMetrics}},
#' \code{\link[scuttle:addPerCellQC]{addPerCellQC}},
#' \code{\link[scuttle:addPerFeatureQC]{addPerFeatureQC}},
#' \code{\link[scuttle:quickPerCellQC]{quickPerCellQC}}
#'
#' @export
setGeneric("summarizeSE",
           signature = c("x"),
           function(x, ...)
             standardGeneric("summarizeSE")
)

#' @rdname summaries-basic
#' @export
setMethod("summarizeSE",
          signature = c(x = "SummarizedExperiment"),
          function(x, abund_values = "counts"){
            .check_rel_neg(x, abund_values)
            sample.summary <- .get_summary_col_data(x, abund_values)
            feature.summary <- .get_summary_rowl_data(x, abund_values)
            return(list("samples" = sample.summary, "features" = feature.summary))
          }
)

################################ HELP FUNCTIONS summarizeSE ####################

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
.get_summary_rowl_data <- function(x, abund_values){
  # should check and extract assay
  assay.x <- .get_assay(x, abund_values)
  summary_row_data <- tibble(total = nrow(assay.x),
                             singletons = .get_singletons(assay.x),
                             per_sample_avg = mean(colSums2(assay.x != 0)),
                             median_counts = median(colSums2(assay.x)),
                             mean_counts = mean(colSums2(assay.x)),
                             stdev_counts = sd(colSums2(assay.x)))
  return(summary_row_data)
}



# Get singletons from assay matrix
#' @importFrom DelayedMatrixStats rowSums2
#' @importFrom tibble tibble
.get_singletons<-function(x){
  length(rowSums2(x)[rowSums2(x) == 1])
}

# Check is relative or negative values
#' @importFrom DelayedMatrixStats colSums2
.check_rel_neg <- function(x, abund_values){
  assay.x <- .get_assay(x, abund_values)
    if(any(colSums2(assay.x) < 1) | any(colSums2(assay.x) < 0)){
    stop("There are samples that sum to 1 or less counts. ",
         "Try to supply raw counts",
         call. = FALSE)
  }
}
