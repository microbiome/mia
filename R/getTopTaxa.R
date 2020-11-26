#' Accessing most relevant features
#'
#' To query a \code{SummarizedExperiment} for interesting features, several
#' functions are available.
#'
#' @param x A
#'  \code{\link[=SummarizedExperiment-class]{SummarizedExperiment}} object.
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
#' \code{getTopTaxa} extracts the most \code{top} abundant \dQuote{FeatureID}s
#' in a \code{\link[=SummarizedExperiment-class]{SummarizedExperiment}} object.
#'
#' @return
#' For \code{getTopTaxa}: A vector of the most \code{top} abundant
#' \dQuote{FeatureID}s
#'
#' @seealso
#' \code{\link[=getPrevalence]{getPrevalentTaxa}}
#'
#' @name getTopTaxa
#'
#' @author
#' Sudarshan A. Shetty
#'
#' @examples
#' data(GlobalPatterns)
#' top_taxa <- getTopTaxa(GlobalPatterns,
#'                        method="mean",
#'                        top=5,
#'                        abund_values="counts")
#' top_taxa
NULL

#' @rdname getTopTaxa
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

#' @rdname getTopTaxa
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
        .check_abund_values(abund_values, x)
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

