#' Getter / setter for relative abundance data
#'
#' \code{relabundance} is a getter/setter for relative abundance stored in the
#' assay slot \sQuote{relabundance} of a
#' \code{\link[TreeSummarizedExperiment:TreeSummarizedExperiment-class]{TreeSummarizedExperiment}}
#' object
#'
#' @param x a
#'   \code{\link[TreeSummarizedExperiment:TreeSummarizedExperiment-class]{TreeSummarizedExperiment}}
#'   object
#'
#' @param abund_values A single character value for selecting the
#'   \code{\link[SummarizedExperiment:SummarizedExperiment-class]{assay}}
#'   to use for calculating the relative abundance.
#' @param name A single character value for selecting the
#'   \code{\link[SummarizedExperiment:SummarizedExperiment-class]{assay}}
#'   to stor the calculated relative abundance.
#' @param value a \code{matrix} to store as the the \sQuote{relabundance} assay
#' @param ... optional arguments not used currently.
#'
#' @return
#' For \code{relabundance} the matrix stored with the name
#' \dQuote{relabundance}.
#' For \code{relAbundanceCounts} a modified \code{x} containing the relative
#' abundances as an assay defined by \code{name}.
#'
#'
#' @name relabundance
#'
#' @export
#'
#' @examples
#' data(GlobalPatterns)
#' GlobalPatterns <- relAbundanceCounts(GlobalPatterns)
#' relabundance(GlobalPatterns)
NULL

#' @rdname relabundance
setGeneric("relabundance", signature = c("x"),
           function(x, ...) standardGeneric("relabundance"))
#' @rdname relabundance
setGeneric("relabundance<-", signature = c("x"),
           function(x, value) standardGeneric("relabundance<-"))
#' @rdname relabundance
setGeneric("relAbundanceCounts", signature = c("x"),
           function(x, abund_values = "counts", name = "relabundance")
               standardGeneric("relAbundanceCounts"))

#' @rdname relabundance
#' @importFrom SummarizedExperiment assays
#' @export
setMethod("relabundance",signature = c(x = "SummarizedExperiment"),
    function(x){
        assays(x)[["relabundance"]]
    }
)

#' @rdname relabundance
#' @importFrom SummarizedExperiment assays<-
#' @export
setReplaceMethod("relabundance", signature = c(x = "SummarizedExperiment"),
    function(x, value){
        assays(x)[["relabundance"]] <- value
        x
    }
)

.calc_rel_abund <- function(mat){
    sweep(mat, 2, colSums(mat), "/")
}

#' @rdname relabundance
#' @importFrom SummarizedExperiment assay assay<-
#' @export
setMethod("relAbundanceCounts",signature = c(x = "SummarizedExperiment"),
    function(x, abund_values = "counts", name = "relabundance"){
        # input check

        if(!.is_non_empty_string(name)){
            stop("'name' must be a single non-empty character value.",
                 call. = FALSE)
        }
        .check_abund_values(abund_values, x)
        #
        assay(x, name) <- .calc_rel_abund(assay(x, abund_values))
        x
    }
)
