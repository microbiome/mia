#' decontam functions
#'
#' The \code{decontam} functions \code{isContaminant} and
#' \code{isNotContaminant} are made available for
#' \code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#' objects.
#'
#' @param seqtab,x
#'   a \code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#'
#' @param abund_values A single character value for selecting the
#'   \code{\link[SummarizedExperiment:SummarizedExperiment-class]{assay}}
#'   to use.
#'
#' @param name A name for the column of the colData in which the contaminant
#'   information should be stored.
#'
#' @param concentration \code{NULL} or a single \code{character} value. Defining
#'   a column with numeric values from the \code{colData} to use as
#'   concentration information. (default: \code{concentration = NULL})
#'
#' @param control \code{NULL} or a single \code{character} value. Defining a
#'   column with logical values from the \code{colData} to define control and
#'   non-control samples. (default: \code{control = NULL})
#'
#' @param batch \code{NULL} or a single \code{character} value. Defining a
#'   column with values interpretable as a factor from the \code{colData} to use
#'   as batch information. (default: \code{batch = NULL})
#'
#' @param ... arguments passed onto
#'   \code{\link[decontam:isContaminant]{decontam:isContaminant}} or
#'   \code{\link[decontam:isNotContaminant]{decontam:isNotContaminant}}.
#'   Currently these are \code{method} and \code{batch.combine}.
#'
#' @param threshold numeric scalar. See
#'   \code{\link[decontam:isContaminant]{decontam:isContaminant}} or
#'   \code{\link[decontam:isNotContaminant]{decontam:isNotContaminant}}
#'
#' @param normalize,detailed logical scalar. See
#'   \code{\link[decontam:isContaminant]{decontam:isContaminant}} or
#'   \code{\link[decontam:isNotContaminant]{decontam:isNotContaminant}}
#'
#' @param ... 
#' \itemize{
#'   \item{for \code{isContaminant}/ \code{isNotContaminant}: }{arguments
#'     passed on to \code{\link[decontam:isContaminant]{decontam:isContaminant}}
#'     or \code{\link[decontam:isNotContaminant]{decontam:isNotContaminant}}}
#'   \item{for \code{addContaminantQC}/\code{addNotContaminantQC}: }{arguments
#'     passed on to \code{isContaminant}/ \code{isNotContaminant}}
#' }
#'
#' @return for \code{isContaminant}/ \code{isNotContaminant} a \code{DataFrame}
#'   or for \code{addContaminantQC}/\code{addNotContaminantQC} a modified object
#'   of \code{class(x)}
#'
#' @name isContaminant
#'
#' @importFrom decontam isContaminant isNotContaminant
#'
#' @seealso
#' \code{\link[decontam:isContaminant]{decontam:isContaminant}},
#' \code{\link[decontam:isNotContaminant]{decontam:isNotContaminant}}
#'
#' @examples
#' data(esophagus)
#' # setup of some mock data
#' colData(esophagus)$concentration <- c(1,2,3)
#' colData(esophagus)$control <- c(FALSE,FALSE,TRUE)
#'
#' isContaminant(esophagus,
#'               method = "frequency",
#'               concentration = "concentration")
#' esophagus <- addContaminantQC(esophagus,
#'                               method = "frequency",
#'                               concentration = "concentration")
#' colData(esophagus)
#'
#' isNotContaminant(esophagus, control = "control")
#' esophagus <- addNotContaminantQC(esophagus, control = "control")
#' colData(esophagus)
NULL

#' @rdname isContaminant
#' @export
setMethod("isContaminant", signature = c(seqtab = "SummarizedExperiment"),
    function(seqtab,
             abund_values = "counts",
             name = "isContaminant",
             concentration = NULL,
             control = NULL,
             batch = NULL,
             threshold = 0.1,
             normalize = TRUE,
             detailed = TRUE,
             ...){
        # input check
        .check_assay_present(abund_values, seqtab)
        if(!.is_a_string(name)){
            stop("'name' must be single character value.",call. = FALSE)
        }
        if(!is.numeric(threshold) || length(threshold) != 1L){
            stop("'threshold' must be single numeric value.", call. = FALSE)
        }
        if(!.is_a_bool(normalize)){
            stop("'normalize' must be TRUE or FALSE.", call. = FALSE)
        }
        if(!.is_a_bool(detailed)){
            stop("'detailed' must be TRUE or FALSE.", call. = FALSE)
        }
        #
        if(!is.null(concentration)){
            concentration <- retrieveCellInfo(seqtab, by = concentration,
                                              search = "colData")$value
            if(!is.numeric(concentration)){
                stop("'concentration' must define a column of colData() ",
                     "containing numeric values.",
                     call. = FALSE)
            }
        }
        if(!is.null(control)){
            control <- retrieveCellInfo(seqtab, by = control,
                                        search = "colData")$value
            if(!is.logical(control)){
                stop("'control' must define a column of colData() ",
                     "containing logical values.",
                     call. = FALSE)
            }
        }
        if(!is.null(batch)){
            batch <- retrieveCellInfo(seqtab, by = batch,
                                      search = "colData")$value
            batch <- factor(batch, sort(unique(batch)))
        }
        mat <- assay(seqtab,abund_values)
        contaminant <- isContaminant(t(mat),
                                     conc = concentration,
                                     neg = control,
                                     batch = batch,
                                     threshold = threshold,
                                     normalize = normalize,
                                     detailed =  detailed,
                                     ...)
        if(is.data.frame(contaminant)){
            contaminant <- DataFrame(contaminant)
        }
        contaminant
    }
)

#' @rdname isContaminant
#' @export
setMethod("isNotContaminant", signature = c(seqtab = "SummarizedExperiment"),
    function(seqtab,
             abund_values = "counts",
             name = "isNotContaminant",
             control = NULL,
             threshold = 0.5,
             normalize = TRUE,
             detailed = FALSE,
             ...){
        # input check
        .check_assay_present(abund_values, seqtab)
        if(!.is_a_string(name)){
            stop("'name' must be single character value.",call. = FALSE)
        }
        if(!is.numeric(threshold) || length(threshold) != 1L){
            stop("'threshold' must be single numeric value.", call. = FALSE)
        }
        if(!.is_a_bool(normalize)){
            stop("'normalize' must be TRUE or FALSE.", call. = FALSE)
        }
        if(!.is_a_bool(detailed)){
            stop("'detailed' must be TRUE or FALSE.", call. = FALSE)
        }
        #
        if(!is.null(control)){
            control <- retrieveCellInfo(seqtab, by = control,
                                        search = "colData")$value
            if(!is.logical(control)){
                stop("'control' must define a column of colData() ",
                     "containing logical values.",
                     call. = FALSE)
            }
        }
        mat <- assay(seqtab,abund_values)
        not_contaminant <- isNotContaminant(t(mat),
                                            neg = control,
                                            threshold = threshold,
                                            normalize = normalize,
                                            detailed =  detailed,
                                            ...)
        if(is.data.frame(not_contaminant)){
            not_contaminant <- DataFrame(not_contaminant)
        }
        not_contaminant
    }
)

#' @rdname isContaminant
#' @export
setGeneric("addContaminantQC", signature = c("x"),
           function(x, name = "isContaminant", ...)
               standardGeneric("addContaminantQC"))

#' @rdname isContaminant
#' @export
setMethod("addContaminantQC", signature = c("SummarizedExperiment"),
    function(x, name = "isContaminant", ...){
        contaminant <- isContaminant(x, ...)
        rowData(x)[[name]] <- contaminant
        x
    }
)

#' @rdname isContaminant
#' @export
setGeneric("addNotContaminantQC", signature = c("x"),
           function(x, name = "isNotContaminant", ...)
               standardGeneric("addNotContaminantQC"))

#' @rdname isContaminant
#' @export
setMethod("addNotContaminantQC", signature = c("SummarizedExperiment"),
    function(x, name = "isNotContaminant", ...){
        not_contaminant <- isNotContaminant(x, ...)
        rowData(x)[[name]] <- not_contaminant
        x
    }
)
