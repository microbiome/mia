#' decontam functions
#'
#' The \code{decontam} functions \code{isContaminant} and
#' \code{isNotContaminant} are made available for
#' \code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#' objects.
#'
#' @inheritParams getDissimilarity
#' @inheritParams getDominant
#'
#' @param seqtab,x a
#' \code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#'
#' @param name \code{Character scalar}. A name for the column of the
#' \code{colData} where results will be stored.
#' (Default: \code{"contaminant"} or \code{"not_contaminant"})
#'
#' @param concentration \code{Character scalar} or \code{NULL}. Defining
#'   a column with numeric values from the \code{colData} to use as
#'   concentration information. (Default: \code{NULL})
#'
#' @param control \code{Character scalar} or \code{NULL}. Defining a
#'   column with logical values from the \code{colData} to define control and
#'   non-control samples. (Default: \code{NULL})
#'
#' @param batch \code{Character scalar} or \code{NULL}. Defining a
#'   column with values interpretable as a factor from the \code{colData} to use
#'   as batch information. (Default: \code{NULL})
#'
#' @param detailed \code{Logical scalar}. If \code{TRUE}, the return value is a
#'   data.frame containing diagnostic information on the contaminant decision.
#'   If FALSE, the return value is a logical vector containing the binary
#'   contaminant classifications. (Default: \code{TRUE})
#'
#' @param ... arguments passed onto
#'   \code{\link[decontam:isContaminant]{decontam:isContaminant}} or
#'   \code{\link[decontam:isNotContaminant]{decontam:isNotContaminant}}.
#'   Currently these are \code{method} and \code{batch.combine}.
#'
#' @param threshold  \code{Numeric scalar}.. See
#'   \code{\link[decontam:isContaminant]{decontam:isContaminant}} or
#'   \code{\link[decontam:isNotContaminant]{decontam:isNotContaminant}}
#'
#' @param normalize \code{Logical scalar}. See
#'   \code{\link[decontam:isContaminant]{decontam:isContaminant}} or
#'   \code{\link[decontam:isNotContaminant]{decontam:isNotContaminant}}
#'
#' @param ...
#' \itemize{
#'   \item for \code{isContaminant}/ \code{isNotContaminant}: arguments
#'     passed on to \code{\link[decontam:isContaminant]{decontam:isContaminant}}
#'     or \code{\link[decontam:isNotContaminant]{decontam:isNotContaminant}}
#'   \item for \code{addContaminantQC}/\code{addNotContaminantQC}: arguments
#'     passed on to \code{isContaminant}/ \code{isNotContaminant}
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
#' # setup of some mock data just for example
#' colData(esophagus)$concentration <- c(1, 2, 3)
#' colData(esophagus)$control <- c(FALSE, FALSE, TRUE)
#'
#' isContaminant(
#'     esophagus,
#'     method = "frequency",
#'     concentration = "concentration"
#'     )
#' esophagus <- addContaminantQC(
#'     esophagus,
#'     method = "frequency",
#'     concentration = "concentration"
#'     )
#' rowData(esophagus)
#'
#' isNotContaminant(esophagus, control = "control")
#' esophagus <- addNotContaminantQC(esophagus, control = "control")
#' rowData(esophagus)
#'
NULL

#' @rdname isContaminant
#' @export
setMethod("isContaminant", signature = c(seqtab = "SummarizedExperiment"),
    function(
        seqtab,
        assay.type = assay_name, assay_name = "counts",
        concentration = NULL,
        control = NULL,
        batch = NULL,
        threshold = 0.1,
        normalize = TRUE,
        detailed = TRUE,
        ...){
        # input check
        .check_assay_present(assay.type, seqtab)
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
        # Get data
        concentration <- .get_concentration(seqtab, concentration)
        control <- .get_control(seqtab, control)
        batch <- .get_batch(seqtab, batch)
        mat <- assay(seqtab,assay.type)
        args <- list(
            seqtab = t(mat),
            conc = concentration,
            neg = control,
            batch = batch,
            threshold = threshold,
            normalize = normalize,
            detailed =  detailed
        )
        # We do not pass this arguments as they are controlled already by
        # concentration and control, respectively.
        args <- c(args, list(...)[ !names(list(...)) %in% c("conc", "neg") ])
        # Run analysis
        contaminant <- do.call(isContaminant, args)
        contaminant <- .wrangle_contaminant_result(contaminant, args)
        return(contaminant)
    }
)

#' @rdname isContaminant
#' @export
setMethod("isNotContaminant", signature = c(seqtab = "SummarizedExperiment"),
    function(
        seqtab,
        assay.type = assay_name, assay_name = "counts",
        control = NULL,
        threshold = 0.5,
        normalize = TRUE,
        detailed = FALSE,
        ...){
        # input check
        .check_assay_present(assay.type, seqtab)
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
        # Get data
        control <- .get_control(seqtab, control)
        mat <- assay(seqtab,assay.type)
        args <- list(
            seqtab = t(mat),
            neg = control,
            threshold = threshold,
            normalize = normalize,
            detailed =  detailed
        )
        # Do not pass this parameter as it is already controlled by control.
        args <- c(args, list(...)[ !names(list(...)) %in% c("neg") ])
        # Run analysis
        not_contaminant <- do.call(isNotContaminant, args)
        not_contaminant <- .wrangle_contaminant_result(not_contaminant, args)
        return(not_contaminant)
    }
)

#' @rdname isContaminant
#' @export
setMethod("addContaminantQC", signature = c("SummarizedExperiment"),
    function(x, name = "contaminant", ...){
        if(!.is_a_string(name)){
            stop("'name' must be single character value.", call. = FALSE)
        }
        contaminant <- isContaminant(x, ...)
        x <- .add_decontam_res(x, contaminant, name, type = "decontam_")
        return(x)
    }
)

#' @rdname isContaminant
#' @export
setMethod("addNotContaminantQC", signature = c("SummarizedExperiment"),
    function(x, name = "not_contaminant", ...){
        if(!.is_a_string(name)){
            stop("'name' must be single character value.", call. = FALSE)
        }
        not_contaminant <- isNotContaminant(x, ...)
        x <- .add_decontam_res(x, not_contaminant, name, type = "decontam_not_")
        return(x)
    }
)

################################ HELP FUNCTIONS ################################

# This function retrieves concentration from colData.
.get_concentration <- function(x, concentration, ...){
    if( !(is.null(concentration) || .is_a_string(concentration)) ){
        stop("'concentration' must be NULL or a single character value.",
             call. = FALSE)
    }
    if(!is.null(concentration)){
        concentration <- retrieveCellInfo(
            x, by = concentration, search = "colData")$value
        if( !is.numeric(concentration) ){
            stop("'concentration' must define a column of colData() ",
                "containing numeric values.", call. = FALSE)
        }
    }
    return(concentration)
}

# This function retrieves batch info from colData
.get_batch <- function(x, batch, ...){
    if(!is.null(batch)){
        batch <- retrieveCellInfo(
            x, by = batch, search = "colData")$value
        batch <- factor(batch, sort(unique(batch)))
    }
    return(batch)
}

# This function retrieves control samples from colData
.get_control <- function(x, control, ...){
    if( !(is.null(control) || .is_a_string(control)) ){
        stop("'control' must be NULL or a single character value.",
            call. = FALSE)
    }
    if(!is.null(control)){
        control <- retrieveCellInfo(
            x, by = control, search = "colData")$value
        if( !is.logical(control) ){
            stop("'control' must define a column of colData() ",
                "containing logical values.", call. = FALSE)
        }
    }
    return(control)
}

# This function wrangles the result of analysis to final, returned format
.wrangle_contaminant_result <- function(res, args){
    # The result can be a data.frame if user wanted detailed info. In that
    # situation, we convert it to DF
    if(is.data.frame(res)){
        res <- DataFrame(res)
    }
    # Add analysis arguments to attributes
    attr(res, "metadata") <- args[ !names(args) %in% c("seqtab")]
    return(res)
}

# This function adds decontam results to TreeSE
.add_decontam_res <- function(x, res, name, type){
    # save metadata
    add_metadata <- attr(res, "metadata")
    attr(res, "metadata") <- NULL
    names(add_metadata) <- paste0(type, names(add_metadata))
    # Add to rowData. The result can be either a vector or DF.
    if( is(res, "DataFrame") ){
        values <- as.list(res)
        name <- paste0(name, "_", names(values))
        x <- .add_values_to_colData(
            x, values = unname(values), name = name, MARGIN = 1L)
    } else{
        rowData(x)[[name]] <- res
    }
    metadata(x) <- c(metadata(x), add_metadata)
    return(x)
}