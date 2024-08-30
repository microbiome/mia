#' Estimate divergence
#'
#' Estimate divergence against a given reference sample.
#' 
#' @inheritParams addDissimilarity
#' 
#' @param x a \code{\link{SummarizedExperiment}} object.
#'   
#' @param assay_name Deprecated. Use \code{assay.type} instead.
#'
#' @param reference \code{Character scalar}. A column name from
#' \code{colData(x)} or either \code{"mean"} or \code{"median"}.
#' (Default: \code{"median"})
#'
#' @param ... optional arguments
#' 
#' @return \code{x} with additional \code{\link{colData}} named \code{*name*}
#' 
#' @details
#'
#' Microbiota divergence (heterogeneity / spread) within a given sample
#' set can be quantified by the average sample dissimilarity or beta
#' diversity with respect to a given reference sample.
#'
#' This measure is sensitive to sample size.
#' Subsampling or bootstrapping can be applied to equalize sample sizes
#' between comparisons.
#' 
#' @seealso
#' \code{\link[scater:plotColData]{plotColData}}
#' \itemize{
#'   \item{\code{\link[mia:estimateRichness]{estimateRichness}}}
#'   \item{\code{\link[mia:estimateEvenness]{estimateEvenness}}}
#'   \item{\code{\link[mia:estimateDominance]{estimateDominance}}}
#' }
#' 
#' @name addDivergence
#' @export
#' 
#' @examples
#' data(GlobalPatterns)
#' tse <- GlobalPatterns
#' 
#' # By default, reference is median of all samples. The name of column where results
#' # is "divergence" by default, but it can be specified. 
#' tse <- addDivergence(tse)
#' 
#' # The method that are used to calculate distance in divergence and 
#' # reference can be specified. Here, euclidean distance and dist function from 
#' # stats package are used. Reference is the first sample. It is recommended
#' # to add reference to colData.
#' tse[["reference"]] <- rep(colnames(tse)[[1]], ncol(tse))
#' tse <- addDivergence(
#'     tse, name = "divergence_first_sample", 
#'     reference = "reference",
#'     method = "euclidean")
#' 
#' # Reference can also be median or mean of all samples. 
#' # By default, divergence is calculated by using median. Here, mean is used.
#' tse <- addDivergence(tse, name = "divergence_average", reference = "mean")
#' 
#' # All three divergence results are stored in colData.
#' colData(tse)
#' 
NULL

#' @rdname addDivergence
#' @export
setGeneric(
    "addDivergence",signature = c("x"),
    function(x, name = "divergence", ...)
    standardGeneric("addDivergence"))

#' @rdname addDivergence
#' @export
setMethod("addDivergence", signature = c(x="SummarizedExperiment"),
    function(x, name = "divergence", ...){
        ################### Input check ###############
        # Check name
        if( !.is_a_string(name) ){
            stop("'name' must be a non-empty character value.",
                 call. = FALSE)
        }
        ################# Input check end #############
        # Calculate values
        res <- getDivergence(x, ...)
        # Add them to colData
        x <- .add_values_to_colData(x, list(res), name)
        return(x)
    }
)

#' @rdname addDivergence
#' @export
setGeneric("getDivergence", signature = c("x"),
    function(
        x, assay.type = assay_name, assay_name = "counts", reference = "median",
        method = "bray", ...)
    standardGeneric("getDivergence"))

#' @rdname addDivergence
#' @export
setMethod("getDivergence", signature = c(x="SummarizedExperiment"),
    function(
        x, assay.type = assay_name, assay_name = "counts", 
        reference = "median", method = "bray", ...){
        ################### Input check ###############
        # Check assay.type
        .check_assay_present(assay.type, x)
        # Check reference
        ref_type <- .get_reference_type(reference, x)
        if( is.null(ref_type) ){
            stop(
                "'reference' must be a column from colData or either 'mean' ",
                "or 'median'.", call. = FALSE)
        }
        # If there are no colnmaes, add them. They are not added to returned
        # vaues.
        if( is.null(colnames(x)) ){
            colnames(x) <- paste0("sample_", seq_len(ncol(x)))
        }
        ################# Input check end #############
        # Get assay and referencess
        mat <- .get_matrix_and_reference(x, assay.type, reference, ref_type)
        reference <- mat[[2]]
        mat <- mat[[1]]
        # Calculate sample-wise divergence
        res <- .calc_divergence(mat, reference, method, ...)
        # Get only values and ensure that their order is correct
        res <- res[match(colnames(x), res[["sample"]]), "value"]
        return(res)
        }
)
############################## HELP FUNCTIONS ##############################

# This function returns reference typr.
# reference must be a column from colData, or either "median" or "mean".
# We also support providing a numeric vector or single sample name, but
# those those are not recommended.
.get_reference_type <- function(reference, x){
    is_col <- .is_a_string(reference) && reference %in% colnames(colData(x)) &&
        all(x[[reference]] %in% colnames(x))
    is_mean_or_median <- .is_a_string(reference) && reference %in% c(
        "mean", "median")
    is_num_vector <- is.numeric(reference) && length(reference) == nrow(x)
    is_char_vector <- is.character(reference) && length(reference) == ncol(x) &&
        all(reference %in% colnames(x))
    is_sample <- .is_a_string(reference) && reference %in% colnames(x)
    #
    res <- NULL
    if( is_col ){
        res <- "colData_column"
    } else if(is_mean_or_median){
        res <- reference
    } else if( is_num_vector ){
        res <- "num_vector"
    } else if( is_char_vector ){
        res <- "char_vector"
    } else if( is_sample ){
        res <- "sample"
    }
    return(res)
}

# This function gets the abundance table along with reference information
.get_matrix_and_reference <- function(
        x, assay.type, reference, ref_type,
        ref.name = "temporal_reference_for_divergence"){
    # Get assay
    mat <- assay(x, assay.type)
    # If reference type is median or mean, calculate it
    if( ref_type %in% c("median", "mean") ){
        reference <- apply(mat, 1, ref_type)
    }
    # In case of numeric values, add them to matrix
    if( ref_type %in% c("num_vector", "median", "mean") ){
        reference <- matrix(reference)
        colnames(reference) <- ref.name
        mat <- cbind(mat, reference)
        reference <- ref.name
    }
    # In case of colData variable, get reference samples from there
    if( ref_type %in% c("colData_column") ){
        reference <- x[[reference]]
    }
    # If the reference is only one sample, replicate it to cover all samples
    if( .is_a_string(reference) ){
        reference <- rep(reference, ncol(mat))
    }
    # Return a list with matrix and reference samples for each sample
    res <- list(mat, reference)
    return(res)
}

# For each sample-pair, this function calculates dissimilarity.
#' @importFrom dplyr mutate
.calc_divergence <- function(mat, reference, method, ...){
    # Create sample-pair data.frame
    reference <- data.frame(sample = colnames(mat), reference = reference)
    # Calculate dissimilarity between samples
    mat <- t(mat)
    mat <- getDissimilarity(mat, method, ...)
    mat <- as.matrix(mat)
    # Assign values to reference table
    reference <- reference |>
        mutate(value = mat[cbind(sample, reference)])
    return(reference)
}
