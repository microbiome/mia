#' Transform assay
#'
#' Variety of transformations for abundance data, stored in \code{assay}.
#' See details for options.
#'
#' @inheritParams getDominant
#' @inheritParams getDissimilarity
#'   
#' @param method \code{Character scalar}. Specifies the transformation
#'   method.
#' 
#' @param MARGIN \code{Character scalar}. Determines whether the
#'   transformation is applied sample (column) or feature (row) wise.
#'   (Default: \code{"samples"})
#' 
#' @param pseudocount \code{Logical scalar} or \code{numeric scalar}. 
#'   When \code{TRUE}, automatically adds half of the minimum positive 
#'   value of \code{assay.type} (missing values ignored by default: 
#'   \code{na.rm = TRUE}).
#'   When FALSE, does not add any pseudocount (pseudocount = 0).
#'   Alternatively, a user-specified numeric value can be added as pseudocount.
#'   (Default: \code{FALSE}).
#' 
#' @param name \code{Character scalar}. A name for the column of the 
#'   \code{colData} where results will be stored. (Default: \code{"method"})
#'
#' @param ... additional arguments passed on to \code{vegan:decostand}:
#' \itemize{
#'   \item \code{reference}: \code{Character scalar}. use to
#'   to fill reference sample's column in returned assay when calculating alr. 
#'   (Default: \code{NA})
#'   \item \code{ref_vals} Deprecated. Use \code{reference} instead.
#' }
#' @details
#'
#' \code{transformAssay} function provides a variety of options for
#' transforming abundance data. The transformed data is calculated and stored
#' in a new \code{assay}.
#'
#' The \code{transformAssay} provides sample-wise (column-wise) or feature-wise
#' (row-wise) transformation to the abundance table
#' (assay) based on specified \code{MARGIN}.
#' 
#' The available transformation methods include:
#'
#' \itemize{
#' 
#' \item 'alr', 'chi.square', 'clr', 'frequency', 'hellinger', 'log', 'normalize', 'pa', 'rank', 'rclr'
#' 'relabundance', 'rrank', 'standardize', 'total': please refer to 
#' \code{\link[vegan:decostand]{decostand}} for details.
#' 
#' \item 'log10': log10 transformation can be used for reducing the skewness
#' of the data.
#' \deqn{log10 = \log_{10} x}{%
#' log10 = log10(x)}
#' where \eqn{x} is a single value of data.
#' 
#' \item 'log2': log2 transformation can be used for reducing the skewness of
#' the data.
#' \deqn{log2 = \log_{2} x}{%
#' log2 = log2(x)}
#' where \eqn{x} is a single value of data.
#'
#' }
#'
#' @return
#' \code{transformAssay} returns the input object \code{x}, with a new 
#' transformed abundance table named \code{name} added in the \code{\link{assay}}.
#'
#' @name transformAssay
#' @export
#'
#' @author Leo Lahti and Tuomas Borman. Contact: \url{microbiome.github.io}
#'
#' @examples
#' data(esophagus, package="mia")
#' tse <- esophagus
#'
#' # By specifying 'method', it is possible to apply different transformations, 
#' # e.g. compositional transformation.
#' tse <- transformAssay(tse, method = "relabundance")
#' 
#' # The target of transformation can be specified with "assay.type"
#' # Pseudocount can be added by specifying 'pseudocount'.
#' 
#' # Perform CLR with smallest positive value as pseudocount
#' tse <- transformAssay(tse, assay.type = "relabundance", method = "clr", 
#'                      pseudocount = TRUE
#'                      )
#'                       
#' head(assay(tse, "clr"))
#' 
#' # With MARGIN, you can specify the if transformation is done for samples or
#' # for features. Here Z-transformation is done feature-wise.
#' tse <- transformAssay(tse, method = "standardize", MARGIN = "features")
#' head(assay(tse, "standardize"))
#' 
#' # Name of the stored table can be specified.
#' tse <- transformAssay(tse, method="hellinger", name="test")
#' head(assay(tse, "test"))
#'
#' # pa returns presence absence table.
#' tse <- transformAssay(tse, method = "pa")
#' head(assay(tse, "pa"))
#' 
#' # rank returns ranks of taxa.
#' tse <- transformAssay(tse, method = "rank")
#' head(assay(tse, "rank"))
#'
#' # In order to use other ranking variants, modify the chosen assay directly:
#' assay(tse, "rank_average", withDimnames = FALSE) <- colRanks(
#'     assay(tse, "counts"), ties.method = "average", preserveShape = TRUE)  
#' 
NULL

#' @rdname transformAssay
#' @export
setGeneric("transformAssay", signature = c("x"),
           function(x,
                    assay.type = "counts", assay_name = NULL,
                    method = c("alr", "chi.square", "clr", "frequency",
                               "hellinger", "log", "log10", "log2", "max",
                               "normalize", "pa", "range", "rank", "rclr",
                               "relabundance", "rrank", "standardize", "total",
                               "z"),
                    MARGIN = "samples",
                    name = method,
                    pseudocount = FALSE,		    
                    ...)
               standardGeneric("transformAssay"))

#' @rdname transformAssay
#' @export
setMethod("transformAssay", signature = c(x = "SummarizedExperiment"),
    function(x,
             assay.type = "counts", assay_name = NULL,
             method = c("alr", "chi.square", "clr", "frequency", "hellinger",
                        "log", "log10", "log2", "max", "normalize", "pa",
                        "range", "rank", "rclr", "relabundance", "rrank",
                        "standardize", "total", "z"),
             MARGIN = "samples",
             name = method,
             pseudocount = FALSE,
             ...){
        # Input check

        if (!is.null(assay_name)) {
            .Deprecated(old="assay_name", new="assay.type", "Now assay_name is deprecated. Use assay.type instead.")
	    assay.type <- assay_name
        }

        # Check assay.type
        .check_assay_present(assay.type, x)

        # Check name
        if(!.is_non_empty_string(name) ||
           name == assay.type){
            stop("'name' must be a non-empty single character value and be ",
                 "different from `assay.type`.",
                 call. = FALSE)
        }
        # Check method
        # If method is not single string, user has not specified transform method,
        # or has given e.g. a vector
        if(!.is_non_empty_string(method)){
            stop("'method' must be a non-empty single character value.",
                 call. = FALSE)
        }
        method <- match.arg(method, several.ok = FALSE)
        # Check that MARGIN is 1 or 2
        MARGIN <- .check_MARGIN(MARGIN)
        # Check pseudocount
        if( !.is_a_bool(pseudocount) && !(is.numeric(pseudocount) && length(pseudocount) == 1 && pseudocount >= 0) ){
            stop("'pseudocount' must be TRUE, FALSE or a number equal to or greater than 0.",
                 call. = FALSE)
        }
        # Input check end

        # Get the method and abundance table
        method <- match.arg(method)
        assay <- assay(x, assay.type)
        
        # Apply pseudocount, if it is not 0
        assay <- .apply_pseudocount(assay, pseudocount, ...)
        # Store pseudocount value and set attr equal to NULL
        pseudocount <- attr(assay, "pseudocount")
        attr(assay, "pseudocount") <- NULL
        
        # Calls help function that does the transformation
        # Help function is different for mia and vegan transformations
        if( method %in% c("log10", "log2") ){
            transformed_table <- .apply_transformation(
                assay, method, MARGIN, ...)
        } else{
            transformed_table <- .apply_transformation_from_vegan(
                assay, method, MARGIN, ...)
        }
        
        # Add pseudocount info to transformed table
        attr(transformed_table, "parameters")$pseudocount <- pseudocount
        
        # Assign transformed table to assays
        assay(x, name, withDimnames=FALSE) <- transformed_table
        x
    }
)

###########################HELP FUNCTIONS####################################
##############################.apply_transformation#############################
# Help function for transformAssay, takes abundance table
# as input and returns transformed table. This function utilizes mia's
# transformation functions.
.apply_transformation <- function(assay, method, MARGIN, ...){

    # Transpose if MARGIN is row
    if( MARGIN == 1L ){
        assay <- t(assay)
    }

    # Function is selected based on the "method" variable
    FUN <- switch(method,
                  log10 = .calc_log,
                   log2 = .calc_log,
    )

    # Get transformed table
    transformed_table <- do.call(
        FUN, list(mat = assay, method = method, ...) )

    # Transpose back to normal if MARGIN is row
    if( MARGIN == 1L ){
        transformed_table <- t(transformed_table)
    }

    # Add method and margin to attributes
    attr(transformed_table, "mia") <- method
    attr(transformed_table, "parameters")$margin <- MARGIN
    return(transformed_table)
}

########################.apply_transformation_from_vegan########################
# Help function for transformAssay, takes abundance
# table as input and returns transformed table. This function utilizes vegan's
# transformation functions.
.apply_transformation_from_vegan <- function(mat, method, MARGIN, reference = ref_vals,
                                            ref_vals = NA, ...){
    # Input check
    # Check reference
    if( length(reference) != 1 ){
        stop("'reference' must be a single value specifying the ",
             "values of the reference sample.",
             call. = FALSE)
    }
    # Input check end
    
    # Adjust method if mia-specific alias was used
    method <- ifelse(method == "relabundance", "total", method)
    
    if (method == "z") {
        .Deprecated(old="z", new="standardize")
    }
    method <- ifelse(method == "z", "standardize", method)
    
    # If method is ALR, vegan drops one column/sample, because it is used
    # as a reference. To work with TreeSE, reference sample must be added back.
    # Get the original order of samples/features
    orig_dimnames <- dimnames(mat)

    # Call vegan::decostand and apply transformation
    transformed_table <- vegan::decostand(mat, method = method, MARGIN = MARGIN, ...)
    
    # Add reference sample back if ALR
    if( method %in% c("alr") ){
        transformed_table <- .adjust_alr_table(
            mat = transformed_table, orig_dimnames = orig_dimnames,
            reference = reference)
    }
    # If table is transposed (like in chi.square), transpose back
    if(identical(rownames(transformed_table), colnames(mat)) &&
       identical(colnames(transformed_table), rownames(mat)) &&
       ncol(transformed_table) != ncol(mat) &&
       nrow(transformed_table != nrow(mat))){
        transformed_table <- t(transformed_table)
    }
    return(transformed_table)
}

####################################.calc_log###################################
# This function applies log transformation to abundance table.
.calc_log <- function(mat, method, ...){
    # If abundance table contains zeros or negative values, gives an error, because
    # it is not possible to calculate log from zeros. Otherwise, calculates log.
    if ( any(mat < 0, na.rm = TRUE) ){
        stop("The assay contains negative values and ", method,
             " transformation is being applied without pseudocount.",
            "`pseudocount` must be specified manually.", call. = FALSE)
    } else if ( any(mat == 0, na.rm = TRUE) ){
        stop("The assay contains zeroes and ", method,
             " transformation is being applied without pseudocount.",
             "`pseudocount` must be set to TRUE.", call. = FALSE)
    }
    # Calculate log2 or log10 abundances
    if(method == "log2"){
        mat <- log2(mat)
    } else{
        mat <- log10(mat)
    }
    # Add parameter to attributes
    attr(mat, "parameters") <- list()
    return(mat)
}

#################################.calc_rel_abund################################
# This function is for other functions to use internally.
.calc_rel_abund <- function(mat, ...){
    mat <- .apply_transformation_from_vegan(
        mat, method = "relabundance", MARGIN = 2)
    return(mat)
}

###############################.adjust_alr_table################################
# vegan::decostand returns ALR transformed abundance table without reference
# sample. Because in TreeSE all assays must have same row and column names,
# the reference sample is assigned back to transformed abundance table.
.adjust_alr_table <- function(mat, orig_dimnames, reference){
    # Store attributes
    attributes <- attributes(mat)
    # Get original and current sample/feature names and dimensions of reference
    # based on what dimensions misses names
    MARGIN <- ifelse(length(orig_dimnames[[1]]) == nrow(mat), 2, 1)
    orig_names <- if(MARGIN == 1){orig_dimnames[[1]]} else {orig_dimnames[[2]]}
    current_names <- if(MARGIN == 1){rownames(mat)} else {colnames(mat)}
    nrow <- ifelse(MARGIN == 1, 1, nrow(mat))
    ncol <- ifelse(MARGIN == 1, ncol(mat), 1)
    # Get the name of reference sample/feature and names of other dimension
    reference_name <- setdiff(orig_names, current_names)
    var_names <- if(MARGIN == 1){colnames(mat)} else {rownames(mat)}
    if(MARGIN == 1){
        ref_dimnames <- list(reference_name, var_names)
    } else {
        ref_dimnames <- list(var_names, reference_name)
        }
    # Reference sample as NAs or with symbols that are specified by user
    reference_sample <- matrix(reference, nrow = nrow, ncol = ncol,  
                               dimnames = ref_dimnames)
    # Add reference sample/feature
    if(MARGIN == 1){
        mat <- rbind(mat, reference_sample)
        # Preserve the original order
        mat <- mat[orig_names, ]
    } else {
        mat <- cbind(mat, reference_sample)
        # Preserve the original order
        mat <- mat[, orig_names]
    }
    # Add those attributes that were related to calculation
    attributes(mat) <- c(attributes(mat),
                         attributes[ !names(attributes) %in%
                                         c("dim", "dimnames") ])
    return(mat)
}

###############################.apply_pseudocount###############################
# This function applies pseudocount to abundance table.
.apply_pseudocount <- function(mat, pseudocount, na.rm = TRUE, ...){
    if( .is_a_bool(pseudocount) ){
        # If pseudocount TRUE and some NAs, a warning is issued
        if ( pseudocount && any(is.na(mat)) ){
            warning("The assay contains missing values (NAs). These will be 
                ignored in pseudocount calculation.", call. = FALSE)
        }    
        # If pseudocount TRUE but some negative values, numerical pseudocount needed
        if ( pseudocount && any(mat < 0, na.rm = TRUE) ){
            stop("The assay contains negative values. ",
                 "'pseudocount' must be specified manually.", call. = FALSE)
        }
        # If pseudocount TRUE, set it to  half of non-zero minimum value
        # else set it to zero.
        # Get min value
        value <- min(mat[mat > 0], na.rm = na.rm)
        value <- value/2
        pseudocount <- ifelse(pseudocount, value, 0)
        # Report pseudocount if positive value
        if ( pseudocount > 0 ){
            message("A pseudocount of ", pseudocount, " was applied.")
        }
    }
    # Give warning if pseudocount should not be added
    # Case 1: only positive values
    if( pseudocount != 0 && all(mat > 0, na.rm = TRUE) ){
        warning(
            "The assay contains only positive values. ",
            "Applying a pseudocount may be unnecessary.", call. = FALSE)
    }
    # Case 2: some negative values
    if( pseudocount != 0 && any(mat < 0, na.rm = TRUE) ){
        warning(
            "The assay contains some negative values. ",
            "Applying a pseudocount may produce meaningless data.",
            call. = FALSE)
    }
    # Add pseudocount
    mat <- mat + pseudocount
    # Set attr equal to pseudocount
    attr(mat, "pseudocount") <- pseudocount
    return(mat)
}
