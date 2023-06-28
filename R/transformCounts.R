#' Transform Counts
#'
#' Variety of transformations for abundance data, stored in \code{assay}.
#' See details for options.
#'
#' @param x A
#'   \code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#'    object.
#'
#' @param assay.type A single character value for selecting the
#'   \code{\link[SummarizedExperiment:SummarizedExperiment-class]{assay}} to be
#'   transformed.
#'
#' @param assay_name a single \code{character} value for specifying which
#'   assay to use for calculation.
#'   (Please use \code{assay.type} instead. At some point \code{assay_name}
#'   will be disabled.)
#'   
#' @param method A single character value for selecting the transformation
#'   method.
#' 
#' @param MARGIN A single character value for specifying whether the
#'   transformation is applied sample (column) or feature (row) wise.
#'   (By default: \code{MARGIN = "samples"})
#'
#' @param name A single character value specifying the name of transformed
#'   abundance table.
#' 
#' @param pseudocount NULL or numeric value deciding whether pseudocount is
#'   added. The numeric value specifies the value of pseudocount.
#'   Recommended default choices for counts and relative abundance assay
#'   \code{pseudocount = 1} and \code{pseudocount = min(assay[assay>0])}, respectively.
#'   (By default: \code{pseudocount = 0})
#'
#' @param ... additional arguments passed on to \code{vegan:decostand}:
#' \itemize{
#'   \item{\code{ref_vals}:} {A single value which will be used to fill 
#'   reference sample's column in returned assay when calculating alr. 
#'   (default: \code{ref_vals = NA})}
#' }
#' @details
#'
#' These \code{transformCount} function provides a variety of options for transforming abundance data.
#' The transformed data is calculated and stored in a new \code{assay}. The previously available
#' wrappers transformSamples, transformFeatures
#' ZTransform, and relAbundanceCounts have been deprecated.
#'
#' The \code{transformCounts} provides sample-wise (column-wise) or feature-wise
#' (row-wise) transformation to the abundance table
#' (assay) based on specified \code{MARGIN}.
#' 
#' The available transformation methods include:
#'
#' \itemize{
#' 
#' \item{'alr'}{ Additive log ratio (alr) transformation, please refer to 
#' \code{\link[vegan:decostand]{decostand}} for details.}
#' 
#' \item{'chi.square'}{ Chi square transformation, please refer to 
#' \code{\link[vegan:decostand]{decostand}} for details.}
#' 
#' \item{'clr'}{ Centered log ratio (clr) transformation, please refer to 
#' \code{\link[vegan:decostand]{decostand}} for details.}
#'
#' \item{'frequency'}{ Frequency transformation, please refer to 
#' \code{\link[vegan:decostand]{decostand}} for details.}
#' 
#' \item{'hellinger'}{ Hellinger transformation, please refer to 
#' \code{\link[vegan:decostand]{decostand}} for details.}
#' 
#' \item{'log'}{ Logarithmic transformation, please refer to 
#' \code{\link[vegan:decostand]{decostand}} for details.}
#' 
#' \item{'log10'}{ log10 transformation can be used for reducing the skewness
#' of the data.
#' \deqn{log10 = \log_{10} x}{%
#' log10 = log10(x)}
#' where \eqn{x} is a single value of data.}
#' 
#' \item{'log2'}{ log2 transformation can be used for reducing the skewness of
#' the data.
#' \deqn{log2 = \log_{2} x}{%
#' log2 = log2(x)}
#' where \eqn{x} is a single value of data.}
#' 
#' \item{'normalize'}{ Make margin sum of squares equal to one. Please refer to 
#' \code{\link[vegan:decostand]{decostand}} for details.}
#' 
#' \item{'pa'}{ Transforms table to presence/absence table. Please refer to 
#' \code{\link[vegan:decostand]{decostand}} for details.}
#'
#' \item{'rank'}{ Rank transformation, please refer to 
#' \code{\link[vegan:decostand]{decostand}} for details.}
#' 
#' \item{'rclr'}{ Robust clr transformation, please refer to 
#' \code{\link[vegan:decostand]{decostand}} for details.}
#' 
#' \item{'relabundance'}{ Relative transformation (alias for 'total'), please refer to 
#' \code{\link[vegan:decostand]{decostand}} for details.}
#' 
#' \item{'rrank'}{ Relative rank transformation, please refer to 
#' \code{\link[vegan:decostand]{decostand}} for details.}
#' 
#' \item{'standardize'}{ Scale 'x' to zero mean and unit variance (alias for
#' 'z'), please refer to \code{\link[vegan:decostand]{decostand}} for details.}
#' 
#' \item{'total'}{ Divide by margin total (alias for
#' 'relabundance'), please refer to 
#' \code{\link[vegan:decostand]{decostand}} for details.}
#' 
#' \item{'z'}{ Z transformation (alias for 'standardize'),
#' please refer to \code{\link[vegan:decostand]{decostand}} for details.}
#'
#' }
#'
#' @return
#' \code{transformCounts} returns the input object \code{x}, with a new 
#' transformed abundance table named \code{name} added in the \code{\link{assay}}.
#'
#' @name transformCounts
#' @export
#'
#' @author Leo Lahti and Tuomas Borman. Contact: \url{microbiome.github.io}
#'
#' @examples
#' data(esophagus, package="mia")
#' x <- esophagus
#'
#' # By specifying 'method', it is possible to apply different transformations, 
#' # e.g. compositional transformation.
#' x <- transformCounts(x, method = "relabundance")
#' 
#' # The target of transformation can be specified with "assay.type"
#' # Pseudocount can be added by specifying 'pseudocount'.
#' 
#' # Get pseudocount; here smallest positive value
#' mat <- assay(x, "relabundance") 
#' pseudonumber <- min(mat[mat>0])
#' # Perform CLR
#' x <- transformCounts(x, assay.type = "relabundance", method = "clr", 
#'                      pseudocount = pseudonumber
#'                      )
#'                       
#' head(assay(x, "clr"))
#' 
#' # With MARGIN, you can specify the if transformation is done for samples or
#' # for features. Here Z-transformation is done feature-wise.
#' x <- transformCounts(x, method = "z", MARGIN = "features")
#' head(assay(x, "z"))
#' 
#' # Name of the stored table can be specified.
#' x <- transformCounts(x, method="hellinger", name="test")
#' head(assay(x, "test"))
#'
#' # pa returns presence absence table.
#' x <- transformCounts(x, method = "pa")
#' head(assay(x, "pa"))
#' 
#' # rank returns ranks of taxa.
#' x <- transformCounts(x, method = "rank")
#' head(assay(x, "rank"))
#'
#' # In order to use other ranking variants, modify the chosen assay directly:
#' assay(x, "rank_average", withDimnames = FALSE) <- colRanks(assay(x, "counts"), 
#'                                                            ties.method="average", 
#'                                                            preserveShape = TRUE)  
#' 
NULL

#' @rdname transformCounts
#' @aliases transformCounts
#' @export
setGeneric("transformSamples", signature = c("x"),
           function(x,
                    assay.type = "counts", assay_name = NULL,
                    method = c("alr", "chi.square", "clr", "frequency", "hellinger",
                               "log", "log10", "log2", "normalize", "pa",
                               "rank", "rclr", "relabundance", "rrank",
                               "total"),
                    name = method,
                    ...
                    )
                    standardGeneric("transformSamples"))

#' @rdname transformCounts
#' @aliases transformCounts
#' @export
setMethod("transformSamples", signature = c(x = "SummarizedExperiment"),
    function(x,
            assay.type = "counts", assay_name = NULL,
            method =  c("alr", "chi.square", "clr", "frequency", "hellinger",
                        "log", "log10", "log2", "normalize", "pa",
                        "rank", "rclr", "relabundance", "rrank",
                        "total"),
            name = method,
            pseudocount = 0,
            ...
            ){
        .Deprecated("transformCounts")
        # Input check
        # Check method
        # If method is not single string, user has not specified transform method,
        # or has given e.g. a vector
        if(!.is_non_empty_string(method)){
            stop("'method' must be a non-empty single character value.",
                 call. = FALSE)
        }
        method <- match.arg(method, several.ok = FALSE)
        # Input check end

        # Call general transformation function with MARGIN specified
        x <- transformCounts(x = x, assay.type = assay.type,
                             method = method, MARGIN = "samples", name = name, ...)
        return(x)
    }
)

#' @rdname transformCounts
#' @aliases transformSamples
#' @export
setGeneric("transformCounts", signature = c("x"),
           function(x,
                    assay.type = "counts", assay_name = NULL,
                    method = c("alr", "chi.square", "clr", "frequency",
                               "hellinger", "log", "log10", "log2", "max",
                               "normalize", "pa", "range", "rank", "rclr",
                               "relabundance", "rrank", "standardize", "total",
                               "z"),
                    MARGIN = "samples",
                    name = method,
                    pseudocount = 0,		    
                    ...)
               standardGeneric("transformCounts"))

#' @rdname transformCounts
#' @aliases transformSamples
#' @export
setMethod("transformCounts", signature = c(x = "SummarizedExperiment"),
    function(x,
             assay.type = "counts", assay_name = NULL,
             method = c("alr", "chi.square", "clr", "frequency", "hellinger",
                        "log", "log10", "log2", "max", "normalize", "pa",
                        "range", "rank", "rclr", "relabundance", "rrank",
                        "standardize", "total", "z"),
             MARGIN = "samples",
             name = method,
             pseudocount = 0,
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
        if( !(length(MARGIN) == 1L && MARGIN %in%
              c("samples", "features", "columns", "col", "row")) ){
            stop("'MARGIN' must be 'samples' or 'features'.",
                 call. = FALSE)
        }
        # Check pseudocount
        if( !(is.numeric(pseudocount) && length(pseudocount) == 1 && pseudocount >= 0) ){
            stop("'pseudocount' must be a non-negative single numeric value.",
                 call. = FALSE)
        }
        # Input check end

        # Get the method and abundance table
        method <- match.arg(method)
        assay <- assay(x, assay.type)
    
        # Apply pseudocount, if it is not 0
        if( pseudocount != 0 ){
            assay <- .apply_pseudocount(assay, pseudocount)
        }
	
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

###############################transformFeatures################################

#' @rdname transformCounts
#' @export
setGeneric("transformFeatures", signature = c("x"),
           function(x,
                    assay.type = "counts", assay_name = NULL,
                    method = c("frequency", "log", "log10", "log2", "max",
                               "pa", "range", "standardize", "z"),
                    name = method,
                    pseudocount = 0,
                    ...)
               standardGeneric("transformFeatures"))

#' @rdname transformCounts
#' @export
setMethod("transformFeatures", signature = c(x = "SummarizedExperiment"),
    function(x,
             assay.type = "counts", assay_name = NULL,
             method = c("frequency", "log", "log10", "log2", "max",
                        "pa", "range", "standardize", "z"),
             name = method,
             pseudocount = 0,
             ...){

        .Deprecated("transformCounts")

        # Input check
        # Check method
        # If method is not single string, user has not specified transform method,
        # or has given e.g. a vector
        if(!.is_non_empty_string(method)){
          stop("'method' must be a non-empty single character value.",
               call. = FALSE)
        }
        method <- match.arg(method, several.ok = FALSE)
        # Input check end

        # Call general transformation function with MARGIN specified
        x <- transformCounts(x = x, assay.type = assay.type, method = method,
                             MARGIN = "features", name = name, ...)
        return(x)
  }
)
##################################Z-TRANSFORM###################################

#' @rdname transformCounts
#' @export
setGeneric("ZTransform", signature = c("x"),
           function(x, MARGIN = "features", ...)
             standardGeneric("ZTransform"))

#' @rdname transformCounts
#' @importFrom SummarizedExperiment assay assay<-
#' @export
setMethod("ZTransform", signature = c(x = "SummarizedExperiment"),
          function(x, ...){
	    .Deprecated("transformCounts")
            transformCounts(x, method = "z", MARGIN = "features", ...)
          }
)

###############################relAbundanceCounts###############################

#' @rdname transformCounts
#' @export
setGeneric("relAbundanceCounts", signature = c("x"),
            function(x, ...)
            standardGeneric("relAbundanceCounts"))

#' @rdname transformCounts
#' @importFrom SummarizedExperiment assay assay<-
#' @export
setMethod("relAbundanceCounts", signature = c(x = "SummarizedExperiment"),
    function(x, ...){
	.Deprecated("transformCounts")    
        transformCounts(x, method = "relabundance", MARGIN = "samples", ...)
    }
)

###########################HELP FUNCTIONS####################################
##############################.apply_transformation#############################
# Help function for transformSamples and transformFeatures, takes abundance table
# as input and returns transformed table. This function utilizes mia's
# transformation functions.
.apply_transformation <- function(assay, method, MARGIN, ...){

    # Transpose if MARGIN is row
    if( MARGIN %in% c("features", "row") ){
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
    if( MARGIN %in% c("features", "row") ){
        transformed_table <- t(transformed_table)
    }

    # Add method and margin to attributes
    attr(transformed_table, "mia") <- method
    attr(transformed_table, "parameters")$margin <- MARGIN
    return(transformed_table)
}

########################.apply_transformation_from_vegan########################
# Help function for transformSamples and transformFeatures, takes abundance
# table as input and returns transformed table. This function utilizes vegan's
# transformation functions.
.apply_transformation_from_vegan <- function(mat, method, MARGIN, ref_vals = NA, ...){
    # Input check
    # Check ref_vals
    if( length(ref_vals) != 1 ){
        stop("'ref_vals' must be a single value specifying the ",
             "values of the reference sample.",
             call. = FALSE)
    }
    # Input check end

    # Adjust MARGIN for vegan. It requires MARGIN in numeric format
    MARGIN <- ifelse(MARGIN %in% c("samples", "columns", "col", 2), 2, 1)
    # Adjust method if mia-specific alias was used
    method <- ifelse(method == "relabundance", "total", method)
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
            ref_vals = ref_vals)
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
    # If abundance table contains zeros, gives an error, because it is not
    # possible to calculate log from zeros. If there is no zeros, calculates log.
    if (any(mat <= 0, na.rm = TRUE)) {
        stop("Abundance table contains zero or negative values and ",
            method, " transformation is being applied without pseudocount.\n",
            "Try to add pseudocount (default choice pseudocount = 1 for count ",
            "assay; or pseudocount = min(x[x>0]) for relabundance assay).",
            call. = FALSE)
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
.adjust_alr_table <- function(mat, orig_dimnames, ref_vals){
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
    reference_sample <- matrix(ref_vals, nrow = nrow, ncol = ncol,  
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
.apply_pseudocount <- function(mat, pseudocount){
    # Give warning if pseudocount should not be added
    if( all(mat>0) ){
        warning("The abundance table contains only positive values. ",
                "A pseudocount is not encouraged to apply.", call. = FALSE)
    }
    # Add pseudocount.
    mat <- mat + pseudocount
    return(mat)
}
