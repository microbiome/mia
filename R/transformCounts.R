#' Transform Counts
#'
#' Variety of transformations for abundance data, stored in \code{assay}.
#' See details for options.
#'
#' @param x A
#'   \code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#'    object.
#'
#' @param assay_name A single character value for selecting the
#'   \code{\link[SummarizedExperiment:SummarizedExperiment-class]{assay}} to be
#'   transformed.
#'
#' @param abund_values a single \code{character} value for specifying which
#'   assay to use for calculation.
#'   (Please use \code{assay_name} instead. At some point \code{abund_values}
#'   will be disabled.)
#'   
#' @param method A single character value for selecting the transformation
#'   method.
#'
#' @param name A single character value specifying the name of transformed
#'   abundance table.
#'
#' @param pseudocount NULL or numeric value deciding whether pseudocount is
#'   added. The numeric value specifies the value of pseudocount.
#'
#' @param threshold A numeric value for setting threshold for pa transformation.
#'   By default it is 0. (Only used for \code{method = "pa"})
#'
#' @param ... additional arguments
#'
#' @details
#'
#' These functions provide a variety of options for transforming abundance data.
#' The transformed data is calculated and stored in a new \code{assay}.
#'
#' Available wrapper functions include:
#'
#' \itemize{
#'
#' \item{\code{transformSamples}} sample-wise (column-wise) transformation.
#' Alias for \code{transformCounts}.
#' 
#' \item{\code{transformFeatures}} feature-wise (row-wise) transformation.
#'
#' \item{\code{ZTransform}} Shortcut for Z-transformation.
#' 
#' \item{\code{relAbundanceCounts}} Shortcut for fetching relative abundance table.
#'
#' }
#' 
#' Altogether, \code{transformCounts} or \code{transformSamples} and \code{transformFeatures}
#' apply transformations to the abundance table (assay). The available transformation methods include:
#'
#' \itemize{
#'
#' \item{'clr'}{ Centered log ratio (clr) transformation aims to remove
#' compositionality effect; it is also used to
#' skewness and to center data.  
#' 
#' If the data contains zeros, pseudocount (commonly the smallest 
#' positive value of the data) must be added since clr is a logarithmic
#' transformation that only allows positive values.
#' (See e.g. Gloor et al. 2017.)
#'
#' \deqn{clr = log_{10}\frac{x{g(x)}} = log_{10}x - log_{10}\mu}{%
#' clr = log10(x/g(x)) = log10 x - log10 \mu}
#' where \eqn{x} is a single value, g(x) is geometric mean of
#' sample-wise values, and \eqn{\mu} is an arithmetic mean of 
#' sample-wise values.}
#'
#' \item{'rclr'}{ rclr or robust clr is similar to regular clr. Problem of regular
#' clr is that logarithmic transformations lead to undefined values when zeros
#' are present in the data. In rclr, values are divided by geometric mean
#' of observed taxa and zero values are not taken into account. Zero values will
#' stay as zeroes. Because of high-dimensionality of data, rclr's geometric mean of 
#' observed taxa is a good approximation to the true geometric mean.
#' (For details, see Martino et al. 2019.).
#'
#' }
#' 
#' \item{'hellinger'}{ Hellinger transformation can be used to reduce the impact of
#' extreme data points. It can be utilize for clustering or ordination analysis.
#' (See e.g. Legendre & Gallagher 2001.)
#'
#' \deqn{hellinger = \sqrt{\frac{x}{x_{tot}}}}{%
#' hellinger = sqrt(x/x_tot)}
#' where \eqn{x} is a single value and \eqn{x_{tot}} is the sum of
#' all values}
#' 
#' \item{'log10'}{ log10 transformation can be used for reducing the skewness of the data.
#'
#' \deqn{log10 = \log_{10} x}{%
#' log10 = log10(x)}
#' where \eqn{x} is a single value of data.}
#' 
#' \item{'pa'}{ Transforms table to presence/absence table. All abundances higher
#' than \eqn{\epsilon} are transformed to 1 (present), otherwise 0 (absent). 
#' By default, threshold is 0.}
#' 
#' \item{'rank'}{ Rank returns ranks of taxa. For each sample, the least abundant 
#' taxa get lower value and more abundant taxa bigger value. The implementation is 
#' based on the colRanks function with ties.method="first".}
#' 
#' \item {'relabundance'}{ Transforms abundances to relative. Generally, all microbiome
#' data are compositional. That is, e.g., because all measuring instruments have their capacity limits.
#' To make results comparable with other results, values must be relative. (See e.g. Gloor et al. 2017.)
#'
#' \deqn{relabundance = \frac{x}{x_{tot}}}{%
#' relabundance = x/x_tot}
#' where \eqn{x} is a single value and \eqn{x_{tot}} is the sum of
#' all values.}
#'
#' \item{'z'}{ Z-transformation, Z score transformation, or Z-standardization normalizes
#' the data by shifting (to mean \eqn{\mu}) and scaling (to standard deviation \eqn{\sigma}).
#' Z-transformation can be done with function \code{ZTransform}. It is done per rows (features / taxa),
#' unlike most other transformations. This is often preceded by log10p or clr transformation.
#' In other words, single value is standardized with respect of feature's values.
#'
#' \deqn{z = \frac{x - \mu}{\sigma}}{%
#' z = (x - µ)/σ}
#' where \eqn{x} is a single value, \eqn{\mu} is the mean of the feature, and
#' \eqn{\sigma} is the standard deviation of the feature.}
#'
#' }
#'
#'
#' @references
#' Gloor GB, Macklaim JM, Pawlowsky-Glahn V & Egozcue JJ (2017)
#' Microbiome Datasets Are Compositional: And This Is Not Optional.
#' Frontiers in Microbiology 8: 2224. doi: 10.3389/fmicb.2017.02224
#'
#' Legendre P & Gallagher ED (2001)
#' Ecologically meaningful transformations for ordination of species data.
#' Oecologia 129: 271-280.
#' 
#' Martino C, Morton JT, Marotz CA, Thompson LR, Tripathi A, Knight R & Zengler K
#' (2019) A Novel Sparse Compositional Technique Reveals Microbial Perturbations.
#' mSystems 4: 1. doi: 10.1128/mSystems.00016-19
#'
#'
#' @return
#' \code{transformCounts}, \code{transformSamples}, \code{transformFeatures}, 
#' \code{relAbundanceCounts}, and \code{ZTransform} return \code{x} with additional, 
#' transformed abundance table named \code{name} in the \code{\link{assay}}.
#'
#' @name transformCounts
#' @export
#'
#' @author Leo Lahti and Tuomas Borman. Contact: \url{microbiome.github.io}
#'
#' @examples
#' data(esophagus)
#' x <- esophagus
#'
#' # By specifying 'method', it is possible to apply different transformations, 
#' # e.g. compositional transformation.
#' x <- transformSamples(x, method = "relabundance")
#' 
#' # The target of transformation can be specified with "assay_name"
#' # Pseudocount can be added by specifying 'pseudocount'.
#' 
#' # Get pseudocount; here smallest positive value
#' mat <- assay(x, "relabundance") 
#' pseudonumber <- min(mat[mat>0])
#' # Perform CLR
#' x <- transformSamples(x, assay_name = "relabundance", method = "clr", 
#'                       pseudocount = pseudonumber
#'                       )
#'                       
#' head(assay(x, "clr"))
#' 
#' # Name of the stored table can be specified.
#' x <- transformSamples(x, method="hellinger", name="test")
#' head(assay(x, "test"))
#'
#' # pa returns presence absence table. With 'threshold', it is possible to set the
#' # threshold to a desired level. By default, it is 0.
#' x <- transformSamples(x, method = "pa", threshold = 35)
#' head(assay(x, "pa"))
#' 
#' # rank returns ranks of taxa. It is calculated column-wise, i.e., per sample
#' # and using the ties.method="first" from the colRanks function
#' x <- transformSamples(x, method = "rank")
#' head(assay(x, "rank"))
#' 
#' # transformCounts is an alias for transformSamples
#' x <- transformCounts(x, method = "relabundance", name = "test2")
#' head(assay(x, "test2"))
#'
#' # In order to use other ranking variants, modify the chosen assay directly:
#' assay(x, "rank_average", withDimnames = FALSE) <- colRanks(assay(x, "counts"), 
#'                                                            ties.method="average", 
#'                                                            preserveShape = TRUE)  
#'                                                            
#' # If you want to do the transformation for features, you can do that by using
#' x <- transformFeatures(x, method = "log10", name = "log10_features", pseudocount = 1)
#' head(assay(x, "log10_features"))
#'
#' # Z-transform can be done for features by using shortcut function
#' x <- ZTransform(x)
#' head(assay(x, "z"))
#' 
#' # For visualization purposes it is sometimes done by applying CLR for samples,
#' # followed by Z transform for taxa
#' x <- ZTransform(transformCounts(x, method="clr", assay_name = "counts", pseudocount = 1))
#'
#' # Relative abundances can be also calculated with the dedicated
#' # relAbundanceCounts function.
#' x <- relAbundanceCounts(x)
#' head(assay(x, "relabundance"))
NULL

#' @rdname transformCounts
#' @aliases transformCounts
#' @export
setGeneric("transformSamples", signature = c("x"),
           function(x,
                    assay_name = abund_values, abund_values = "counts",
                    method = c("alr", "clr", "rclr", "hellinger", "log10", "pa", 
                               "rank", "relabundance"),
                    name = method,
                    pseudocount = NULL,
                    threshold = 0,
                    ...
                    )
                    standardGeneric("transformSamples"))

#' @rdname transformCounts
#' @aliases transformCounts
#' @export
setMethod("transformSamples", signature = c(x = "SummarizedExperiment"),
    function(x,
            assay_name = abund_values, abund_values = "counts",
            method = c("alr", "clr", "rclr", "hellinger", "log10", "pa", 
                       "rank", "relabundance"),
            name = method,
            pseudocount = NULL,
            threshold = 0,
            ...
            ){
        # Input check
        # Check assay_name
        .check_assay_present(assay_name, x)
        # Check name
        if(!.is_non_empty_string(name) ||
           name == assay_name){
            stop("'name' must be a non-empty single character value and be ",
                 "different from `assay_name`.",
                 call. = FALSE)
        }
        # Check method
        # If method is not single string, user has not specified transform method,
        # or has given e.g. a vector
        if(!.is_non_empty_string(method)){
          stop("'method' must be a non-empty single character value.",
               call. = FALSE)
        }
        method <- match.arg(method)

        # Gets the abundance table
        assay <- assay(x, assay_name)
        # Calls help function that does the transformation
        transformed_table <- .apply_transformation(assay, method, pseudocount, threshold, ...)
        # Assign transformed table to assays
        assay(x, name, withDimnames=FALSE) <- transformed_table
        x
    }
)

#' @rdname transformCounts
#' @aliases transformSamples
#' @export
setGeneric("transformCounts", signature = c("x"),
           function(x,
                    assay_name = abund_values, abund_values = "counts",
                    method = c("alr", "clr", "rclr", "hellinger", "log10", "pa", 
                               "rank", "relabundance"),
                    name = method,
                    pseudocount = NULL,
                    threshold = 0,
                    ...)
               standardGeneric("transformCounts"))

#' @rdname transformCounts
#' @aliases transformSamples
#' @export
setMethod("transformCounts", signature = c(x = "SummarizedExperiment"),
    function(x,
             assay_name = abund_values, abund_values = "counts",
             method = c("alr", "clr", "rclr", "hellinger", "log10", "pa", 
                        "rank", "relabundance"),
             name = method,
             pseudocount = NULL,
             threshold = 0,
             ...){
        transformSamples(x, 
                       assay_name = assay_name,
                       method = method,
                       name = name,
                       pseudocount = pseudocount,
                       threshold = threshold, ...)
    }
)

###############################transformFeatures################################

#' @rdname transformCounts
#' @export
setGeneric("transformFeatures", signature = c("x"),
           function(x,
                    assay_name = abund_values, abund_values = "counts",
                    method = c("log10", "pa", "z"),
                    name = method,
                    pseudocount = NULL,
                    threshold = 0)
               standardGeneric("transformFeatures"))

#' @rdname transformCounts
#' @export
setMethod("transformFeatures", signature = c(x = "SummarizedExperiment"),
    function(x,
             assay_name = abund_values, abund_values = "counts",
             method = c("log10", "pa", "z"),
             name = method,
             pseudocount = NULL,
             threshold = 0){
        # Input check
        # Check assay_name
        .check_assay_present(assay_name, x)
        # Check name
        if(!.is_non_empty_string(name) ||
           name == assay_name){
            stop("'name' must be a non-empty single character value and be ",
                 "different from `assay_name`.",
                 call. = FALSE)
        }
        # Check method
        # If method is not single string, user has not specified transform method,
        # or has given e.g. a vector
        if(!.is_non_empty_string(method)){
          stop("'method' must be a non-empty single character value.",
               call. = FALSE)
        }
        method <- match.arg(method)
      
        # Gets the abundance table, and transposes it
        assay <- t(assay(x, assay_name))
        # Calls help function that does the transformation
        transformed_table <- .apply_transformation(assay, method, pseudocount, threshold)
        # Transposes transformed table to right orientation
        transformed_table <- t(transformed_table)
        # Assign transformed table to assays
        assay(x, name, withDimnames=FALSE) <- transformed_table
        x
  }
)
##################################Z-TRANSFORM###################################

#' @rdname transformCounts
#' @export
setGeneric("ZTransform", signature = c("x"),
           function(x, ...)
             standardGeneric("ZTransform"))

#' @rdname transformCounts
#' @importFrom SummarizedExperiment assay assay<-
#' @export
setMethod("ZTransform", signature = c(x = "SummarizedExperiment"),
          function(x, ...){
            transformFeatures(x, method = "z", ...)
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
setMethod("relAbundanceCounts",signature = c(x = "SummarizedExperiment"),
    function(x, ...){
        transformSamples(x, method = "relabundance", ...)
    }
)

###########################HELP FUNCTIONS####################################

# Help function for transformSamples and transformFeatures, takes abundance table
# as input and returns transformed table.
.apply_transformation <- function(assay, method, pseudocount, threshold, ...){
    # Input check
    # Check pseudocount
    if( !( is.null(pseudocount) || 
           (length(pseudocount) == 1L && is.numeric(pseudocount)) ) ){
        stop("'pseudocount' must be NULL or a single numeric value.",
             call. = FALSE)
    } 
    # Check threshold
    if(!is.numeric(threshold)){
        stop("'threshold' must be a numeric value, and it can be used ",
             "only with transformation method 'pa'.",
             call. = FALSE)
    }
    # Input check end
    
    # apply pseudocount, if it is numeric
    if( is.numeric(pseudocount) ){
        assay <- .apply_pseudocount(assay, pseudocount)
    }
    # Get transformed table
    transformed_table <-
        .get_transformed_table(assay = assay,
                               method = method,
                               threshold = threshold,
                               ...)
    return(transformed_table)
}

# Chooses which transformation function is applied
.get_transformed_table <- function(assay, method, threshold, ...){
    # Function is selected based on the "method" variable
    FUN <- switch(method,
                  relabundance = .calc_rel_abund,
                  log10 = .calc_log10,
                  pa = .calc_pa,
                  hellinger = .calc_hellinger,
                  clr = .calc_clr,
                  rank = .calc_rank,
                  z = .calc_ztransform,
                  rclr = .calc_clr,
                  alr = .calc_clr
                  )

    # Does the function call, arguments are "assay" abundance table and "pseudocount"
    do.call(FUN,
            list(mat = assay,
                 threshold = threshold,
                 method = method,
                 ...))
}

#' @importFrom DelayedMatrixStats colSums2
.calc_rel_abund <- function(mat, ...){
    mat <- sweep(mat, 2, colSums2(mat, na.rm = TRUE), "/")
    return(mat)
}

.calc_log10 <- function(mat, ...){
    # If abundance table contains zeros, gives an error, because it is not possible
    # to calculate log from zeros. If there is no zeros, calculates log.
    if (any(mat <= 0, na.rm = TRUE)) {
        stop("Abundance table contains zero or negative values and ",
            "log10 transformation is being applied without pseudocount.\n",
            "Try to add pseudocount (default choice pseudocount = 1 for count ",
            "assay; or pseudocount = min(x[x>0]) for relabundance assay).",
            call. = FALSE)
    }
    mat <- log10(mat)
    return(mat)
}

.calc_pa <- function(mat, threshold, ...){
    # If value is over zero, gets value 1. If value is zero, gets value 0.
    mat <- (mat > threshold) - 0L
    return(mat)
}

.calc_hellinger <- function(mat, ...){
    # If there is negative values, gives an error.
    if (any(mat < 0, na.rm = TRUE)) {
        stop("Abundance table contains negative values and hellinger ",
            "transformation is being applied without (suitable) pseudocount.\n",
            "Try to add pseudocount (default choice pseudocount = 1 for count ",
            "assay; or pseudocount = min(x[x>0]) with relabundance assay).",
            call. = FALSE)
    }
    # Gets the relative abundance
    mat <- .calc_rel_abund(mat)
    # Takes square root
    mat <- sqrt(mat)
    return(mat)
}

#' @importFrom DelayedMatrixStats colSums2 colMeans2
.calc_clr <- function(mat, method, reference = 1, MARGIN = 1, reference_values = NA, ...){
    # Calculate colSums
    colsums <- colSums2(mat, na.rm = TRUE)
    # Check that they are equal; affects the result of CLR. CLR expects a fixed
    # constant
    if( round(max(colsums)-min(colsums), 3) != 0  ){
        warning("All the total abundances of samples do not sum-up to a fixed constant.\n",
                "Please consider to apply, e.g., relative transformation in prior to ",
                "CLR transformation.",
                call. = FALSE)
    }
    # Utilize function from vegan (check vegan documentation!!!!!!!!!!!) there ar erelative abundance....
    mat <- t(mat)
    
    if( method == "alr" ){
        # Check reference_values
        if( length(reference_values) != 1 ){
            stop("'reference_values' must be a single value specifying the values of ",
                 "the reference sample.",
                 call. = FALSE)
        }
        # Reference sample
        reference_name <- rownames(mat)[reference]
        # Get the order of samples
        col_index <- seq_len(nrow(mat))
        
        mat <- vegan::decostand(mat, method = method, reference = reference, MARGIN = MARGIN)
        # Reference sample as NAs
        reference_sample <- matrix(reference_values, nrow = 1, ncol = ncol(mat),  
                                   dimnames = list(reference_name, colnames(mat)) )
        # Mat index
        mat_index <- col_index[-reference]
        
        # Add reference sample
        mat <- rbind(mat, reference_sample)
        
        # Preserve the original order
        mat <- mat[ match(col_index, c(mat_index, reference)), ]
        
    } else{
        mat <- vegan::decostand(mat, method = method)
    }
    mat <- t(mat)
    return(mat)
}

#' @importFrom DelayedMatrixStats colRanks
.calc_rank <- function(mat, ...){
    # Converts, e.g., DelayedArrays to matrix
    mat <- as.matrix(mat)
    # For every sample, finds ranks of taxa.
    # Column-wise, NAs are kept as NAs, and ties get the minimum rank value.
    # Transpose ensures that dimensions of matrix are right.
    mat <- colRanks(mat, ties.method = "first", preserveShape = TRUE)
    return(mat)
}

#' @importFrom DelayedMatrixStats colMeans2 colSds
.calc_ztransform <- function(mat, ...){
    # Converts, e.g., DelayedArrays to matrix
    mat <- as.matrix(mat)
    # Z transformation column-wise
    # Centers the feature data. After that, divides with
    # the standard deviation of feature.
    cm <- colMeans2(mat, na.rm = TRUE)
    csd <- colSds(mat, na.rm = TRUE)
    # Transposes the table, otherwise calculation below would be done in wrong 
    # direction, e.g, cm and csd for samples would be applied to rows. 
    mat <- t(mat)
    mat <- (mat - cm)/csd
    # Transposes the table
    mat <- t(mat)
    return(mat)
}

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
