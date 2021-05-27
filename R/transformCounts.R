#' Transform Counts
#'
#' These functions provide a variety of options for transforming abundance data.
#' By using these functions, transformed table is calculated and stored in \code{assay}. 
#' \code{transformSamples} does the transformation sample-wise, i.e., column-wise. 
#' It is alias for \code{transformCounts}. \code{transformFeatures}  does the transformation 
#' feature-wise, i.e., row-wise. \code{ZTransform} is a shortcut for Z-transformation.
#' \code{relAbundanceCounts} is a shortcut for fetching relative abundance table.
#'
#' @param x A
#'   \code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#'    object.
#'
#' @param abund_values A single character value for selecting the
#'   \code{\link[SummarizedExperiment:SummarizedExperiment-class]{assay}} to be
#'   transformed.
#'
#' @param method A single character value for selecting the transformation
#'   method.
#'
#' @param name A single character value specifying the name of transformed
#'   abundance table.
#'
#' @param pseudocount FALSE or numeric value deciding whether pseudocount is
#'   added. Numerical value specifies the value of pseudocount. (Only used for 
#'   methods \code{method = "clr"}, \code{method = "hellinger"}, or
#'   \code{method = "log10"})
#'
#' @param threshold A numeric value for setting threshold for pa transformation.
#'   By default it is 0. (Only used for \code{method = "pa"})
#'
#' @param ... additional arguments
#'
#' @details
#' \code{transformCounts} or \code{transformSamples} and \code{transformFeatures}
#' applies transformation to abundance table. Provided transformation methods include:
#'
#' \itemize{
#' 
#' \item{'clr'}{ Centered log ratio (clr) transformation can be used for reducing the
#' skewness of data and for centering it. (See e.g. Gloor et al. 2017.)
#'
#' \deqn{clr = log_{10}x_{r} - log_{10}µ_{r}}{%
#' clr = log10 x_r - log10 µ_r}
#' where \eqn{x_{r}} is a single relative value, \eqn{\mu_{r}} is
#' mean relative value".}
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
#' \deqn{log10 = \log_10 x}{%
#' log10 = log10(x)}
#' where \eqn{x} is a single value of data.}
#' 
#' \item{'pa'}{ Transforms table to presence/absence table. All abundances higher
#' than \eqn{\epsilon} are transformed to 1 (present), otherwise 0 (absent). By default, threshold is 0.}
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
#' \deqn{z = \frac{x + \mu}{\sigma}}{%
#' z = (x + µ)/σ}
#' where \eqn{x} is a single value, \eqn{\mu} is the mean of the feature, and
#' \eqn{\sigma} is the standard deviation of the feature.}
#'
#' }
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
#' # By specifying, it is possible to apply different transformations, e.g. clr transformation.
#' # Pseudocount can be added by specifying 'pseudocount'.
#' x <- transformSamples(x, method="clr", pseudocount=1)

#' head(assay(x, "clr"))
#'
#' # Also, the target of transformation
#' # can be specified with "abund_values".
#' x <- transformSamples(x, method="relabundance")
#' x <- transformSamples(x, method="clr", abund_values="relabundance", 
#'                         pseudocount = min(assay(x, "relabundance")[assay(x, "relabundance")>0]))
#' x2 <- transformSamples(x, method="clr", abund_values="counts", pseudocount = 1)
#' head(assay(x, "clr"))
#'
#' # Different pseudocounts used by default for counts and relative abundances
#' x <- transformSamples(x, method="relabundance")
#' mat <- assay(x, "relabundance"); 
#' pseudonumber <- min(mat[mat>0])
#' x <- transformSamples(x, method="clr", abund_values = "relabundance", pseudocount=pseudonumber)
#' x <- transformSamples(x, method="clr", abund_values = "counts", pseudocount=1)
#'
#' # Name of the stored table can be specified.
#' x <- transformSamples(x, method="hellinger", name="test")
#' head(assay(x, "test"))
#'
#' # pa returns presence absence table. With 'threshold', it is possible to set the
#' # threshold to desired level. By default, it is 0.
#' x <- transformSamples(x, method="pa", threshold=35)
#' head(assay(x, "pa"))
#' 
#' # rank returns ranks of taxa. It is calculated column-wise, i.e., per sample
#' # and using the ties.method="first" from the colRanks function
#' x <- transformSamples(x, method="rank")
#' head(assay(x, "rank"))
#' 
#' # transformCounts is an alias for transformSamples
#' x <- transformCounts(x, method="relabundance", name="test2")
#' head(assay(x, "test2"))
#'
#' # In order to use other ranking variants, modify the chosen assay directly:
#' assay(x, "rank_average", withDimnames = FALSE) <- colRanks(assay(x, "counts"), 
#'                                                            ties.method="average", 
#'                                                            preserveShape = TRUE)  
#'                                                            
#' # If you want to do the transformation for features, you can do that by using
#' x <- transformFeatures(x, method="log10", name="log10_features", pseudocount = 1)
#' head(assay(x, "log10_features"))
#'
#' # Z-transform can be done for features by using shortcut function
#' x <- ZTransform(x)
#' head(assay(x, "z"))
#' 
#' # For visualization purposes it is sometimes done CLR for samples, followed by Z transform for taxa
#' x <- ZTransform(transformCounts(x, method="clr", abund_values = "counts", pseudocount = 1))
#'
#' # Relative abundances can be also calculate with the dedicated
#' # relAbundanceCounts function.
#' x <- relAbundanceCounts(x)
#' head(assay(x, "relabundance"))
NULL

#' @rdname transformCounts
#' @aliases transformSamples
#' @export
setGeneric("transformCounts", signature = c("x"),
           function(x,
                    abund_values = "counts",
                    method = c("clr", "hellinger", "log10", "pa", "rank", "relabundance"),
                    name = method,
                    pseudocount = FALSE,
                    threshold = 0)
             standardGeneric("transformCounts"))

#' @rdname transformCounts
#' @aliases transformSamples
#' @export
setMethod("transformCounts", signature = c(x = "SummarizedExperiment"),
          function(x,
                   abund_values = "counts",
                   method = c("clr", "hellinger", "log10", "pa", "rank", "relabundance"),
                   name = method,
                   pseudocount = FALSE,
                   threshold = 0){
        transformSamples(x, 
                         abund_values = abund_values,
                         method = method,
                         name = name,
                         pseudocount = pseudocount,
                         threshold = threshold)
          }
)

#' @rdname transformCounts
#' @aliases transformCounts
#' @export
setGeneric("transformSamples", signature = c("x"),
           function(x,
                    abund_values = "counts",
                    method = c("clr", "hellinger", "log10", "pa", "rank", "relabundance"),
                    name = method,
                    pseudocount = FALSE,
                    threshold = 0)
                    standardGeneric("transformSamples"))

#' @rdname transformCounts
#' @aliases transformCounts
#' @export
setMethod("transformSamples", signature = c(x = "SummarizedExperiment"),
    function(x,
            abund_values = "counts",
            method = c("clr", "hellinger", "log10", "pa", "rank", "relabundance"),
            name = method,
            pseudocount = FALSE,
            threshold = 0){
        # Input check
        # Check abund_values
        .check_assay_present(abund_values, x)
        # Check name
        if(!.is_non_empty_string(name) ||
           name == abund_values){
            stop("'name' must be a non-empty single character value and be ",
                 "different from `abund_values`.",
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
        assay <- assay(x, abund_values)
        # Calls help function that does the transformation
        transformed_table <- .apply_transformation(assay, method, pseudocount, threshold)
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
                    abund_values = "counts",
                    method = c("log10", "pa", "z"),
                    name = method,
                    pseudocount = FALSE,
                    threshold = 0)
               standardGeneric("transformFeatures"))

#' @rdname transformCounts
#' @export
setMethod("transformFeatures", signature = c(x = "SummarizedExperiment"),
    function(x,
             abund_values = "counts",
             method = c("log10", "pa", "z"),
             name = method,
             pseudocount = FALSE,
             threshold = 0){
        # Input check
        # Check abund_values
        .check_assay_present(abund_values, x)
        # Check name
        if(!.is_non_empty_string(name) ||
           name == abund_values){
            stop("'name' must be a non-empty single character value and be ",
                 "different from `abund_values`.",
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
        assay <- t(assay(x, abund_values))
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
        transformCounts(x, method = "relabundance", ...)
    }
)

###########################HELP FUNCTIONS####################################

# Help function for transformSamples and transformFeatures, takes abundance table
# as input and returns transformed table.
.apply_transformation <- function(assay, method, pseudocount, threshold, ...){
    # Input check
    # Check pseudocount
    if(length(pseudocount) != 1L || 
       !(pseudocount == FALSE || is.numeric(pseudocount))){
        stop("'pseudocount' must be FALSE or a single numeric value.",
             call. = FALSE)
    } else if(is.numeric(pseudocount)){
        if (method == "relabundance" && pseudocount > 0){
            warning("Relative abundances vary in [0, 1]; adding a",
                    "pseudocount > 0 on relabundance will cause ",
                    "non-sensicale results. Recommended to cross-check ",
                    "that the pseudocount choice is correct and intended.",
                    call. = FALSE)
        }
        if(pseudocount == 0){
            pseudocount <- FALSE
        }
    }
    # Check threshold
    if(!is.numeric(threshold)){
        stop("'threshold' must be a numeric value, and it can be used ",
             "only with transformation method 'pa'.",
             call. = FALSE)
    }
    # Input check end
    
    # apply pseudocount
    abund <- .apply_pseudocount(assay, pseudocount)
    # Get transformed table
    transformed_table <-
        .get_transformed_table(assay = abund,
                               method = method,
                               threshold = threshold)
    return(transformed_table)
}

# Chooses which transformation function is applied
.get_transformed_table <- function(assay, method, threshold){
    # Function is selected based on the "method" variable
    FUN <- switch(method,
                  relabundance = .calc_rel_abund,
                  log10 = .calc_log10,
                  pa = .calc_pa,
                  hellinger = .calc_hellinger,
                  clr = .calc_clr,
                  rank = .calc_rank,
                  z = .calc_ztransform)

    # Does the function call, arguments are "assay" abundance table and "pseudocount"
    do.call(FUN,
            list(mat = assay,
                 threshold = threshold))
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

#' @importFrom DelayedMatrixStats colMeans2
.calc_clr <- function(mat, ...){
    mat <- .calc_rel_abund(mat)
    # If there is negative values, gives an error.
    if (any(mat <= 0, na.rm = TRUE)) {
        stop("Abundance table contains zero or negative values and ",
             "clr-transformation is being applied without (suitable) ",
             "pseudocount. \n",
             "Try to add pseudocount (default choice pseudocount = 1 for ",
             "count assay; or pseudocount = min(x[x>0]) with relabundance ",
             "assay).",
             call. = FALSE)
    }
    # In every sample, calculates the log of individual entries. After that calculates
    # the sample-specific mean value and subtracts every entries' value with that.
    clog <- log(mat)
    clogm <- colMeans2(clog)
    mat <- t(t(clog) - clogm)
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
    # If "pseudocount" is not FALSE, it is numeric value specified by user. 
    # Then add pseudocount.
    if(!pseudocount == FALSE){
        mat <- mat + pseudocount
    }
    return(mat)
}
