#' Transform Counts
#'
#' These functions provide a variety of options for transforming abundance data.
#' By using \code{transformCounts}, transformed table is in \code{assay}. By using
#' specific \code{ZTransform} function, Z-transformation can be applied for features.
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
#'   methods \code{method = "log10"}, \code{method = "hellinger"} or 
#'   \code{method = "clr"})
#'
#' @param threshold A numeric value for setting threshold for pa transformation.
#'   By default it is 0. (Only used for \code{method = "pa"})
#'
#' @param ... additional arguments
#'
#' @details
#' \code{transformCounts} applies transformation to abundance table.
#' Provided transformation methods include:
#'
#' \itemize{
#' \item {'relabundance'}{ Transforms abundances to relative. Generally, all microbiome
#' data are compositional. That is, e.g., because all measuring instruments have their capacity limits.
#' To make results comparable with other results, values must be relative. (See e.g. Gloor et al. 2017.)
#'
#' \deqn{relabundance = \frac{x}{x_{tot}}}{%
#' relabundance = x/x_tot}
#' where \eqn{x} is a single value and \eqn{x_{tot}} is the sum of
#' all values.}
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
#' \item{'Z'}{ Z-transformation, Z score transformation, or Z-standardization normalizes
#' the data by shifting (to mean \eqn{\mu}) and scaling (to standard deviation \eqn{\sigma}).
#' Z-transformation can be done with function \code{ZTransform}. It is done per rows (features / taxa),
#' unlike most other transformations. This is often preceded by log10p or clr transformation.
#' In other words, single value is standardized with respect of feature's values.
#'
#' \deqn{Z = \frac{x + \mu}{\sigma}}{%
#' Z = (x + µ)/σ}
#' where \eqn{x} is a single value, \eqn{\mu} is the mean of the feature, and
#' \eqn{\sigma} is the standard deviation of the feature.}
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
#' \item{'clr'}{ Centered log ratio (clr) transformation can be used for reducing the
#' skewness of data and for centering it. (See e.g. Gloor et al. 2017.)
#'
#' \deqn{clr = log_{10}x_{r} - log_{10}µ_{r}}{%
#' clr = log10 x_r - log10 µ_r}
#' where \eqn{x_{r}} is a single relative value, \eqn{\mu_{r}} is
#' mean relative value".}
#' 
#' \item{'rank'}{ Rank returns ranks of taxa. For each sample, the least abundant 
#' taxa get lower value and more abundant taxa bigger value. The implementation is 
#' based on the colRanks function with ties.method="first".}
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
#' \code{transformCounts}, \code{relAbundanceCounts}, and \code{ZTransform} return
#' \code{x} with additional, transformed abundance table named \code{*name*} in
#' the \code{\link{assay}}.
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
#' x <- transformCounts(x, method="clr", pseudocount=1)
#' head(assay(x, "clr"))
#'
#' # Also, the target of transformation
#' # can be specified with "abund_values".
#' x <- transformCounts(x, method="relabundance")
#' x <- transformCounts(x, method="clr", abund_values="relabundance", 
#'                         pseudocount = min(assay(x, "relabundance")[assay(x, "relabundance")>0]))
#' x2 <- transformCounts(x, method="clr", abund_values="counts", pseudocount = 1)
#' head(assay(x, "clr"))
#'
#' # Different pseudocounts used by default for counts and relative abundances
#' x <- transformCounts(x, method="relabundance")
#' mat <- assay(x, "relabundance"); 
#' pseudonumber <- min(mat[mat>0])
#' x <- transformCounts(x, method="clr", abund_values = "relabundance", pseudocount=pseudonumber)
#' x <- transformCounts(x, method="clr", abund_values = "counts", pseudocount=1)
#'
#' # Name of the stored table can be specified. 
#' x <- transformCounts(x, method="hellinger", name="test")
#' head(assay(x, "test"))
#'
#' # pa returns presence absence table. With 'threshold', it is possible to set the
#' # threshold to desired level. By default, it is 0.
#' x <- transformCounts(x, method="pa", threshold=35)
#' head(assay(x, "pa"))
#' 
#' # rank returns ranks of taxa. It is calculated column-wise, i.e., per sample
#' # and using the ties.method="first" from the colRanks function
#' x <- transformCounts(x, method="rank")
#' head(assay(x, "rank"))
#'
#' # In order to use other ranking variants, modify the chosen assay directly:
#' assay(x, "rank_average", withDimnames = FALSE) <- t(colRanks(assay(x, "counts"), ties.method="average"))
#'
#' # Z-transform can be done for features, not for samples as in the other transformations
#' x <- ZTransform(x)
#' head(assay(x, "ZTransform"))
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
#' @export
setGeneric("transformCounts", signature = c("x"),
           function(x,
                    abund_values = "counts",
                    method = c("relabundance", "log10", "pa", "hellinger", "clr", "rank"),
                    name = method,
                    pseudocount = FALSE,
                    threshold = 0)
               standardGeneric("transformCounts"))

#' @rdname transformCounts
#' @export
setMethod("transformCounts", signature = c(x = "SummarizedExperiment"),
    function(x,
             abund_values = "counts",
             method = c("relabundance", "log10", "pa", "hellinger", "clr", "rank"),
             name = method,
             pseudocount = FALSE,
             threshold = 0){
        # Input check
        # Check abund_values
        .check_assay_present(abund_values, x)

        # Check method
        # If method is not single string, user has not specified transform method,
        # or has given e.g. a vector
        if(!.is_non_empty_string(method)){
            stop("'method' must be a non-empty single character value. \n",
                 "Give one method from the following list: \n",
                 "'relabundance', 'log10', 'pa', 'hellinger', 'clr', 'rank'",
                 call. = FALSE)
        }
        method <- match.arg(method)

        # Check name
        if(!.is_non_empty_string(name) ||
           name == abund_values){
            stop("'name' must be a non-empty single character value and be ",
                 "different from `abund_values`.",
                 call. = FALSE)
        }

        # Check pseudocount
        if(!(pseudocount==FALSE || is.numeric(pseudocount))){
            stop("'pseudocount' must be FALSE or numeric value.",
                 call. = FALSE)
        }

        # Check threshold
        if(!is.numeric(threshold)){
            stop("'threshold' must be a numeric value, and it can be used ",
                 "only with transformation method 'pa'.",
                 call. = FALSE)
        }

        # apply pseudocount
        abund <- .apply_pseudocount(assay(x, abund_values), pseudocount)

        # Get transformed table
        transformed_table <-
            .get_transformed_table(assay = abund,
                                   method = method,
                                   threshold = threshold)

        # Assign transformed table to assays
        assay(x, name, withDimnames=FALSE) <- transformed_table
        x
    }
)

##################################Z-TRANSFORM###################################

#' @rdname transformCounts
#' @export
setGeneric("ZTransform", signature = c("x"),
           function(x,
                    abund_values = "counts",
                    name = "ZTransform",
                    pseudocount = FALSE)
               standardGeneric("ZTransform"))


#' @rdname transformCounts
#' @export
setMethod("ZTransform", signature = c(x = "SummarizedExperiment"),
    function(x,
             abund_values = "counts",
             name = "ZTransform",
             pseudocount = FALSE){

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

        # Check pseudocount
        if(!(pseudocount==FALSE || is.numeric(pseudocount))){
          stop("'pseudocount' must be FALSE or numeric value.",
               call. = FALSE)
        }

        # apply pseudocount
        mat <- .apply_pseudocount(assay(x, abund_values), pseudocount)

        # Get transformed table
        transformed_table <- .calc_ztransform(mat = mat)

        # Assign transformed table to assays
        assay(x, name) <- transformed_table
        x
    }
)

###############################relAbundanceCounts###############################

#' @rdname transformCounts
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


# Chooses which transformation function is applied
.get_transformed_table <- function(assay, method, threshold){
    # Function is selected based on the "method" variable
    FUN <- switch(method,
                  relabundance = .calc_rel_abund,
                  log10 = .calc_log10,
                  pa = .calc_pa,
                  hellinger = .calc_hellinger,
                  clr = .calc_clr,
                  rank = .calc_rank)

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
             "log10 transformation is being applied without pseudocount. ",
             "Try to add pseudocount (default choice pseudocount = 1 for count assay;
	      or pseudocount = min(x[x>0]) for relabundance assay).",
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
             "transformation is being applied without (suitable) pseudocount. ",
             "Try to add pseudocount (default choice pseudocount = 1 for count assay;
	      or pseudocount = min(x[x>0]) with relabundance assay).",
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
             "pseudocount.",
             "Try to add pseudocount (default choice pseudocount = 1 for count assay;
	      or pseudocount = min(x[x>0]) with relabundance assay).",
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
    # For every sample, finds ranks of taxa.
    # Column-wise, NAs are kept as NAs, and ties get the minimum rank value.
    # Transpose ensures that dimensions of matrix are right.
    mat <- t(colRanks(mat, ties.method="first"))
    return(mat)
}

#' @importFrom DelayedMatrixStats rowMeans2 rowSds
.calc_ztransform <- function(mat){

    # Z transform for features
    # Centers the feature data. After that, divides with
    # the standard deviation of feature.
    rm <- rowMeans2(mat, na.rm = TRUE)
    rsd <- rowSds(mat, na.rm = TRUE)
    mat <- (mat - rm)/rsd
    return(mat)
}

.apply_pseudocount <- function(mat, pseudocount){
    # If "pseudocount" is not FALSE, it is numeric value specified by user. 
    # Then add pseudocount.
    if(!pseudocount==FALSE){
        mat <- mat + pseudocount
    }
    mat
}
