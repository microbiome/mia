#' Transform Abundance
#'
#' These functions applies transformation to abundance table. By using
#' \code{transformAbundance}, transformed table is in \code{assay}. By using
#' specific \code{ZTransform} function, Z-transformation can be applied for features.
#'
#' @param x
#' A \code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#' object.
#'
#' @param abund_values
#' A single character value for selecting the
#'   \code{\link[SummarizedExperiment:SummarizedExperiment-class]{assay}}
#'   to be transformed.
#'
#' @param transform
#' A single character value for selecting the transformation method.
#'
#' @param name
#' A single character value specifying the name of transformed abundance table.
#'
#' @param pseudocount
#' FALSE or numeric value deciding whether pseudocount is added. Numerical
#' value specifies the value of pseudocount.
#'
#' @param scalingfactor
#' A numeric value for scaling. Values are multiplied with specified value.
#'
#' @details
#' \code{transformAbundance} applies transformation to abundance table.
#' Transformation methods include
#'
#' \itemize{
#'
#'     \item{"relabundance" }{Transforms abundances to relative. Generally, all microbiome
#'     data are compositional. That is, e.g., because all measuring instruments have their capacity limits.
#'     To make results comparable with other results, values must be relative.
#'     ($\frac{x}{x_{tot}}$, where $x$ is a single value and $x_{tot}$ is the sum of
#'     all values.) (Gloor et al. 2017.)}
#'
#'     \item{"log10" }{log10 transformation can be used for reducing the skewness of the data.
#'     ($log_{10}x$, where $x$ is a single value of data.)}
#'
#'     \item{"pa" }{Transforms table to presence/absence table. If value is over 0,
#'     then value is 1. If value is 0, then value is 0.)}
#'
#'     \item{"Z" }{Z-transformation or Z-standardization can be used for normalizing the data.
#'     In function \code{transformAbundance}, it is done per-sample like in all other transformations
#'     that are provided by this function. As a special case, Z-transformation can also be done
#'     per-features by using function \code{ZTransform}. However, Z-transformation
#'     done for features can give misleading results.
#'     ($\frac{x + µ}{σ}$, where $x$ is a single value, $µ$ is the mean of the sample, and
#'     $σ$ is the standard deviation of the sample.)}
#'
#'     \item{"hellinger" }{Hellinger transformation can be used for reducing the impact of
#'     extreme data points. It can be utilize for clustering or ordination analysis.
#'     ($\sqrt{\frac{x}{x_{tot}}}$, where $x$ is a single value and $x_{tot}$ is the sum of
#'     all values.) (Legendre & Gallagher 2001.)}
#'
#'     \item{"clr" }{Centered log ratio (clr) transformation can be used for reducing the
#'     skewness of data and for centering it.
#'     ($log_{10}x_{r}) - log_{10}µ_{r}$, where $x_{r}$ is a single relative value, $µ_{r}$ is
#'     the mean of relative values of whole sample.) (Gloor et al. 2017.)}
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
#' \code{transformAbundance} and \code{ZTransform} return \code{x} with additional, transformed
#' abundance table named \code{*name*} in the \code{\link{assay}}.
#'
#' @seealso
#' \itemize{
#'   \item{\code{\link[mia:relAbundanceCounts]{relAbundanceCounts}}}
#' }
#'
#' @name transformAbundance
#' @export
#'
#' @author Leo Lahti and Tuomas Borman. Contact: \url{microbiome.github.io}
#'
#' @examples
#' data(esophagus)
#' x <- esophagus
#'
#' # By specifying, it is possible to apply different transformations, e.g. clr transformation.
#' # Pseudocount can be added by giving the value TRUE. Then pseudocount is 1.
#' x <- transformAbundance(x, transform="clr", pseudocount=TRUE)
#' assays(x)$clr
#'
#' # Name of the stored table can be specified. Also, the target of transform
#' # can be specified with "abund_values". In addition to pseudocount=TRUE, it is
#' # also possible to specify it by giving numeric value.
#' x <- transformAbundance(x, transform="hellinger", name="test", pseudocount=5)
#' assays(x)$test
#' x <- transformAbundance(x, transform="Z", abund_values="test")
#' assays(x)$Z
#' # Z-transform can also be done for features
#' x <- ZTransform(x, pseudocount=TRUE)
#' assays(x)$ZTransform
#'
NULL

#' @rdname transformAbundance
#' @export
setGeneric("transformAbundance", signature = c("x"),
           function(x,
                    abund_values = "counts",
                    transform = c("relabundance", "log10", "pa", "Z", "hellinger", "clr"),
                    name = transform,
                    pseudocount = FALSE,
                    scalingFactor = 1)
               standardGeneric("transformAbundance"))


#' @rdname transformAbundance
#' @export
setMethod("transformAbundance", signature = c(x = "SummarizedExperiment"),
          function(x,
                   abund_values = "counts",
                   transform = c("relabundance", "log10", "pa", "Z", "hellinger", "clr"),
                   name = transform,
                   pseudocount = FALSE,
                   scalingFactor = 1){

              # Input check
              # Check abund_values
              .check_abund_values(abund_values, x)

              # Check transform
              # If transform is not single string, user has not probably specified transform method,
              # or has given a vector
              if(!.is_non_empty_string(name)){
                  stop("'transform' must be a non-empty single character value.
                       Give one method from the following list:
                       'relabundance', 'log10', 'pa', 'Z', 'hellinger', 'clr'",
                       call. = FALSE)
              }
              transform <- match.arg(transform)

              # Check name
              if(!.is_non_empty_string(name)){
                  stop("'name' must be a non-empty single character value.",
                       call. = FALSE)
              }

              # Check pseudocount
              if(!(pseudocount==FALSE || is.numeric(pseudocount))){
                  stop("'pseudocount' must be FALSE or numeric value.",
                       call. = FALSE)
              }

              # Check scalingFactor
              if(!is.numeric(scalingFactor)){
                  stop("'scalingFactor' must be a numeric value.",
                       call. = FALSE)
              }

              # Get transformed table
              transformed_table <- .get_transformed_table(assay = assay(x, abund_values),
                                                          transform = transform,
                                                          pseudocount = pseudocount,
                                                          scalingFactor = scalingFactor)

              # Assign transformed table to assays
              assay(x, name) <- transformed_table

              return(x)

          }
)

##################################Z-TRANSFORM###################################

#' @rdname transformAbundance
#' @export
setGeneric("ZTransform", signature = c("x"),
           function(x,
                    abund_values = "counts",
                    name = "ZTransform",
                    pseudocount = FALSE,
                    scalingFactor = 1)
               standardGeneric("ZTransform"))


#' @rdname transformAbundance
#' @export
setMethod("ZTransform", signature = c(x = "SummarizedExperiment"),
          function(x,
                   abund_values = "counts",
                   name = "ZTransform",
                   pseudocount = FALSE,
                   scalingFactor = 1){

              # Input check
              # Check abund_values
              .check_abund_values(abund_values, x)

              # Check name
              if(!.is_non_empty_string(name)){
                  stop("'name' must be a non-empty single character value.",
                       call. = FALSE)
              }

              # Check pseudocount
              if(!(pseudocount==FALSE || is.numeric(pseudocount))){
                  stop("'pseudocount' must be FALSE or numeric value.",
                       call. = FALSE)
              }

              # Check scalingFactor
              if(!is.numeric(scalingFactor)){
                  stop("'scalingFactor' must be a numeric value.",
                       call. = FALSE)
              }

              # Get transformed table
              transformed_table <- .get_ztransformed_table(assay = assay(x, abund_values),
                                                          pseudocount = pseudocount,
                                                          scalingFactor = scalingFactor)

              # Assign transformed table to assays
              assay(x, name) <- transformed_table

              return(x)
          }
)

###########################HELP FUNCTIONS####################################


# Chooses which transformation function is applied
.get_transformed_table <- function(assay, transform, pseudocount, scalingFactor){

    # If "pseudocount" is not FALSE, it is numeric value specified by user. Then add pseudocount.
    if(!pseudocount==FALSE){
        assay <- assay + pseudocount
    }

    # Multiply values with scalingFactor. By default, scalingFactor is 1, so no changes are made
    assay <- assay * scalingFactor

    # Function is selected based on the "transform" variable
    FUN <- switch(transform,
                  relabundance = .get_relabundance_table,
                  Z = .get_z_table,
                  log10 = .get_log10_table,
                  pa = .get_pa_table,
                  hellinger = .get_hellinger_table,
                  clr = .get_clr_table)

    # Does the function call, arguments are "assay" abundance table and "pseudocount"
    do.call(FUN,
            list(assay = assay))
}

.get_relabundance_table <- function(assay){

    # If assay contains missing values, replace their values so it is possible to calculate
    # Stores the indices of missing values
    missingValues <- which(is.na(assay))
    # Missing values are replaced with value 0
    assay[missingValues] <- 0

    # Calculates the relative abundances. Uses internal function from relabundance.R.
    mat <- .calc_rel_abund(assay)

    # If there were missing values, add them back to data
    mat[missingValues] <- NA

    return(mat)
}

.get_z_table <- function(assay){

    # Log10 can not be calculated if there is zero
    if (any(assay == 0, na.rm = TRUE)) {
        stop("Abundance table contains zero and Z transformation
                is being applied without pseudocount. Try Z with
                pseudocount=TRUE.")
    }

    # Gets the log transformed table
    mat <- .get_log10_table(assay)

    # Performs z transformation for samples
    mat <- apply(mat, 2, function(x) {
        (x - mean(x, na.rm=TRUE))/sd(x, na.rm=TRUE)
    })

    return(mat)
}

.get_log10_table <- function(assay){

    # If abundance table contains zeros, gives an error, because it is not possible
    # to calculate log from zeros. If there is no zeros, calculates log.
    if (any(assay == 0, na.rm = TRUE)) {
        stop("Abundance table contains zero and log10 transformation
                is being applied without pseudocount. Try log10 with
                pseudocount = 1 or other numeric value.")
    } else {
        mat <- log10(assay)
    }
    return(mat)
}

.get_pa_table <- function(assay){

    # If value is over zero, gets value 1. If value is zero, gets value 0.
    mat <- ifelse(assay > 0, 1, 0)

    return(mat)
}

.get_hellinger_table <- function(assay){

    # Gets the relative abundance
    mat <- .get_relabundance_table(assay)

    # Takes square root
    mat <- sqrt(mat)

    return(mat)
}

.get_clr_table <- function(assay){

    # If there is negative values, gives an error.
    if (any(assay < 0, na.rm = TRUE)) {
        stop("Abundance table contains negative values and clr transformation
                is being applied without pseudocount. Try clr with
                pseudocount=TRUE.")
    }
    # If abundance table contains zeros, gives an error
    if (any(assay == 0, na.rm = TRUE)) {
        stop("Abundance table contains zero and clr transformation
                is being applied without pseudocount. Try clr with
                pseudocount = 1 or other numeric value.")
    }

    mat <- .get_relabundance_table(assay)

    # Stores row and col names
    row_names <- rownames(mat)
    col_names <- colnames(mat)

    # In every sample, calculates the log of individual entries. After that calculates
    # the sample-specific mean value and subtracts every entries' value with that.
    mat <- apply(mat, 2, function(x) {
        log(x) - mean(log(x), na.rm=TRUE)
    })

    # Adds sample names to calculated table
    rownames(mat) <- row_names
    # Adds entry names to calculated table
    colnames(mat) <- col_names

    return(mat)
}

.get_ztransformed_table <- function(assay, pseudocount, scalingFactor){

    # If "pseudocount" is not FALSE, it is numeric value specified by user. Then add pseudocount.
    if(!pseudocount==FALSE){
        assay <- assay + pseudocount
    }

    # Multiply values with scalingFactor. By default, scalingFactor is 1, so no changes are made
    assay <- assay * scalingFactor


    # Log10 can not be calculated if there is zero
    if (any(assay == 0, na.rm = TRUE)) {
        stop("Abundance table contains zero and Z transformation
            is being applied without pseudocount. Try ZTransform with
            pseudocount = 1 or other numeric value.")
    }

    # Gets the log transformed table
    mat <- .get_log10_table(assay)

    # Z transform for features. Centers the feature data. After that, divides with
    # the standard deviation of feature.
    trans <- t(scale(t(mat)))
    # Saves all features that are undetectable in every sample i.e. are NA
    undetectables <- which(rowMeans(is.na(trans)) == 1)

    # If there are some features that are undetectable
    # i.e. "undetectable" has length over 0, and table have zeros
    if (length(undetectables) > 0 & min(mat, na.rm=TRUE) == 0) {

        warning("Some features were not detectable. In all samples, signal was
            under detectable limit.")

        # Some features are undetectable in all samples, they are NA when scaled.
        # 0 is assigned to those features.

        trans[names(undetectables), ] <- 0

        mat <- trans

        # Deletes extra information
        attr(mat,"scaled:scale") <- NULL
        attr(mat,"scaled:center") <- NULL
    }

    return(mat)

}
