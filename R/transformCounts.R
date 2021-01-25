#' Transform Counts
#'
#' These functions applies transformation to abundance table. By using
#' \code{transformCounts}, transformed table is in \code{assay}. By using
#' specific \code{ZTransform} function, Z-transformation can be applied for features.
#' \code{relAbundanceCounts} is a shortcut for fetching relative abundance table.
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
#' @param method
#' A single character value for selecting the transformation method.
#'
#' @param name
#' A single character value specifying the name of transformed abundance table.
#'
#' @param pseudocount
#' FALSE or numeric value deciding whether pseudocount is added. Numerical
#' value specifies the value of pseudocount.
#'
#' @details
#' \code{transformCounts} applies transformation to abundance table.
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
#'     In function \code{transformCounts}, it is done per-sample like in all other transformations
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
#' \code{transformCounts} and \code{ZTransform} return \code{x} with additional, transformed
#' abundance table named \code{*name*} in the \code{\link{assay}}.
#' For \code{relAbundanceCounts} a modified \code{x} containing the relative
#' abundances as an assay defined by \code{name}.
#'
#' @seealso
#' \itemize{
#'   \item{\code{\link[mia:relAbundanceCounts]{relAbundanceCounts}}}
#' }
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
#' # assays(x)$clr
#'
#' # Name of the stored table can be specified. Also, the target of transformation
#' # can be specified with "abund_values".
#' x <- transformCounts(x, method="hellinger", name="test", pseudocount=5)
#' # assays(x)$test
#' x <- transformCounts(x, method="Z", abund_values="test")
#' # assays(x)$Z
#' # Z-transform can also be done for features
#' x <- ZTransform(x, pseudocount=1)
#' # assays(x)$ZTransform
#'
#' # Relative abundances can be also fetched with function relAbundanceCounts.
#' x <- relAbundanceCounts(x)
#' #assays(x)$relabundance
#'
NULL

#' @rdname transformCounts
#' @export
setGeneric("transformCounts", signature = c("x"),
           function(x,
                    abund_values = "counts",
                    method = c("relabundance", "log10", "pa", "Z", "hellinger", "clr"),
                    name = method,
                    pseudocount = FALSE)
               standardGeneric("transformCounts"))


#' @rdname transformCounts
#' @export
setMethod("transformCounts", signature = c(x = "SummarizedExperiment"),
          function(x,
                   abund_values = "counts",
                   method = c("relabundance", "log10", "pa", "Z", "hellinger", "clr"),
                   name = method,
                   pseudocount = FALSE){

              # Input check
              # Check abund_values
              .check_abund_values(abund_values, x)

              # Check method
              # If method is not single string, user has not specified transform method,
              # or has given e.g. a vector
              if(!.is_non_empty_string(method)){
                  stop("'method' must be a non-empty single character value.
                       Give one method from the following list:
                       'relabundance', 'log10', 'pa', 'Z', 'hellinger', 'clr'",
                       call. = FALSE)
              }
              method <- match.arg(method)

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

              # If the method is "relabundance", relAbundanceCounts is used
              if(method == "relabundance"){

                  # Get relative abundances stored in assays
                  x <- relAbundanceCounts(x,
                                          abund_values,
                                          name,
                                          pseudocount)
              }
              # If the method is something else, internal function is used
              else{
              # Get transformed table
              transformed_table <- .get_transformed_table(assay = assay(x, abund_values),
                                                          method = method,
                                                          pseudocount = pseudocount)

              # Assign transformed table to assays
              assay(x, name) <- transformed_table

              }


              return(x)

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

              # Get transformed table
              transformed_table <- .get_ztransformed_table(assay = assay(x, abund_values),
                                                          pseudocount = pseudocount)

              # Assign transformed table to assays
              assay(x, name) <- transformed_table

              return(x)
          }
)

###############################relAbundanceCounts###############################

#' @rdname transformCounts
setGeneric("relAbundanceCounts", signature = c("x"),
           function(x,
                    abund_values = "counts",
                    name = "relabundance",
                    pseudocount = FALSE)
               standardGeneric("relAbundanceCounts"))

#' @rdname transformCounts
#' @importFrom SummarizedExperiment assay assay<-
#' @export
setMethod("relAbundanceCounts",signature = c(x = "SummarizedExperiment"),
          function(x,
                   abund_values = "counts",
                   name = "relabundance",
                   pseudocount = FALSE){

              # Input check
              # Check abund_values
              .check_abund_values(abund_values, x)

              # Check name
              if(!.is_non_empty_string(name)){
                  stop("'name' must be a single non-empty character value.",
                       call. = FALSE)
              }

              # Check pseudocount
              if(!(pseudocount==FALSE || is.numeric(pseudocount))){
                  stop("'pseudocount' must be FALSE or numeric value.",
                       call. = FALSE)
              }

              # If "pseudocount" is not FALSE, it is numeric value specified by user. Then add pseudocount.
              if(!pseudocount==FALSE){
                  assay(x, abund_values) <- assay(x, abund_values) + pseudocount
              }


              # Get and store relabundance table
              assay(x, name) <- .calc_rel_abund(assay(x, abund_values))

              return(x)
          }
)

###########################HELP FUNCTIONS####################################


# Chooses which transformation function is applied
.get_transformed_table <- function(assay, method, pseudocount){

    # If "pseudocount" is not FALSE, it is numeric value specified by user. Then add pseudocount.
    if(!pseudocount==FALSE){
        assay <- assay + pseudocount
    }

    # Function is selected based on the "method" variable
    FUN <- switch(method,
                  Z = .get_z_table,
                  log10 = .get_log10_table,
                  pa = .get_pa_table,
                  hellinger = .get_hellinger_table,
                  clr = .get_clr_table)

    # Does the function call, arguments are "assay" abundance table and "pseudocount"
    do.call(FUN,
            list(assay = assay))
}

.calc_rel_abund <- function(assay){

    # If assay contains missing values, replace their values so it is possible to calculate
    # relative abundances. Otherwise, whole column would get NAs.
    # Stores the indices of missing values
    missingValues <- which(is.na(assay))
    # Missing values are replaced with value 0
    assay[missingValues] <- 0

    # Calculates the relative abundances. Uses internal function from relabundance.R.
    mat <- sweep(assay, 2, colSums(assay), "/")

    # If there were missing values, add them back to data
    mat[missingValues] <- NA

    return(mat)
}

.get_z_table <- function(assay){

    # Log10 can not be calculated if there is zero
    if (any(assay == 0, na.rm = TRUE)) {
        stop("Abundance table contains zero and Z transformation
                is being applied without pseudocount. Try Z with
                pseudocount = 1 or other numeric value.")
    }
    # If there is negative values, gives an error.
    if (any(assay < 0, na.rm = TRUE)) {
        stop("Abundance table contains negative values and Z-transformation
                is being applied without pseudocount. Try z with
                pseudocount = 1 or other numeric value.")
    }

    # Gets the log transformed table
    mat <- .get_log10_table(assay)

    # Performs z transformation
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
        }
    if (any(assay < 0, na.rm = TRUE)) {
            stop("Abundance table contains negative values and log10 transformation
                is being applied without pseudocount. Try log10 with
                pseudocount = 1 or other numeric value.")
    }

    mat <- log10(assay)

    return(mat)
}

.get_pa_table <- function(assay){

    # If value is over zero, gets value 1. If value is zero, gets value 0.
    mat <- ifelse(assay > 0, 1, 0)

    return(mat)
}

.get_hellinger_table <- function(assay){

    # If there is negative values, gives an error.
    if (any(assay < 0, na.rm = TRUE)) {
        stop("Abundance table contains negative values and hellinger transformation
                is being applied without pseudocount. Try hellinger with
                pseudocount = 1 or other numeric value.")
    }

    # Gets the relative abundance
    mat <- .calc_rel_abund(assay)

    # Takes square root
    mat <- sqrt(mat)

    return(mat)
}

.get_clr_table <- function(assay){

    # If there is negative values, gives an error.
    if (any(assay < 0, na.rm = TRUE)) {
        stop("Abundance table contains negative values and clr transformation
                is being applied without pseudocount. Try clr with
                pseudocount = 1 or other numeric value.")
    }
    # If abundance table contains zeros, gives an error
    if (any(assay == 0, na.rm = TRUE)) {
        stop("Abundance table contains zero and clr transformation
                is being applied without pseudocount. Try clr with
                pseudocount = 1 or other numeric value.")
    }

    mat <- .calc_rel_abund(assay)

    # In every sample, calculates the log of individual entries. After that calculates
    # the sample-specific mean value and subtracts every entries' value with that.
    mat <- apply(mat, 2, function(x) {
        log(x) - mean(log(x), na.rm=TRUE)
    })

    return(mat)
}

.get_ztransformed_table <- function(assay, pseudocount){

    # If "pseudocount" is not FALSE, it is numeric value specified by user. Then add pseudocount.
    if(!pseudocount==FALSE){
        assay <- assay + pseudocount
    }

    # Log10 can not be calculated if there is zero
    if (any(assay == 0, na.rm = TRUE)) {
        stop("Abundance table contains zero and Z-transformation
            is being applied without pseudocount. Try ZTransform with
            pseudocount = 1 or other numeric value.")
    }
    # If there is negative values, gives an error.
    if (any(assay < 0, na.rm = TRUE)) {
        stop("Abundance table contains negative values and Z-transformation
                is being applied without pseudocount. Try ZTransform with
                pseudocount = 1 or other numeric value.")
    }

    # Gets the log transformed table
    mat <- .get_log10_table(assay)

    # Transpoese the table
    trans <- t(mat)

    # Z transform for features. Centers the feature data. After that, divides with
    # the standard deviation of feature.
    # Performs z transformation
    trans <- apply(trans, 2, function(x) {
        (x - mean(x, na.rm=TRUE))/sd(x, na.rm=TRUE)
    })

    # Transposes the table back to its original form
    trans <- t(trans)

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

    }

    return(mat)

}

