#' Subsample Counts
#' 
#' \code{subsampleCounts} will randomly subsample counts in 
#' \code{SummarizedExperiment} and return the a modified object in which each 
#' sample has same number of total observations/counts/reads. 
#'
#' @details
#' Although the subsampling approach is highly debated in microbiome research, 
#' we include the \code{subsampleCounts} function because there may be some 
#' instances where it can be useful.
#' Note that the output of \code{subsampleCounts} is not the equivalent as the 
#' input and any result have to be verified with the original dataset.
#'
#' Subsampling/Rarefying may undermine downstream analyses and have unintended
#' consequences. Therefore, make sure this normalization is appropriate for
#' your data.
#'
#' To maintain the reproducibility, please define the seed using set.seed() 
#' before implement this function.
#'
#' @param x A \code{SummarizedExperiment} object.
#'
#' @param assay.type A single character value for selecting the
#'   \code{SummarizedExperiment} \code{assay} used for random subsampling. 
#'   Only counts are useful and other transformed data as input will give 
#'   meaningless output.
#'   
#' @param assay_name a single \code{character} value for specifying which
#'   assay to use for calculation.
#'   (Please use \code{assay.type} instead. At some point \code{assay_name}
#'   will be disabled.)
#'   
#' @param min_size A single integer value equal to the number of counts being 
#'   simulated this can equal to lowest number of total counts 
#'   found in a sample or a user specified number. 
#'   
#' @param replace Logical Default is \code{TRUE}. The default is with 
#'   replacement (\code{replace=TRUE}). 
#'   See \code{\link[phyloseq:rarefy_even_depth]{phyloseq::rarefy_even_depth}}
#'   for details on implications of this parameter.   
#' 
#' @param name A single character value specifying the name of transformed
#'   abundance table.
#' 
#' @param verbose Logical Default is \code{TRUE}. When \code{TRUE} an additional 
#'   message about the random number used is printed.
#' 
#' @param ... additional arguments not used
#' 
#' @references
#' McMurdie PJ, Holmes S. Waste not, want not: why rarefying microbiome data 
#' is inadmissible. PLoS computational biology. 2014 Apr 3;10(4):e1003531.
#' 
#' Gloor GB, Macklaim JM, Pawlowsky-Glahn V & Egozcue JJ (2017)
#' Microbiome Datasets Are Compositional: And This Is Not Optional.
#' Frontiers in Microbiology 8: 2224. doi: 10.3389/fmicb.2017.02224
#' 
#' Weiss S, Xu ZZ, Peddada S, Amir A, Bittinger K, Gonzalez A, Lozupone C, 
#' Zaneveld JR, Vázquez-Baeza Y, Birmingham A, Hyde ER. Normalization and 
#' microbial differential abundance strategies depend upon data characteristics. 
#' Microbiome. 2017 Dec;5(1):1-8.
#' 
#' @return \code{subsampleCounts} return \code{x} with subsampled data.
#' 
#' @author Sudarshan A. Shetty and Felix G.M. Ernst
#' 
#' @name subsampleCounts
#'  
#' @examples
#' # When samples in TreeSE are less than specified min_size, they will be
#' # removed. If after subsampling features are not present in any of the
#' # samples, they will be removed.
#' data(GlobalPatterns)
#' tse <- GlobalPatterns
#' set.seed(4759)
#' tse.subsampled <- subsampleCounts(tse, 
#'                                   min_size = 60000, 
#'                                   name = "subsampled")
#' tse.subsampled
#' dim(tse)
#' dim(tse.subsampled)
#' 
NULL

#' @rdname subsampleCounts
#' @aliases rarifyCounts
#' @export
setGeneric("subsampleCounts", signature = c("x"),
    function(x, assay.type = assay_name, assay_name = "counts", 
        min_size = min(colSums2(assay(x, assay.type))),
        replace = TRUE,
        name = "subsampled", verbose = TRUE, ...)
    standardGeneric("subsampleCounts"))

#' @importFrom SummarizedExperiment assay assay<-
#' @importFrom DelayedMatrixStats colSums2 rowSums2
#' @rdname subsampleCounts
#' @aliases rarifyCounts
#' @export
setMethod("subsampleCounts", signature = c(x = "SummarizedExperiment"),
    function(x, assay.type = assay_name, assay_name = "counts", 
            min_size = min(colSums2(assay(x, assay.type))), replace = TRUE, 
            name = "subsampled", verbose = TRUE, ...){
        # Input check
        # CHeck that assay name is correct and that assay is counts table.
        .check_assay_present(assay.type, x)
        if( any(assay(x, assay.type) %% 1 != 0) ){
            warning("assay contains non-integer values. Only counts table ",
                    "is applicable...")
        }
        # Check that verbose and replace are boolean values
        if( !.is_a_bool(verbose) ){
            stop("'verbose' must be TRUE or FALSE.", call. = FALSE)
        }
        if( !.is_a_bool(replace) ){
            stop("`replace` must be TRUE or FALSE.", call. = FALSE)
        } 
        # Check name of new assay
        if( !.is_non_empty_string(name) || name == assay.type ){
            stop("'name' must be a non-empty single character value and be ",
                "different from 'assay.type'.", call. = FALSE)
        }
        # Check min_size. It must be single positive integer value.
        if(!is.numeric(min_size) || length(min_size) != 1 ||
            as.integer(min_size) != min_size && min_size <= 0  ){
            stop("min_size needs to be a positive integer value.")
        }
        # Input check end
        
        # min_size determines the number of reads subsampled from samples.
        # This means that every samples should have at least min_size of reads.
        # If they do not have, drop those samples at this point.
        min_reads <- colSums2(assay(x, assay.type)) < min_size
        if( any(min_reads) ){
            # Get those sample names that we are going to remove due to too
            # small number of reads
            rmsams <- colnames(x)[ min_reads ]
            # Remove sample(s)
            newtse <- x[, !colnames(x) %in% rmsams]
            # Return NULL, if no samples were found after subsampling
            if( ncol(x) == 0 ){
                stop("No samples were found after subsampling. Consider ",
                    "lower 'min_size'.", call. = FALSE)
            }
            # Give message which samples were removed
            if( verbose ){
                message(
                    length(rmsams), " samples removed because they contained ",
                    "fewer reads than `min_size`.")
            }
            
        }
        # Subsample specified assay.
        newassay <- apply(assay(x, assay.type), 2, .subsample_assay,
                        min_size=min_size, replace=replace)
        # Add rownames to new assay. The returned value from .subsample_assay
        # is a vector that do not have feature names.
        rownames(newassay) <- rownames(x)
        # remove features not present in any samples after subsampling
        feat_inc <- rowSums2(newassay) > 0
        newassay <- newassay[feat_inc, ]
        # Give message if some features were dropped
        if( verbose && any(!feat_inc) ){
            message(
                sum(!feat_inc), " features removed because they are not ",
                "present in all samples after subsampling."
                )
        }
        # Subset the TreeSE based on new feature-set
        x <- x[rownames(newassay),]
        # Add new assay to TreeSE
        assay(x, name, withDimnames = FALSE) <- newassay
        # Add info on min_size to metadata
        x <- .add_values_to_metadata(
            x, 
            "subsampleCounts_min_size",
            min_size)
        return(x)
    }
)


# Modified Sub sampling function from phyloseq internals
.subsample_assay <- function(x, min_size, replace){
    # Create replacement species vector
    rarvec <- numeric(length(x))  
    # Perform the sub-sampling. Suppress warnings due to old R compat issue.
    # Also, make sure to avoid errors from x summing to zero, 
    # and there are no observations to sample.
    # The initialization of rarvec above is already sufficient.
    if(sum(x) <= 0){
        # Protect against, and quickly return an empty vector, 
        # if x is already an empty count vector
        return(rarvec)
    }
    if(replace){
        # resample with replacement
        obsvec <- seq_along(x)
        prob <- x
    } else {
        # resample without replacement
        obsvec <- mapply(rep_len, x = seq_along(x), length.out = x)
        obsvec <- unlist(obsvec, use.names = FALSE)
        # use `sample` for subsampling. Hope that obsvec doesn't overflow.
        prob <- NULL
    }
    # Do the sampling of features from the single sample
    suppressWarnings(subsample <- sample(
        obsvec,
        min_size,
        replace = replace,
        prob = prob))
    # Tabulate the results (these are already named by the order in `x`)
    sstab <- table(subsample)
    # Assign the tabulated random subsample values to the species vector
    rarvec[as(names(sstab), "integer")] <- sstab
    # Return abundance vector. Let replacement happen elsewhere.
    return(rarvec)
}
