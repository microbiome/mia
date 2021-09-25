#' Subsample a \code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#' object
#' 
#' \code{getSubsample} will randomly subsample counts in SE/TSE to return a SE/TSE in
#' which each sample has same number of total observations/counts/reads. 
#'
#' @details
#' Although the subsampling approach is highly debated in microbiome research, 
#' we include the \code{getSubsample} function because in some instances these 
#' can be useful. Note that the output of \code{getSubsample} is a modified SE/TSE.
#'
#' @param x A
#'   \code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#'    object.
#'
#' @param abund_values A single character value for selecting the
#'   \code{\link[SummarizedExperiment:SummarizedExperiment-class]{assay}} used for
#'   random subsampling. Only counts are useful and other transformed data as input
#'   will give meaningless output.
#'   
#' @param min_size A single integer value equal to the number of counts being 
#'   simulated this can equal to lowest number of total counts 
#'   found in a sample or a user specified number. 
#'   
#' @param seed A random number seed for reproducibility of sampling. 
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
#' @return a TSE with subsampled data that is different from input TSE.
#' 
#' @author Sudarshan A. Shetty 
#' 
#' @name getSubsample
#'  
#' @examples
#' # When samples in TSE are less than specified min_size, they will be removed.
#' # If after subsampling features are not present in any of the samples, 
#' # they will be removed.
#' data("GlobalPatterns")
#' tse <- GlobalPatterns
#' tse.subsampled <- getSubsample(tse, min_size = 60000, name = "subsampled")
#' tse.subsampled
#' 
#' dim(tse)
#' 
#' dim(tse.subsampled)
#' 
NULL

#' @rdname getSubsample
#' @export
setGeneric("getSubsample", signature = c("x"),
           function(x, 
                    abund_values = "counts",
                    min_size = min(colSums2(assay(x))),
                    seed = runif(1, 0, .Machine$integer.max),
                    replace = TRUE,
                    name = "rarefied",
                    verbose = TRUE, ...)
             standardGeneric("getSubsample"))


#' @importFrom SummarizedExperiment assay assay<-
#' @importFrom DelayedMatrixStats colSums2 rowSums2
#' @rdname getSubsample
#' @export
setMethod("getSubsample", signature = c(x = "SummarizedExperiment"),
          function(x, 
                   abund_values = "counts",
                   min_size = min(colSums2(assay(x))),
                   seed = runif(1, 0, .Machine$integer.max), 
                   replace = TRUE, 
                   name = "subsampled", 
                   verbose = TRUE, ...){
            # Input check
            .check_assay_present(abund_values, x)
            if(verbose){
              # Print to screen this value
              message("`set.seed(", seed, ")` was used to initialize repeatable random subsampling.")
              message("Please record this for your records so others can reproduce. \n ... \n")
            }
            # Make sure min_size is of length 1.
            if(length(min_size) > 1){
              warning("`min_size` had more than one value. ", 
                      "Using only the first. \n ... \n")
              min_size <- min_size[1]	
            }
            if(min_size <= 0){
              stop("min_size less than or equal to zero. ", 
                   "Need positive sample size to work.")
            }
            # get samples with less than min number of reads
            if(min(colSums2(assay(x, abund_values))) < min_size){
              rmsams <- colnames(x)[colSums2(assay(x, abund_values)) < min_size]
              if(verbose){
                message(length(rmsams), " samples removed ",
                        "because they contained fewer reads than `min_size`.")
              }
              # remove sample(s)
              newtse <- x[, !colnames(x) %in% rmsams]
            }
            newassay <- apply(assay(newtse, abund_values), 2, 
                              .subsample_assay,
                              min_size=min_size, replace=replace)
            rownames(newassay) <- rownames(newtse)
            # remove features not present in any samples after subsampling
            message(paste(length(which(rowSums2(newassay) == 0)), "features", 
                          "removed becasue they are not present in all samples", 
                          "after subsampling.\n"))
            # get features features with non-zero sum across samples.
            keepfeatures <- rownames(newassay[which(rowSums2(newassay) != 0),])
            # add the subsampled assay
            assay(newtse, name, withDimnames=FALSE) <- newassay
            # filter tse to keep only features with non-zero sum across samples
            newtse <- newtse[keepfeatures,]
            newtse
          }
)


# Sub sampling function from phyloseq internals
.subsample_assay <- function(x, min_size, replace=FALSE){
  
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
    suppressWarnings(subsample <- sample(1:length(x), min_size, replace=TRUE, prob=x))
  } else {
    # resample without replacement
    obsvec <- apply(data.frame(featuresi=1:length(x), times=x), 1, function(x){
      rep_len(x["featuresi"], x["times"])
    })
    obsvec <- unlist(obsvec, use.names=FALSE)
    # use `sample` for subsampling. Hope that obsvec doesn't overflow.
    suppressWarnings(subsample <- sample(obsvec, min_size, replace=FALSE))
  }
  # Tabulate the results (these are already named by the order in `x`)
  sstab <- table(subsample)
  # Assign the tabulated random subsample values to the species vector
  rarvec[as(names(sstab), "integer")] <- sstab
  # Return abundance vector. Let replacement happen elsewhere.
  return(rarvec)
}
