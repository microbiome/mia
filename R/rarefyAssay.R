#' Subsample Counts
#' 
#' \code{rarefyAssay} randomly subsamples counts within a 
#' \code{SummarizedExperiment} object and returns a new 
#' \code{SummarizedExperiment} containing the original assay and the new 
#' subsampled assay.
#'
#' @details
#' Although the subsampling approach is highly debated in microbiome research, 
#' we include the \code{rarefyAssay} function because there may be some 
#' instances where it can be useful.
#' Note that the output of \code{rarefyAssay} is not the equivalent as the 
#' input and any result have to be verified with the original dataset.
#' To maintain the reproducibility, please define the seed using set.seed() 
#' before implement this function.
#'
#' @param x A
#'   \code{SummarizedExperiment} object.
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
#'   abundance table that will be added to the new \code{SummarizedExperiment}.
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
#' @return \code{rarefyAssay} return \code{x} with subsampled data.
#' 
#' @author Sudarshan A. Shetty and Felix G.M. Ernst
#' 
#' @name rarefyAssay
#'  
#' @examples
#' # When samples in TreeSE are less than specified min_size, they will be removed.
#' # If after subsampling features are not present in any of the samples, 
#' # they will be removed.
#' data(GlobalPatterns)
#' tse <- GlobalPatterns
#' set.seed(123)
#' tse.subsampled <- rarefyAssay(tse, 
#'                                   min_size = 60000, 
#'                                   name = "subsampled" 
#'                                   )
#' tse.subsampled
#' dim(tse)
#' dim(assay(tse.subsampled, "subsampled"))
#' 
NULL

#' @rdname rarefyAssay
#' @export
setGeneric("rarefyAssay", signature = c("x"),
           function(x, assay.type = assay_name, assay_name = "counts", 
                    min_size = min(colSums2(assay(x))),
                    replace = TRUE,
                    name = "subsampled", verbose = TRUE, ...)
               standardGeneric("rarefyAssay"))

#' @importFrom SummarizedExperiment assay assay<-
#' @importFrom DelayedMatrixStats colSums2 rowSums2
#' @rdname rarefyAssay
#' @export
setMethod("rarefyAssay", signature = c(x = "SummarizedExperiment"),
          function(x, assay.type = assay_name, assay_name = "counts", 
                   min_size = min(colSums2(assay(x))),
                    replace = TRUE, 
                   name = "subsampled", verbose = TRUE, ...){
              
              warning("Subsampling/Rarefying may undermine downstream analyses ",
                      "and have unintended consequences. Therefore, make sure ",
                      "this normalization is appropriate for your data.",
                      call. = FALSE)
              .check_assay_present(assay.type, x)
              if(any(assay(x, assay.type) %% 1 != 0)){
                  warning("assay contains non-integer values. Only counts table ",
                          "is applicable...", call. = FALSE)
              }
              if(any(assay(x, assay.type) < 0)){
                  stop("assay contains strictly-negative values. Only counts ",
                      "table is applicable...", call. = FALSE)
              }
              if(!is.logical(verbose)){
                  stop("`verbose` has to be logical i.e. TRUE or FALSE", 
                      call. = FALSE)
              }
            
              if(!is.logical(replace)){
                  stop("`replace` has to be logical i.e. TRUE or FALSE",
                      call. = FALSE)
              } 
              # Check name
              if(!.is_non_empty_string(name) ||
                 name == assay.type){
                  stop("'name' must be a non-empty single character value and be ",
                      "different from `assay.type`.",
                      call. = FALSE)
              }
              #set.seed(seed)
              # Make sure min_size is of length 1.
              if(length(min_size) > 1){
                  stop("`min_size` had more than one value. ", 
                      "Specifiy a single integer value.", call. = FALSE)
                  min_size <- min_size[1]    
              }
              if(!is.numeric(min_size) || 
                 as.integer(min_size) != min_size && min_size <= 0){
                  stop("min_size needs to be a positive integer value.",
                      call. = FALSE)
              }
              # get samples with less than min number of reads
              if(min(colSums2(assay(x, assay.type))) < min_size){
                  rmsams <- colnames(x)[colSums2(assay(x, assay.type)) < min_size]
                  # Return NULL, if no samples were found after subsampling
                  if( !any(!colnames(x) %in% rmsams) ){
                      stop("No samples were found after subsampling.",
                          call. = FALSE)
                  }
                  if(verbose){
                      message(length(rmsams), " samples removed ",
                              "because they contained fewer reads than `min_size`.")
                  }
                  # remove sample(s)
                  newtse <- x[, !colnames(x) %in% rmsams]
              } else {
                  newtse <- x
              }
              newassay <- apply(assay(newtse, assay.type), 2, 
                                .subsample_assay,
                                min_size=min_size, replace=replace)
              rownames(newassay) <- rownames(newtse)
              # remove features not present in any samples after subsampling
              message(paste(length(which(rowSums2(newassay) == 0)), "features", 
                            "removed because they are not present in all samples", 
                            "after subsampling."))
              newassay <- newassay[rowSums2(newassay)>0,]
              newtse <- newtse[rownames(newassay),]
              assay(newtse, name, withDimnames=FALSE) <- newassay
              newtse <- .add_values_to_metadata(newtse, 
                                                "subsampleCounts_min_size",
                                                min_size)
              return(newtse)
          }
)

## Modified Sub sampling function from phyloseq internals
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
    suppressWarnings(subsample <- sample(obsvec,
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
