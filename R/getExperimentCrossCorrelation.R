#' Calculate cross-correlation
#' 
#' @param x A
#'   \code{\link[MultiAssayExperiment:MultiAssayExperiment-class]{MultiAssayExperiment}} or
#'   \code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#'   object.
#'   
#' @param experiment1 A single character or numeric value for selecting the experiment 1
#'    from \code{experiment(x)} of \code{MultiassayExperiment} object. 
#'    (By default: \code{experiment1 = 1})
#'    
#' @param experiment2 A single character or numeric value for selecting the experiment 2
#'    from\code{experiment(x)} of \code{MultiAssayExperiment} object or 
#'    \code{altExp(x)} of \code{SummarizedExperiment} object.
#'    (By default: \code{experiment2 = 2} when input is \code{MAE} and 
#'    \code{experiment2 = NULL} when input is \code{TreeSE})
#'    
#' @param abund_values1 A single character value for selecting the
#'   \code{\link[SummarizedExperiment:SummarizedExperiment-class]{assay}} of 
#'   experiment 1 to be transformed. (By default: \code{abund_values1 = "counts"})
#'   
#' @param abund_values2 A single character value for selecting the
#'   \code{\link[SummarizedExperiment:SummarizedExperiment-class]{assay}} of 
#'   experiment 2 to be transformed. (By default: \code{abund_values2 = "counts"})
#'   
#' @param method A single character value for selecting association method 
#'    ('pearson', or 'spearman' for continuous; categorical for discrete)
#'     (By default: \code{method = "spearman"})
#' 
#' @param mode A single character value for selecting output format 
#'    Must be 'table' or 'matrix'.  (By default: \code{mode = "table"})
#' 
#' @param p_adj_method A single character value for selecting adjustment method of
#'    p-values. Passed to \code{p.adjust} function. 
#'    Check available methods from \code{help(p.adjust))}.
#'    (By default: \code{p_adj_method = "fdr"})
#' 
#' @param p_adj_threshold A single numeric value (in [0, 1]) for selecting 
#'    adjusted p-value threshold. (By default: \code{p_adj_threshold = 0.05})
#' 
#' @param cor_threshold A single numeric absolute value (in [0, 1]) for selecting 
#'    correlation threshold to include features. (By default: \code{cor_threshold = NULL})
#' 
#' @param sort A single boolean value for selecting whether to sort features or not
#'    in result. Used method is hierarchical clustering. (By default: \code{sort = FALSE})
#' 
#' @param filter_self_correlations A single boolean value for selecting whether to 
#'    filter out correlations between identical items. Applies only when correlation
#'    between experiment itself is tested, i.e., when input is \code{MAE} and 
#'    \code{experiment1 == experiment2} or when input is \code{MAE} and 
#'    \code{experiment2 = NULL}.
#'    (By default: \code{filter_self_correlations = FALSE})
#' 
#' @param verbose Verbose
#'    
#' @param ... Additional arguments.
#'    
#'    
#' @details
#' Calculate cross-correlations between features of two experiments. 
#'
#' @references
#' Add references here. ###########################################################
#'
#' @return 
#' Matrix or table. ##############################################################
#'
#' @name getExperimentCrossCorrelation
#' @export
#'
#' @author Leo Lahti and Tuomas Borman. Contact: \url{microbiome.github.io}
#'
#' @examples
#' mae <- microbiomeDataSets::peerj32()
#' 
#' # Subset so that less observations / quicker to run, just for example
#' mae[[1]] <- mae[[1]][1:20, 1:10]
#' mae[[2]] <- mae[[2]][1:20, 1:10]
#' # Calculate cross-correlations
#' result <- getExperimentCrossCorrelation(mae, method = "pearson")
#' # Show first 5 entries
#' head(result, 5)
#' 
#' # Same can be done with SummarizedExperiment and altExp
#' # Create TreeSE with altExp
#' tse <- mae[[1]]
#' altExp(tse, "exp2") <- mae[[2]]
#' result <- getExperimentCrossCorrelation(tse, experiment2 = "exp2", method = "pearson")
#' # Show first 5 entries
#' head(result, 5)
#' 
NULL

#' @rdname getExperimentCrossCorrelation
#' @export
setGeneric("getExperimentCrossCorrelation", signature = c("x"),
           function(x,
                    experiment1 = 1,
                    experiment2 = 2,
                    abund_values1 = "counts",
                    abund_values2 = "counts",
                    method = c("categorical", "kendall", "pearson","spearman"),
                    mode = "table",
                    p_adj_method = c("BH", "bonferroni", "BY", "fdr", "hochberg", 
                                     "holm", "hommel", "none"),
                    p_adj_threshold = 0.05,
                    cor_threshold = NULL,
                    sort = FALSE,
                    filter_self_correlations = FALSE,
                    verbose = TRUE,
                    ...)
             standardGeneric("getExperimentCrossCorrelation"))

#' @rdname getExperimentCrossCorrelation
#' @export
setMethod("getExperimentCrossCorrelation", signature = c(x = "MultiAssayExperiment"),
    function(x,
             experiment1 = 1,
             experiment2 = 2,
             abund_values1 = "counts",
             abund_values2 = "counts",
             method = c("categorical", "kendall", "pearson","spearman"),
             mode = "table",
             p_adj_method = c("BH", "bonferroni", "BY", "fdr", "hochberg", 
                              "holm", "hommel", "none"),
             p_adj_threshold = 0.05,
             cor_threshold = NULL,
             sort = FALSE,
             filter_self_correlations = FALSE,
             ...){
        ############################# INPUT CHECK ##############################
        # Check experiment1 and experiment2
        # Negation of "if value is character and can be found from experiments or
        # if value is numeric and is smaller or equal to the list of experiments.
        if( !( is.character(experiment1) && experiment1 %in% names(experiments(x)) || 
            is.numeric(experiment1) && experiment1 <= length(experiments(x)) ) ){
          stop("'experiment1' must be numeric or character value specifying", 
               " experiment in experiment(x).", call. = FALSE)
        }
        if( !( is.character(experiment2) && experiment2 %in% names(experiments(x)) || 
            is.numeric(experiment2) && experiment2 <= length(experiments(x)) ) ){
          stop("'experiment2' must be numeric or character value specifying", 
               " experiment in experiment(x).", call. = FALSE)
        }
        # If experiments are not the same, disable filter_self_correlations
        if( experiment1 != experiment2 ){
          filter_self_correlations <- FALSE
        }
        # Fetch tse objects
        tse1 <- mae[[experiment1]]
        tse2 <- mae[[experiment2]]
        # Check abund_values1 and abund_values2
        .check_assay_present(abund_values1, tse1)
        .check_assay_present(abund_values2, tse2)
        # Check method
        # If method is not single string, user has not specified method,
        # or has given e.g. a vector
        if(!.is_non_empty_string(method)){
          stop("'method' must be a non-empty single character value.",
               call. = FALSE)
        }
        method <- match.arg(method)
        # Check mode
        if( !.is_non_empty_string(mode) && !mode %in% c("matrix", "table") ){
          stop("'mode' must be 'matrix' or 'table'.", call. = FALSE)
        }
        # Check p_adj_method
        p_adj_method <- match.arg(p_adj_method)
        # Check p_adj_threshold
        if( !(is.numeric(p_adj_threshold) && (p_adj_threshold>=0 && p_adj_threshold<=1)  || 
            is.null(cor_threshold) ) ){
          stop("'p_adj_threshold' must be a numeric value [0,1].", call. = FALSE)
        }
        # Check cor_threshold
        if( ifelse( is.numeric(cor_threshold), !(cor_threshold>=0 && cor_threshold<=1), !is.null(cor_threshold)) ){
          stop("'cor_threshold' must be a numeric value greater than or equal to 0.", 
               call. = FALSE)
        }
        # Check sort
        if( !(sort == TRUE || sort == FALSE) ){
          stop("'sort' must be a boolean value.", 
               call. = FALSE)
        }
        # Check filter_self_correlations
        if( !(filter_self_correlations == TRUE || filter_self_correlations == FALSE) ){
          stop("'filter_self_correlations' must be a boolean value.", 
               call. = FALSE)
        }
        ############################ INPUT CHECK END ###########################
        # Fetch assays to correlate
        assay1 <- assay(tse1, abund_values1)
        assay2 <- assay(tse2, abund_values1)
        # Transposes tables to right format
        assay1 <- t(assay1)
        assay2 <- t(assay2)

        # Test if data is in right format
        .cor_test_data_type(assay1, method)
        .cor_test_data_type(assay2, method)
        # Calculate correlations
        if(verbose){
          message("Calculating correlations...")
        }
        result <- .calculate_correlation(assay1, assay2, method, p_adj_method)
        # Do filtering
        if( !is.null(p_adj_threshold) || !is.null(cor_threshold) || filter_self_correlations ){
          if(verbose){
            message("Filtering results...")
          }
          result <- .correlation_filter(result, 
                                        p_adj_threshold,
                                        cor_threshold,
                                        sort,
                                        assay1, 
                                        assay2, 
                                        filter_self_correlations)
          
        }
        # Matrix or table?
        if (mode == "matrix") {
          if(verbose){
            message("Converting matrices into table...")
          }
          result <- .correlation_table_to_matrix(result)
          if(sort){
            if(verbose){
              message("Sorting results...")
            }
            result <- .correlation_sort(result, sort)
          }
        }
        return(result)
    }
)

#' @rdname getExperimentCrossCorrelation
#' @export
setMethod("getExperimentCrossCorrelation", signature = c(x = "SummarizedExperiment"),
    function(x,
             y = NULL,
             abund_values1 = "counts",
             abund_values2 = "counts",
             method = c("categorical", "kendall", "pearson","spearman"),
             mode = "table",
             p_adj_method = c("BH", "bonferroni", "BY", "fdr", "hochberg", 
                              "holm", "hommel", "none"),
             p_adj_threshold = 0.05,
             cor_threshold = NULL,
             sort = FALSE,
             filter_self_correlations = FALSE,
             ...){
        ############################## INPUT CHECK #############################
        if( !(class(y) == "SummarizedExperiment" || class(y) == "TreeSummarizedExperiment" ||
            (is.character(y) && y %in% names(altExps(x)))  || 
            (is.numeric(y) && y <= length(altExps(x))) ||
            is.null(y) ) ){
          stop("'y' must be SE or TreeSE object, or numeric or character value specifying", 
               " experiment in altExps(x) or it must be NULL.", call. = FALSE)
        }
        ############################ INPUT CHECK END ###########################
        # Fetch data sets and create a MAE object
        exp1 <- x
        # If experiment2 is NULL, then experiment1 == experiment2
        if( is.null(y) ){
          exp2 <- exp1
          x <- MultiAssayExperiment::MultiAssayExperiment(experiments = ExperimentList(exp1 = exp1))
          exp2_num <- 1
          
        } else if ( is.character(y) ){
          exp2 <- altExps(x)[[y]]
          x <- MultiAssayExperiment::MultiAssayExperiment(experiments = ExperimentList(exp1 = exp1, exp2 = exp2))
          exp2_num <- 2
        } else {
          exp2 <- y
          x <- MultiAssayExperiment::MultiAssayExperiment(experiments = ExperimentList(exp1 = exp1, exp2 = exp2))
          exp2_num <- 2
        }
        # Call method with MAE object as an input
        getExperimentCrossCorrelation(x,
                                      experiment1 = 1,
                                      experiment2 = exp2_num,
                                      abund_values1,
                                      abund_values2,
                                      method,
                                      mode,
                                      p_adj_method,
                                      p_adj_threshold,
                                      cor_threshold,
                                      sort,
                                      filter_self_correlations)
    }
)

################################ HELP FUNCTIONS ################################
############################## .cor_test_data_type #############################
.cor_test_data_type <- function(assay, method){
    # Different available methods
    numeric_methods <- c("kendall", "pearson","spearman")
    categorical_methods <- c("categorical")
    # Check if method match with values, otherwise give an error.
    # For numeric methods, expect only numeric values. For categorical methods, expect only factors.
    if (method %in% numeric_methods && !is.numeric(assay)) {
      # If there are no numeric values, give an error
      stop("Assay, specified by 'abund_values', of 'experiment1' does not include",
           " numeric values. Choose categorical method for 'method'.",
           call. = FALSE)
    } else if (method %in% categorical_methods && !is.character(assay)) {
      # If there are no factor values, give an error
      stop("Assay, specified by 'abund_values', of 'experiment1' does not include",
           " factor values. Choose numeric method for 'method'.",
           call. = FALSE)
    }
    return(assay)
}

############################# .calculate_correlation ###########################
.calculate_correlation <- function(assay1, assay2, method, p_adj_method){
    # Functions for correlation calculation
    FUN_numeric <- function(feature_pair){
      # Get features
      feature1 <- assay1[ , feature_pair[1]]
      feature2 <- assay2[ , feature_pair[2]]
      temp <- cor.test(feature1, feature2, 
                       method=method, use="pairwise.complete.obs")
      # Take only correlation and p-value
      temp <- c(temp$estimate, temp$p.value)
      return(temp)
    }
    FUN_categorical <- function(feature_pair){
      feature1 <- assay1[ , feature_pair[1]]
      feature2 <- assay2[ , feature_pair[2]]
      # Keep only those samples that have values in both features
      keep <- rowSums(is.na(cbind(feature1, feature2))) == 0
      feature1 <- feature1[keep]
      feature2 <- feature2[keep]
      # Calculate cross-correlation using Goorma and Kruskal tau
      .calculate_gktau(feature1, feature2)
    }
    # Calculate correlations, different methods for numeric and categorical data
    if (method %in% c("kendall", "pearson","spearman")) {
      FUN <- FUN_numeric
    } else {
      FUN <- FUN_categorical
    }
    
    # All the sample pairs
    feature_pairs <- expand.grid(colnames(assay1), colnames(assay2))
    # Calculate correlations
    correlations_and_p_values <- apply(feature_pairs, 1, FUN = FUN)
    # Transpose into the same orientation as feature-pairs
    correlations_and_p_values  <- t(correlations_and_p_values)
    # Give names
    if( ncol(correlations_and_p_values) == 1 ){
      colnames(correlations_and_p_values) <- c("cor")
    }
    else if( ncol(correlations_and_p_values) == 2 ){
      colnames(correlations_and_p_values) <- c("cor", "pval")
    }
    # Combine feature-pair names with correlation and p-values
    correlations_and_p_values <- cbind(feature_pairs, correlations_and_p_values)
    # If there are p_values that are not NA, adjust them
    if( !is.null(correlations_and_p_values$pval) ){
      correlations_and_p_values$p_adj <- p.adjust(correlations_and_p_values$pval,
                                                  method=p_adj_method)
    }
    
    return(correlations_and_p_values)
}

############################## .correlation_filter #############################
.correlation_filter <- function(result,
                                p_adj_threshold,
                                cor_threshold, 
                                sort,
                                assay1, assay2, filter_self_correlations){
    # If if there are no p_values disable p-value threshold
    if( is.null(result$p_adj) ){
      p_adj_threshold <- NULL
    }
    # Which features have significant correlations?
    if ( !is.null(result$p_adj) && !is.null(p_adj_threshold) ) {
      # Get those feature-pairs that have significant correlations
      p_adj_under_th <- result[result$p_adj < p_adj_threshold, ]
      # Subset so that there are only feature1s that have atleast one significant correlation
      result <- result[result[, "Var1"] %in% unique(p_adj_under_th[, "Var1"]), ]
      # Subset so that there are only feature2s that have atleast one significant correlation
      result <- result[result[, "Var2"] %in% unique(p_adj_under_th[, "Var2"]), ]
    }
    # Which features have correlation over correlation threshold?
    if ( !is.null(cor_threshold) ) {
      # Get those feature-pairs that have correlations over threshold
      corr_over_th <- result[abs(result$cor) > cor_threshold | abs(result$cor) < cor_threshold, ]
      # Subset so that there are only feature1s that have atleast one correlation over threshold
      result <- result[result[, "Var1"] %in% unique(corr_over_th[, "Var1"]), ]
      # Subset so that there are only feature2s that have atleast one correlation over threshold
      result <- result[result[, "Var2"] %in% unique(corr_over_th[, "Var2"]), ]
    }
    
    # If there are no significant correlations
    if ( nrow(result) == 0 ) {
      message("No significant correlations with the given criteria\n")
      result <- NULL
    } else{
      # Adjust levels
      result$Var1 <- factor(result$Var1)
      result$Var2 <- factor(result$Var2)

      # Filter self correlations if it's specified
      if ( identical(assay1, assay2) && filter_self_correlations ) {
        # Take only those rows where features differ
        result <- result[result$Var1 != result$Var2, ]
      }
      # Adjust rownames
      rownames(result) <- seq(nrow(result))
    }
    return(result)
}

.correlation_sort <- function(result, sort){
    # Fetch data
    correlations <- result$cor
    p_values <- result$pval
    p_values_adjusted <- result$p_adj
    
    # If sort was specified and there is more than 1 feature left in both feature sets
    if (sort && nrow(correlations) >= 2 && ncol(correlations) >= 2) {
      # Order in visually appealing order
      tmp <- correlations
      rownames(tmp) <- NULL
      colnames(tmp) <- NULL
      # Do hierarchical clustering
      row_index <- hclust(as.dist(1 - cor(t(tmp),
                                          use="pairwise.complete.obs")))$order
      col_index <- hclust(as.dist(1 - cor(tmp,
                                          use="pairwise.complete.obs")))$order
      # Get the order of features from hierarchical clustering
      rownames <- rownames(correlations)[row_index]
      colnames <- colnames(correlations)[col_index]
      
      # Order the correlation matrix  based on order of hierarchical clustering
      correlations <- correlations[row_index, col_index]
      # Add column and rownames
      rownames(correlations) <- rownames
      colnames(correlations) <- colnames
      
      # Oder also p-values if they are not NULL
      if(!is.null(p_values) && !is.null(p_values_adjusted) ){
        p_values <- p_values[row_index, col_index]
        p_values_adjusted <- p_values_adjusted[row_index, col_index]
        # Add column and rownames
        rownames(p_values_adjusted) <- rownames(p_values) <- rownames
        colnames(p_values_adjusted) <- colnames(p_values) <- colnames
      }
    }
    return(list(cor = correlations, 
                pval = p_values, 
                p_adj = p_values_adjusted))
}

######################### .correlation_table_to_matrix #########################
#' @importFrom reshape2 acast
.correlation_table_to_matrix <- function(result){
  cor <- reshape2::acast(result, Var1 ~ Var2, value.var = "cor")
  if( !is.null(result$pval) ){
    pval <- reshape2::acast(result, Var1 ~ Var2, value.var = "pval")
  } else{
    pval <- NULL
  }
  if( !is.null(result$p_adj) ){
    p_adj <- reshape2::acast(result, Var1 ~ Var2, value.var = "p_adj")
  } else{
    p_adj <- NULL
  }
  return(list(cor = cor, pval = pval, p_adj = p_adj))
}

############################### .calculate_gktau ###############################
.calculate_gktau <- function(x, y){
    # First, compute the IxJ contingency table between x and y
    Nij <- table(x, y, useNA="ifany")
    # Next, convert this table into a joint probability estimate
    PIij <- Nij/sum(Nij)
    # Compute the marginal probability estimates
    PIiPlus <- apply(PIij, MARGIN=1, sum)
    PIPlusj <- apply(PIij, MARGIN=2, sum)
    # Compute the marginal variation of y
    Vy <- 1 - sum(PIPlusj^2)
    # Compute the expected conditional variation of y given x
    InnerSum <- apply(PIij^2, MARGIN=1, sum)
    VyBarx <- 1 - sum(InnerSum/PIiPlus)
    # Compute and return Goodman and Kruskal's tau measure
    tau <- (Vy - VyBarx)/Vy
    tau
}
