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
        if( !is.numeric(p_adj_threshold) && p_adj_threshold<0 || p_adj_threshold>1 ){
          stop("'p_adj_threshold' must be a numeric value >= 0.", call. = FALSE)
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
        results <- .calculate_correlation(assay1, assay2, method, p_adj_method)
        correlations <- results$correlations
        p_values <- results$p_values
        p_values_adjusted <- results$p_values_adjusted
        # Do filtering
        result <- .correlation_filter(correlations, 
                                      p_values, 
                                      p_values_adjusted, 
                                      p_adj_threshold,
                                      cor_threshold,
                                      sort,
                                      assay1, assay2, filter_self_correlations)
        # Matrix or table?
        if (mode == "table") {
          result <- .correlation_matrix_to_table(result)
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
  # Create empty matrices
  correlations <- matrix(NA, ncol(assay1), ncol(assay2))
  rownames(correlations) <- colnames(assay1)
  colnames(correlations) <- colnames(assay2)
  p_values <- correlations
  
  # Calculate correlations, different methods for numeric and categorical data
  if (method %in% c("kendall", "pearson","spearman")) {
    # Loop over every feature in assay2. Result is a list (feature1) of lists 
    # (correlations and p_values: individual feature1 vs all the feature2)
    correlations_and_p_values <- apply(assay2, 2, function(yi) {
      # Loop over every feature in assay1
      temp <- apply(assay1, 2, function(xi) {
        # Do correlation test, and store the result
        # to temporary object
          temp2 <- cor.test(xi, yi, 
                          method=method, use="pairwise.complete.obs")
          # Take only correlation and p-value
          temp2 <- c(temp2$estimate, temp2$p.value)
      })
      # Return a list where 1st element includes all the correlation values, and
      # second all thep-values
      list(temp[1,], temp[2,])
    })
    # Store correct names
    feature_names <- names(correlations_and_p_values)
    # Unlist list of lists to list
    correlations_and_p_values <- unlist(correlations_and_p_values, recursive = FALSE)
    # Take only correlations
    correlations <- correlations_and_p_values[seq(1, length(correlations_and_p_values), 2)]
    # Take only p-values
    p_values <- correlations_and_p_values[seq(2, length(correlations_and_p_values), 2)]
    
    # Unlisting changed names because otherwise there would have been duplicated names.
    # Give correct names back
    names(correlations) <- feature_names
    names(p_values) <- feature_names
    
    # Convert first to data frame and then to matrix
    correlations <- as.matrix(as.data.frame(correlations, check.names = FALSE))
    p_values <- as.matrix(as.data.frame(p_values, check.names = FALSE))
    
  } 
  # If method is categorical
  else if (method == "categorical") {
    
    correlations <- apply(assay2, 2, function(yi) {
      # Loop over every feature in assay1
      temp <- apply(assay1, 2, function(xi) {
        
        # Keep only those samples that have values in both features
        keep <- rowSums(is.na(cbind(xi, yi))) == 0
        xi <- xi[keep]
        yi <- yi[keep]
        
        # Calculate cross-correlation using Goorma and Kruskal tau
        .calculate_gktau(xi, yi) 
      })
    })
    p_values <- NULL
    p_values_adjusted <- NULL
  }
  
  # If there are p_values that are not NA, adjust them
  if ( !is.null(p_values) ) {
    # Corrected p-values
    p_values_adjusted <- matrix(p.adjust(p_values, method=p_adj_method), nrow=nrow(p_values))
    dimnames(p_values_adjusted) <- dimnames(p_values)
  } else{
    p_values_adjusted <- NULL
  }
  
  return(list(correlations = correlations, 
              p_values = p_values, 
              p_values_adjusted = p_values_adjusted))
}

############################## .correlation_filter #############################
.correlation_filter <- function(correlations,
                                p_values,
                                p_values_adjusted,
                                p_adj_threshold,
                                cor_threshold, 
                                sort,
                                assay1, assay2, filter_self_correlations){
  if( is.null(p_values_adjusted) ){
    p_adj_threshold <- NULL
  }
  
  # Filter
  if (!is.null(p_adj_threshold) || !is.null(cor_threshold)) {
    # Filter by adjusted p-values and correlations
    features1_p_value <- features1_p_value <- features1_correlation <- features1_correlation <- NULL
    # Which features have significant correlations?
    if (!is.null(p_adj_threshold) ) {
      p_adj_under_th <- p_values_adjusted < p_adj_threshold
      features1_p_value <- rowSums(p_adj_under_th, na.rm = TRUE) > 0
      features2_p_value <- colSums(p_adj_under_th, na.rm = TRUE) > 0
    }
    # Which features have correlation over correlation threshold?
    if (!is.null(cor_threshold)) {
      corr_over_th <- abs(correlations) > cor_threshold | abs(correlations) < cor_threshold
      features1_correlation <- rowSums(corr_over_th, na.rm = TRUE) > 0
      features2_correlation <- colSums(corr_over_th, na.rm = TRUE) > 0
    }
    # Combine results from previous steps
    if (!is.null(p_adj_threshold) && !is.null(cor_threshold)) {
      features1 <- features1_p_value & features1_correlation
      features2 <- features2_p_value & features2_correlation
    } else if (is.null(p_adj_threshold) && !is.null(cor_threshold)) {
      features1 <- features1_correlation
      features2 <- features2_correlation
    } else if (!is.null(p_adj_threshold) && is.null(cor_threshold)) {
      features1 <- features1_p_value
      features2 <- features2_p_value
    }
    
    # If both features have significant correlations
    if (sum(features1) > 0 && sum(features2) > 0) {
      # Get names of those features that were TRUE
      rownames <- rownames(correlations)[features1]
      colnames <- colnames(correlations)[features2]
      # Subset correlations, p_values, and adjusted p_values table
      correlations <- matrix(correlations[features1, features2, drop=FALSE], nrow=sum(features1))
      
      if(!is.null(p_values) && !is.null(p_values_adjusted) ){
        p_values <- matrix(p_values[features1, features2, drop=FALSE], nrow=sum(features1))
        p_values_adjusted <- matrix(p_values_adjusted[features1, features2, drop=FALSE], nrow=sum(features1))
        
        # Add row and column names
        rownames(p_values_adjusted) <- rownames(p_values) <- rownames(correlations) <- rownames
        colnames(p_values_adjusted) <- colnames(p_values) <- colnames(correlations) <- colnames
      } else{
        # Add row and column names
        rownames(correlations) <- rownames
        colnames(correlations) <- colnames
      }
      
      
      
      # If sort was specified and there is more than 1 feature left in both feature sets
      if (sort && sum(features1) >= 2 && sum(features2) >= 2) {
        # Order in visually appealing order
        tmp <- correlations
        rownames(tmp) <- NULL
        colnames(tmp) <- NULL
        # Do hierarchical clustering
        row_index <- hclust(as.dist(1 - cor(t(tmp),
                                       use="pairwise.complete.obs")))$order
        col_index <- hclust(as.dist(1 - cor(tmp,
                                       use="pairwise.complete.obs")))$order
        # Get the order of features from hiearchical clustering
        rownames <- rownames(correlations)[row_index]
        colnames <- colnames(correlations)[col_index]
        
        # Order the tables based on order of hierarchical clustering
        correlations <- correlations[row_index, col_index]
        
        if(!is.null(p_values) && !is.null(p_values_adjusted) ){
          p_values <- p_values[row_index, col_index]
          p_values_adjusted <- p_values_adjusted[row_index, col_index]
          
          # Add column and rownames
          rownames(p_values_adjusted) <- rownames(p_values) <- rownames(correlations) <- rownames
          colnames(p_values_adjusted) <- colnames(p_values) <- colnames(correlations) <- colnames
        } else{
          # Add column and rownames
          rownames(correlations) <- rownames
          colnames(correlations) <- colnames
        }
      }
    } 
    # If there were no significant correlations, give a message
    else {
      message("No significant correlations with the given criteria\n")
      correlations <- p_values <- p_values_adjusted <- NULL
    }
  }
  res <- list(cor=correlations, pval=p_values, p.adj=p_values_adjusted)
  # Filter self correlations if it's specified
  if (nrow(assay1) == nrow(assay2) && ncol(assay1) == ncol(assay2) && filter_self_correlations) {
    diag(res$cor) <- diag(res$pval) <- diag(res$p.adj) <- NA
  }
  return(res)
}

######################### .correlation_matrix_to_table #########################
.correlation_matrix_to_table <- function(res) {
  # Melt correlation table
  ctab <- ID <- NULL
  if (!is.null(res$cor)) {
    ctab <- as.data.frame(res$cor)
    ctab$ID <- rownames(res$cor)
    ctab <- reshape2::melt(ctab, "ID")
    colnames(ctab) <- c("X1", "X2", "Correlation")
    ctab$Correlation <- as.numeric(as.character(ctab$Correlation))
  
    # Melt p-values and add them to the melted correlation table
    # If there are adjusted p_values
    if (!is.null(res$p.adj)) {
      ctab2 <- as.data.frame(res$p.adj)
      ctab2$ID <- rownames(res$p.adj)
      ctab2 <- reshape2::melt(ctab2, "ID")
      colnames(ctab2) <- c("X1", "X2", "p.adj")
      ctab2$p.adj <- as.numeric(as.character(ctab2$p.adj))
      
      ctab <- cbind(ctab, ctab2$p.adj)
      colnames(ctab) <- c("X1", "X2", "Correlation", "p.adj")
      ctab <- ctab[order(ctab$p.adj), ]
      colnames(ctab) <- c("X1", "X2", "Correlation", "p.adj")
      
    } 
    # If there are only p-values that are not adjusted
    else {
      ctab2 <- as.data.frame(res$pval)
      ctab2$ID <- rownames(res$pval)
      ctab2 <- reshape2::melt(ctab2, "ID")
      colnames(ctab2) <- c("X1", "X2", "value")
      ctab2$value <- as.numeric(as.character(ctab2$value))
      
      ctab <- cbind(ctab, ctab2$value)
      ctab <- ctab[order(-abs(ctab$Correlation)), ]
      colnames(ctab) <- c("X1", "X2", "Correlation", "pvalue")
    }
  
    # Convert feature names into characters
    ctab$X1 <- as.character(ctab$X1)
    ctab$X2 <- as.character(ctab$X2)
    # Keep the original order of factor levels
    ctab$X1 <- factor(as.character(ctab$X1), levels=rownames(res$cor))
    ctab$X2 <- factor(as.character(ctab$X2), levels=colnames(res$cor))
    # Remove NAs
    ctab <- ctab[!is.na(ctab$Correlation), ]
    
    # Order the table by p-value
    if ("p.adj" %in% colnames(ctab)) {
      ctab <- ctab[order(ctab$p.adj), ]
    } else if ("pvalue" %in% colnames(ctab)) {
      ctab <- ctab[order(ctab$pvalue), ]
    }
  }
  ctab
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
