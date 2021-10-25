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
#'    \code{altExp(x)} of \code{SummarizedExperiment} object. Alternatively, 
#'    \code{experiment2} can also be \code{SE} object when \code{x} is \code{SE} object.
#'    (By default: \code{experiment2 = 2} when \code{x} is \code{MAE} and 
#'    \code{experiment2 = x} when \code{x} is \code{SE})
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
#'    ('kendall', pearson', or 'spearman' for continuous; 'categorical' for discrete)
#'     (By default: \code{method = "spearman"})
#' 
#' @param mode A single character value for selecting output format 
#'    Must be 'table' or 'matrix'.  (By default: \code{mode = "table"})
#' 
#' @param p_adj_method A single character value for selecting adjustment method of
#'    p-values. Passed to \code{p.adjust} function. 
#'    (By default: \code{p_adj_method = "fdr"})
#' 
#' @param p_adj_threshold A single numeric value (in [0, 1]) for selecting 
#'    adjusted p-value threshold. (By default: \code{p_adj_threshold = 0.05})
#' 
#' @param cor_threshold A single numeric absolute value (in [0, 1]) for selecting 
#'    correlation threshold to include features. (By default: \code{cor_threshold = NULL})
#' 
#' @param sort A single boolean value for selecting whether to sort features or not
#'    in result matrices. Used method is hierarchical clustering. 
#'    Disabled when \code{mode = "table"}. (By default: \code{sort = FALSE})
#' 
#' @param filter_self_correlations A single boolean value for selecting whether to 
#'    filter out correlations between identical items. Applies only when correlation
#'    between experiment itself is tested, i.e., when input is \code{MAE} and 
#'    \code{experiment1 == experiment2} and \code{abund_values1 == abund_values2} 
#'    or when input is two identical \code{SE} objects 
#'    and \code{abund_values1 == abund_values2}. 
#'    (By default: \code{filter_self_correlations = FALSE})
#' 
#' @param verbose A single boolean value for selecting whether to get messages
#'    about progress of calculation.
#'
#' @param ... Additional arguments:
#'    \itemize{
#'        \item{\code{test_significance}}{A single boolean value 
#'        in function \code{getExperimentCrossCorrelation} for selecting 
#'        whether to test significance or not. 
#'        (By default: \code{test_significance = FALSE})}
#'    }
#'    
#'    
#' @details
#' These functions calculates associations between features of two experiments. 
#' \code{getExperimentCrossCorrelation} calculates only associations by default.
#' \code{testForExperimentCrossCorrelation} calculates also significance of 
#' associations.
#'
#' @return 
#' These functions return associations in table or matrix format. In table format,
#' returned value is a data frame that includes  features and associations 
#' (and p-values) in columns. In matrix format, returned value is a one matrix
#' when only associations are calculated. If also significances are tested, then
#' returned value is a list of matrices.
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
#' # Whe mode = matrix, matrix is returned
#' result <- getExperimentCrossCorrelation(tse, y = "exp2", method = "pearson", 
#'                                         mode = "matrix")
#' # Show first 5 entries
#' head(result, 5)
#' 
#' # testForExperimentCorrelation returns also significances
#' # filter_self_correlations = TRUE filters self correlations
#' result <- testForExperimentCrossCorrelation(tse, y = tse, method = "pearson",
#'                                             filter_self_correlations = TRUE)
#' # Show first 5 entries
#' head(result, 5)
#' 
#' # Also getExperimentCrossCorrelation returns significances when 
#' # test_signicance = TRUE
#' result <- getExperimentCrossCorrelation(mae[[1]], y = mae[[2]], method = "pearson",
#'                                         mode = "matrix", test_significance = TRUE)
#' # Returned value is a list of matrices
#' names(result)
NULL

#' @rdname getExperimentCrossCorrelation
#' @export
setGeneric("getExperimentCrossCorrelation", signature = c("x"),
           function(x, ...)
               standardGeneric("getExperimentCrossCorrelation"))

#' @rdname getExperimentCrossCorrelation
#' @export
setMethod("getExperimentCrossCorrelation", signature = c(x = "MultiAssayExperiment"),
    function(x,
           experiment1 = 1,
           experiment2 = 2,
           abund_values1 = "counts",
           abund_values2 = "counts",
           method = c("spearman", "categorical", "kendall", "pearson"),
           mode = "table",
           p_adj_method = c("fdr", "BH", "bonferroni", "BY", "hochberg", 
                            "holm", "hommel", "none"),
           p_adj_threshold = 0.05,
           cor_threshold = NULL,
           sort = FALSE,
           filter_self_correlations = FALSE,
           verbose = TRUE,
           ...){
        .get_experiment_cross_correlation(x,
                                          experiment1 = experiment1,
                                          experiment2 = experiment2,
                                          abund_values1 = abund_values1,
                                          abund_values2 = abund_values2,
                                          method = method,
                                          mode = mode,
                                          p_adj_method = p_adj_method,
                                          p_adj_threshold = p_adj_threshold,
                                          cor_threshold = cor_threshold,
                                          sort = sort,
                                          filter_self_correlations = filter_self_correlations,
                                          verbose = verbose,
                                          ...)
    }
)

#' @rdname getExperimentCrossCorrelation
#' @importFrom MultiAssayExperiment MultiAssayExperiment ExperimentList
#' @importFrom SingleCellExperiment altExps
#' @export
setMethod("getExperimentCrossCorrelation", signature = "SummarizedExperiment",
    function(x, experiment2 = x, ...){
        ############################## INPUT CHECK #############################
        # If y is  SE or TreeSE object
        if( class(experiment2) == "SummarizedExperiment" || 
            class(experiment2) == "TreeSummarizedExperiment" ){}
        # If y is  character specifying name of altExp, 
        else if( is.character(experiment2) && experiment2 %in% names(altExps(x)) ){}
        # If y is numeric value specifying altExp
        else if( is.numeric(experiment2) && experiment2 <= length(altExps(x)) ){} 
        # If y does not match, then give error
        else{
            stop("'y' must be SE or TreeSE object, or numeric or character value specifying", 
                " experiment in altExps(x) or it must be NULL.", call. = FALSE)
        }
        ############################ INPUT CHECK END ###########################
        # Fetch data sets and create a MAE object
        exp1 <- x
        # If experiment2 is character or numeric, it specifies altExp
        if ( is.character(experiment2) || is.numeric(experiment2) ){
            exp2 <- altExps(x)[[experiment2]]
            experiments <- ExperimentList(exp1 = exp1, exp2 = exp2)
            exp2_num <- 2
        } else {
            exp2 <- experiment2
            experiments <- ExperimentList(exp1 = exp1, exp2 = exp2)
            exp2_num <- 2
        }
        x <- MultiAssayExperiment(experiments = experiments)
        # Call method with MAE object as an input
        .get_experiment_cross_correlation(x = x,
                                          experiment1 = 1,
                                          experiment2 = exp2_num,
                                          ...)
    }
)

#' @rdname getExperimentCrossCorrelation
#' @export
setGeneric("testForExperimentCrossCorrelation", signature = c("x"),
           function(x, ...)
               standardGeneric("testForExperimentCrossCorrelation"))

#' @rdname getExperimentCrossCorrelation
#' @export
setMethod("testForExperimentCrossCorrelation", signature = c(x = "ANY"),
          function(x, ...){
              getExperimentCrossCorrelation(x, test_significance = TRUE, ...)
          }
)

############################## MAIN FUNCTIONALITY ##############################
# This function includes all the main functionality. 
.get_experiment_cross_correlation <- function(x,
                                              experiment1 = 1,
                                              experiment2 = 2,
                                              abund_values1 = "counts",
                                              abund_values2 = "counts",
                                              method = c("spearman", "categorical", "kendall", "pearson"),
                                              mode = "table",
                                              p_adj_method = c("fdr", "BH", "bonferroni", "BY", "hochberg", 
                                                               "holm", "hommel", "none"),
                                              p_adj_threshold = 0.05,
                                              cor_threshold = NULL,
                                              sort = FALSE,
                                              filter_self_correlations = FALSE,
                                              verbose = TRUE,
                                              test_significance = FALSE,
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
    # Fetch tse objects
    tse1 <- x[[experiment1]]
    tse2 <- x[[experiment2]]
    # Check that experiments have same amount of samples
    if( ncol(tse1) != ncol(tse2) ){
        stop("Samples must match between experiments.",
             call. = FALSE)
    }
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
    if( !(is.numeric(p_adj_threshold) && 
          (p_adj_threshold>=0 && p_adj_threshold<=1)  || 
          is.null(cor_threshold) ) ){
        stop("'p_adj_threshold' must be a numeric value [0,1].", call. = FALSE)
    }
    # Check cor_threshold
    if(!is.numeric(cor_threshold) && !is.null(cor_threshold)){
        stop("'cor_threshold' must be a numeric value greater than or ",
             "equal to 0.", 
             call. = FALSE)
    } else if(cor_threshold >= 0 && cor_threshold <= 1 && !is.null(cor_threshold)){
        stop("'cor_threshold' must be a numeric value greater than or ",
             "equal to 0.", 
             call. = FALSE)
    }
    # Check sort
    if( !(sort == TRUE || sort == FALSE) ){
        stop("'sort' must be a boolean value.", 
             call. = FALSE)
    }
    # Check filter_self_correlations
    if( !(filter_self_correlations == TRUE || 
          filter_self_correlations == FALSE) ){
        stop("'filter_self_correlations' must be a boolean value.", 
             call. = FALSE)
    }
    # Check test_significance
    if( !(test_significance == TRUE || test_significance == FALSE) ){
        stop("'test_significance' must be a boolean value.", 
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
    result <- .calculate_correlation(assay1, assay2, method, p_adj_method, 
                                     test_significance)
    # Do filtering
    if( !is.null(p_adj_threshold) || 
        !is.null(cor_threshold) || 
        filter_self_correlations ){
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
    if (mode == "matrix" && !is.null(result) ) {
        if(verbose){
            message("Converting table into matrices...")
        }
        result <- .correlation_table_to_matrix(result)
        # If sort was specified and there are more than 1 features
        if(sort && nrow(result$cor) > 1 && ncol(result$cor) > 1 ){
            if(verbose){
                message("Sorting results...")
            }
            result <- .correlation_sort(result)
        }
    }
    # If result includes only one element, return only the element
    if( length(result) == 1 ){
        result <- result[[1]]
    }
    return(result)
}

################################ HELP FUNCTIONS ################################
############################## .cor_test_data_type #############################
# This function tests if values match with chosen method. With numeric methods, 
# numeric values are expected, and with categorical method factor or character 
# values are expected. Otherwise gives an error.

# Input: assay and method
# Output: assay
.cor_test_data_type <- function(assay, method){
    # Different available methods
    numeric_methods <- c("kendall", "pearson","spearman")
    categorical_methods <- c("categorical")
    # Check if method match with values, otherwise give an error.
    # For numeric methods, expect only numeric values. For categorical methods, 
    # expect only factors/characters.
    if (method %in% numeric_methods && !is.numeric(assay)) {
        # If there are no numeric values, give an error
        stop("Assay, specified by 'abund_values', of 'experiment' does not ",
             "include numeric values. Choose categorical method for 'method'.",
             call. = FALSE)
    } else if (method %in% categorical_methods && !is.character(assay)) {
        # If there are no factor values, give an error
        stop("Assay, specified by 'abund_values', of 'experiment' does not ",
             "include factor values. Choose numeric method for 'method'.",
             call. = FALSE)
    }
    return(assay)
}

############################# .calculate_correlation ###########################
# This function calculates correlations between every feature-pair. For numeric
# values uses cor.test() cor() and for categorical values Goodman and Kruskal's tau test.

# Input: Assays that share samples but that have different features and different parameters.
# Output: Correlation table including correlation values (and p-values and adjusted p-values)
#' @importFrom stats p.adjust
.calculate_correlation <- function(assay1, assay2, method, p_adj_method, test_significance){
    # Choose correct method for numeric and categorical data
    if( method %in% c("kendall", "pearson","spearman") ) {
        FUN <- .calculate_correlation_for_numeric_values
    } else if( method == "categorical" ){
        FUN <- .calculate_correlation_for_categorical_values
    }
    
    # Get all the sample pairs
    feature_pairs <- expand.grid(colnames(assay1), colnames(assay2))
    # Calculate correlations
    correlations_and_p_values <- apply(feature_pairs, 1, 
                                       FUN = FUN, 
                                       test_significance = test_significance, 
                                       assay1 = assay1, 
                                       assay2 = assay2,
                                       method = method)
    
    # Convert into data.frame if it is vector, 
    # otherwise transpose into the same orientation as feature-pairs if it's not a vector
    if( is.vector(correlations_and_p_values) ){
        correlations_and_p_values <- data.frame(correlations_and_p_values)
    } else{
        correlations_and_p_values  <- t(correlations_and_p_values)
    }
    # Give names
    if( ncol(correlations_and_p_values) == 1 ){
        colnames(correlations_and_p_values) <- c("cor")
    }
    else if( ncol(correlations_and_p_values) == 2 ){
        colnames(correlations_and_p_values) <- c("cor", "pval")
    }
    # Combine feature-pair names with correlation values and p-values
    correlations_and_p_values <- cbind(feature_pairs, correlations_and_p_values)
    # If there are p_values, adjust them
    if( !is.null(correlations_and_p_values$pval) ){
        correlations_and_p_values$p_adj <- 
            p.adjust(correlations_and_p_values$pval,
                     method = p_adj_method)
    }
    return(correlations_and_p_values)
}

################### .calculate_correlation_for_numeric_values ##################
# This function calculates correlation between feature-pair. If significance test
# is specified, then it uses cor.test(), if not then cor().

# Input: Vector of names that belong to feature pair, and assays.
# Output: Correlation value or list that includes correlation value and p-value.
#' @importFrom stats cor.test cor 
.calculate_correlation_for_numeric_values <- function(feature_pair, test_significance, 
                                                      assay1, assay2, method, ...){
    # Get features
    feature1 <- assay1[ , feature_pair[1]]
    feature2 <- assay2[ , feature_pair[2]]
    # Whether to test significance
    if( test_significance ){
        #calculate correlatiom
        temp <- cor.test(feature1,
                         feature2, 
                         method = method,
                         use = "pairwise.complete.obs")
        # Take only correlation and p-value
        temp <- c(temp$estimate, temp$p.value)
    } else{
        # Calculate only correlation value
        temp <- cor(feature1,
                    feature2, 
                    method = method,
                    use = "pairwise.complete.obs")
    }
    return(temp)
}

################# .calculate_correlation_for_categorical_values ################
# This function calculates correlation between feature-pair. Correlation value is 
# calculated with Goodman and Kruskal's tau test. P-value is calculated with Pearson's
# chi-squared test.

# Input: Vector of names that belong to feature pair, and assays.
# Output: Correlation value or list that includes correlation value and p-value.
.calculate_correlation_for_categorical_values <- 
    function(feature_pair,
             test_significance, 
             assay1,
             assay2,
             ...){
    # Get features
    feature1 <- assay1[ , feature_pair[1]]
    feature2 <- assay2[ , feature_pair[2]]
    # Keep only those samples that have values in both features
    keep <- rowSums2(is.na(cbind(feature1, feature2))) == 0
    feature1 <- feature1[keep]
    feature2 <- feature2[keep]
    # Calculate cross-correlation using Goodman and Kruskal tau
    temp <- .calculate_gktau(feature1,
                             feature2,
                             test_significance = test_significance)
    # Whether to test significance
    if( test_significance ){
        # Take correlation and p-value
        temp <- c(temp$estimate, temp$p.value)
    } else{
        # Take only correlation value
        temp <- temp$estimate
    }
    return(temp)
}

############################## .correlation_filter #############################
# This filters off features that do not have any value under specified thresholds. 

# Input: Correlation table and thresholds
# Output: Filtered correlation table (or NULL if there are no observations after filtering)
.correlation_filter <- function(result,
                                p_adj_threshold,
                                cor_threshold, 
                                sort,
                                assay1,
                                assay2,
                                filter_self_correlations){
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
        return(N)
    }
    # Adjust levels
    result$Var1 <- factor(result$Var1)
    result$Var2 <- factor(result$Var2)
    
    # Filter self correlations if it's specified and assays match with each other
    if ( identical(assay1, assay2) && filter_self_correlations ) {
        # Take only those rows where features differ
        result <- result[result$Var1 != result$Var2, ]
    }
    # Adjust rownames, numbers from 1 to nrow()
    rownames(result) <- seq(nrow(result))
    return(result)
}

############################## .correlation_sort  ##############################
# This sorts correlation, p-value and adjusted p-values matrices based on 
# hierarchical clustering

# Input: List of matrices (cor, p-values and adjusted p-values / matrix (cor)
# Output: Lst of sorted matrices (cor, p-values and adjusted p-values / matrix (cor)
#' @importFrom stats hclust
.correlation_sort <- function(result){
    # Fetch data matrices
    correlations <- result$cor
    p_values <- result$pval
    p_values_adjusted <- result$p_adj
    
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
    colnames(correlations) <- colnames
    rownames(correlations) <- rownames
    
    # Create a result list
    result_list <- list(cor = correlations)
    
    # If p_values is not NULL. order and add to the list
    if( !is.null(p_values) ){
        # Sort the matrix
        p_values <- p_values[row_index, col_index]
        # Add column and rownames
        colnames(p_values) <- colnames
        rownames(p_values) <- rownames
        # Add matrix to result list
        result_list[["pval"]] <- p_values
    }
    # If p_values_adjusted is not NULL. order and add to the list
    if( !is.null(p_values_adjusted) ){
        # Sort the matrix
        p_values_adjusted <- p_values_adjusted[row_index, col_index]
        # Add column and rownames
        colnames(p_values_adjusted) <- colnames
        rownames(p_values_adjusted) <- rownames
        # Add matrix to result list
        result_list[["p_adj"]] <- p_values_adjusted
    }
    return(result_list)
}

######################### .correlation_table_to_matrix #########################
# This function converts correlation table into matrices

# Input: Correlation table
# Output: List of matrices (cor, p-values and adjusted p-values / matrix (cor)
#' @importFrom dplyr select
#' @importFrom tidyr pivot_wider
.correlation_table_to_matrix <- function(result){
    # Correlation matrix is done from Var1, Var2, and cor columns
    # Select correct columns
    cor <- result %>% dplyr::select(Var1, Var2, cor) %>% 
      # Create a tibble, colum names from Var2, values from cor,
      # first column includes Var1
      tidyr::pivot_wider(names_from = Var2, values_from = cor) %>%
      # Convert into data frame
      as.data.frame()
    # Give rownames and remove additional column
    rownames(cor) <- cor$Var1
    cor$Var1 <- NULL
    cor <- as.matrix(cor)
    # Create a result list
    result_list <- list(cor = cor)
    # If p_values exist, then create a matrix and add to the result list
    if( !is.null(result$pval) ){
        pval <- result %>% dplyr::select(Var1, Var2, pval) %>% 
          tidyr::pivot_wider(names_from = Var2, values_from = pval) %>% 
          as.data.frame()
        rownames(pval) <- pval$Var1
        pval$Var1 <- NULL
        pval <- as.matrix(pval)
        result_list[["pval"]] <- pval
    } 
    # If adjusted p_values exist, then create a matrix and add to the result list
    if( !is.null(result$p_adj) ){
        p_adj <- result %>% dplyr::select(Var1, Var2, p_adj) %>% 
          tidyr::pivot_wider(names_from = Var2, values_from = p_adj) %>% 
          as.data.frame()
        rownames(p_adj) <- p_adj$Var1
        p_adj$Var1 <- NULL
        p_adj <- as.matrix(p_adj)
        result_list[["p_adj"]] <- p_adj
    } 
    return(result_list)
}

############################### .calculate_gktau ###############################
# This function calculates association with Goodman and Kruskal's tau test and 
# p-value with  Pearson's chi-squared test

# Input: Two vectors, one represent feature1, other feature2, which share samples
# Output: list of tau and p-value or just tau
#' @importFrom DelayedMatrixStats rowSums2 colSums2
#' @importFrom stats chisq.test
.calculate_gktau <- function(x, y, test_significance = FALSE){
    # First, compute the IxJ contingency table between x and y
    Nij <- table(x, y, useNA="ifany")
    # Next, convert this table into a joint probability estimate
    PIij <- Nij/sum(Nij)
    # Compute the marginal probability estimates
    #PIiPlus <- apply(PIij, MARGIN=1, sum)
    PIiPlus <- rowSums2(PIij)
    #PIPlusj <- apply(PIij, MARGIN=2, sum)
    PIPlusj <- colSums2(PIij)
    # Compute the marginal variation of y
    Vy <- 1 - sum(PIPlusj^2)
    # Compute the expected conditional variation of y given x
    # Because of the previous step, there should not be any NAs
    #InnerSum <- apply(PIij^2, MARGIN=1, sum)
    InnerSum <- rowSums2(PIij^2)
    VyBarx <- 1 - sum(InnerSum/PIiPlus)
    # Compute and return Goodman and Kruskal's tau measure
    tau <- (Vy - VyBarx)/Vy
    
    # If test significance is specified, then calculate significance with chi-squared test.
    # Are these two features independent or not?
    if ( !test_significance ){
        return(list(estimate = tau))
    } 
    # Do the Pearson's chi-squared test
    temp <- chisq.test(x, y)
    # Take the p-value
    p_value <- temp$p.value
    # Result is combination of tau and p-value
    return(list(estimate = tau, p.value = p_value))
}
