#' Calculate correlations between features of two experiments.
#' 
#' @param x A
#'   \code{\link[MultiAssayExperiment:MultiAssayExperiment-class]{MultiAssayExperiment}} or
#'   \code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#'   object.
#'   
#' @param experiment1 A single character or numeric value for selecting the experiment 1
#'    from \code{experiments(x)} of \code{MultiassayExperiment} object. 
#'    (By default: \code{experiment1 = 1})
#'    
#' @param experiment2 A single character or numeric value for selecting the experiment 2
#'    from\code{experiments(x)} of \code{MultiAssayExperiment} object or 
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
#' @param direction A single character value for selecting if association are calculated
#'   row-wise (for features) or column-wise (for samples). Must be \code{"row"} or
#'   \code{"column"}. (By default: \code{direction = "row"}) 
#'   
#' @param method A single character value for selecting association method 
#'    ('kendall', pearson', or 'spearman' for continuous/numeric; 'categorical' for discrete)
#'     (By default: \code{method = "spearman"})
#' 
#' @param mode A single character value for selecting output format 
#'    Available formats are  'table' and 'matrix'.  (By default: \code{mode = "table"})
#' 
#' @param p_adj_method A single character value for selecting adjustment method of
#'    p-values. Passed to \code{p.adjust} function. 
#'    (By default: \code{p_adj_method = "fdr"})
#' 
#' @param p_adj_threshold A single numeric value (from 0 to  1) for selecting 
#'    adjusted p-value threshold for filtering. 
#'    (By default: \code{p_adj_threshold = NULL})
#' 
#' @param cor_threshold A single numeric absolute value (from 0 to 1) for selecting 
#'    correlation threshold for filtering.
#'    (By default: \code{cor_threshold = NULL})
#' 
#' @param sort A single boolean value for selecting whether to sort features or not
#'    in result matrices. Used method is hierarchical clustering. 
#'    Disabled when \code{mode = "table"}. (By default: \code{sort = FALSE})
#' 
#' @param filter_self_correlations A single boolean value for selecting whether to 
#'    filter out correlations between identical items. Applies only when correlation
#'    between experiment itself is tested, i.e., when assays are identical. 
#'    (By default: \code{filter_self_correlations = FALSE})
#' 
#' @param verbose A single boolean value for selecting whether to get messages
#'    about progress of calculation.
#'    
#' @param test_significance A single boolean value for selecting whether to test
#'    statistical significance of associations.
#'    
#' @param show_warnings A single boolean value for selecting whether to show warnings
#'    that might occur when correlations and p-values are calculated.
#'
#' @param paired A single boolean value for specifying if samples are paired or not.
#'    \code{colnames} must match between twp experiments. \code{paired} is disabled
#'    when \code{direction = "row"}. (By default: \code{paired = FALSE})
#'
#' @param ... Additional arguments:
#'    \itemize{
#'        \item{\code{test_significance}}{A single boolean value 
#'        in function \code{getExperimentCrossAssociation} for selecting 
#'        whether to test significance or not. 
#'        (By default: \code{test_significance = FALSE})}
#'        \item{\code{association_FUN}}{A function that is used to calculate (dis-)similarity
#'        between features. Function must take matrix as an input and give numeric
#'        values as an output. Adjust \code{method} and other parameters correspondingly.
#'        Supported functions are, for example, \code{stats::dist} and \code{vegan::vegdist}.}
#'    }
#'    
#' @details
#' These functions calculates associations between features of two experiments. 
#' \code{getExperimentCrossAssociation} calculates only associations by default.
#' \code{testExperimentCrossAssociation} calculates also significance of 
#' associations.
#'
#' @return 
#' These functions return associations in table or matrix format. In table format,
#' returned value is a data frame that includes  features and associations 
#' (and p-values) in columns. In matrix format, returned value is a one matrix
#' when only associations are calculated. If also significances are tested, then
#' returned value is a list of matrices.
#'
#' @name getExperimentCrossAssociation
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
#' result <- getExperimentCrossAssociation(mae, method = "pearson")
#' # Show first 5 entries
#' head(result, 5)
#' 
#' # Same can be done with SummarizedExperiment and altExp
#' # Create TreeSE with altExp
#' tse <- mae[[1]]
#' altExp(tse, "exp2") <- mae[[2]]
#' 
#' # When mode = matrix, matrix is returned
#' result <- getExperimentCrossAssociation(tse, experiment2 = "exp2", method = "pearson", 
#'                                         mode = "matrix")
#' # Show first 5 entries
#' head(result, 5)
#' 
#' # testExperimentCorrelation returns also significances
#' # filter_self_correlations = TRUE filters self correlations
#' # With p_adj_threshold it is possible to filter those features that do no have
#' # any correlations that have p-value under threshold
#' result <- testExperimentCrossAssociation(tse, experiment2 = tse, method = "pearson",
#'                                          filter_self_correlations = TRUE,
#'                                          p_adj_threshold = 0.05)
#' # Show first 5 entries
#' head(result, 5)
#' 
#' # Also getExperimentCrossAssociation returns significances when 
#' # test_signicance = TRUE
#' # Warnings can be suppressed by using show_warnings = FALSE
#' result <- getExperimentCrossAssociation(mae[[1]], experiment2 = mae[[2]], method = "pearson",
#'                                         mode = "matrix", test_significance = TRUE,
#'                                         show_warnings = FALSE)
#'                                         
#' # Returned value is a list of matrices
#' names(result)
#' 
#' # Calculate Bray-Curtis dissimilarity between features
#' result <- getExperimentCrossAssociation(mae[[1]], mae[[2]], 
#'                                         association_FUN = vegan::vegdist, method = "bray")
NULL

#' @rdname getExperimentCrossAssociation
#' @export
setGeneric("getExperimentCrossAssociation", signature = c("x"),
           function(x, ...)
               standardGeneric("getExperimentCrossAssociation"))

#' @rdname getExperimentCrossAssociation
#' @export
setMethod("getExperimentCrossAssociation", signature = c(x = "MultiAssayExperiment"),
    function(x,
           experiment1 = 1,
           experiment2 = 2,
           abund_values1 = "counts",
           abund_values2 = "counts",
           direction = "row",
           method = c("spearman", "categorical", "kendall", "pearson"),
           mode = "table",
           p_adj_method = c("fdr", "BH", "bonferroni", "BY", "hochberg", 
                            "holm", "hommel", "none"),
           p_adj_threshold = NULL,
           cor_threshold = NULL,
           sort = FALSE,
           filter_self_correlations = FALSE,
           verbose = TRUE,
           test_significance = FALSE,
           show_warnings = TRUE,
           paired = FALSE,
           ...){
        .get_experiment_cross_association(x,
                                          experiment1 = experiment1,
                                          experiment2 = experiment2,
                                          abund_values1 = abund_values1,
                                          abund_values2 = abund_values2,
                                          direction = direction,
                                          method = method,
                                          mode = mode,
                                          p_adj_method = p_adj_method,
                                          p_adj_threshold = p_adj_threshold,
                                          cor_threshold = cor_threshold,
                                          sort = sort,
                                          filter_self_correlations = filter_self_correlations,
                                          verbose = verbose,
                                          test_significance = test_significance,
                                          show_warnings = show_warnings,
                                          paired = paired,
                                          ...)
    }
)

#' @rdname getExperimentCrossAssociation
#' @importFrom MultiAssayExperiment MultiAssayExperiment ExperimentList
#' @importFrom SingleCellExperiment altExps
#' @export
setMethod("getExperimentCrossAssociation", signature = "SummarizedExperiment",
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
        .get_experiment_cross_association(x = x,
                                          experiment1 = 1,
                                          experiment2 = exp2_num,
                                          ...)
    }
)

#' @rdname getExperimentCrossAssociation
#' @export
setGeneric("testExperimentCrossAssociation", signature = c("x"),
           function(x, ...)
               standardGeneric("testExperimentCrossAssociation"))

#' @rdname getExperimentCrossAssociation
#' @export
setMethod("testExperimentCrossAssociation", signature = c(x = "ANY"),
          function(x, ...){
              getExperimentCrossAssociation(x, test_significance = TRUE, ...)
          }
)

############################## MAIN FUNCTIONALITY ##############################
# This function includes all the main functionality. 
.get_experiment_cross_association <- function(x,
                                              experiment1 = 1,
                                              experiment2 = 2,
                                              abund_values1 = "counts",
                                              abund_values2 = "counts",
                                              direction = "row",
                                              method = c("spearman", "categorical", "kendall", "pearson"),
                                              mode = "table",
                                              p_adj_method = c("fdr", "BH", "bonferroni", "BY", "hochberg", 
                                                               "holm", "hommel", "none"),
                                              p_adj_threshold = NULL,
                                              cor_threshold = NULL,
                                              sort = FALSE,
                                              filter_self_correlations = FALSE,
                                              verbose = TRUE,
                                              test_significance = FALSE,
                                              show_warnings = TRUE,
                                              paired = FALSE,
                                              ...){
    ############################# INPUT CHECK ##############################
    # Check experiment1 and experiment2
    .test_experiment_of_mae(x, experiment1)
    .test_experiment_of_mae(x, experiment2)
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
    # Check direction
    if( !.is_non_empty_string(direction) && !direction %in% c("row", "column") ){
      stop("'direction' must be 'row' or 'column'.", call. = FALSE)
    }
    # Check method
    # Checked in .calculate_association
    # Check mode
    if( !.is_non_empty_string(mode) && !mode %in% c("matrix", "table") ){
        stop("'mode' must be 'matrix' or 'table'.", call. = FALSE)
    }
    # Check p_adj_method
    p_adj_method <- match.arg(p_adj_method)
    # Check p_adj_threshold
    if( !(is.numeric(p_adj_threshold) && 
          (p_adj_threshold>=0 && p_adj_threshold<=1)  || 
          is.null(p_adj_threshold) ) ){
        stop("'p_adj_threshold' must be a numeric value [0,1].", call. = FALSE)
    }
    # Check cor_threshold
    if( !(is.numeric(cor_threshold) && 
          (cor_threshold>=0 && cor_threshold<=1)  || 
          is.null(cor_threshold) ) ){
      stop("'cor_threshold' must be a numeric value [0,1].", call. = FALSE)
    }
    # Check sort
    if( !.is_a_bool(sort) ){
        stop("'sort' must be a boolean value.", 
             call. = FALSE)
    }
    # Check filter_self_correlations
    if( !.is_a_bool(filter_self_correlations) ){
        stop("'filter_self_correlations' must be a boolean value.", 
             call. = FALSE)
    }
    # Check test_significance
    if( !.is_a_bool(test_significance) ){
        stop("'test_significance' must be a boolean value.", 
             call. = FALSE)
    }
    # Check verbose
    if( !.is_a_bool(verbose) ){
      stop("'verbose' must be a boolean value.", 
           call. = FALSE)
    }
    # Check show_warnings
    if( !.is_a_bool(show_warnings) ){
      stop("'show_warnings' must be a boolean value.", 
           call. = FALSE)
    }
    # Check paired
    if( !.is_a_bool(paired) ){
      stop("'paired' must be a boolean value.", 
           call. = FALSE)
    }
    ############################ INPUT CHECK END ###########################
    # Fetch assays to correlate
    assay1 <- assay(tse1, abund_values1)
    assay2 <- assay(tse2, abund_values2)
    # Transposes tables to right format, if row is specified
    if( direction == "row" ){
      assay1 <- t(assay1)
      assay2 <- t(assay2)
      # Disable paired
      paired <- FALSE
    }
    
    # If significance is not calculated, p_adj_method is NULL
    if( !test_significance ){
      p_adj_method <- NULL
    }
    # Calculate correlations
    result <- .calculate_association(assay1, assay2, method, p_adj_method, 
                                     test_significance, show_warnings, paired, 
                                     verbose, direction, ...)
    # Disable p_adj_threshold if there is no adjusted p-values
    if( is.null(result$p_adj) ){
      p_adj_threshold <- NULL
    }
    # Disable cor_threshold if there is no correlations
    if( is.null(result$cor) ){
      cor_threshold <- NULL
    }
    # Disable filter_self_correlation if assays are not the same
    if( !identical(assay1, assay2) ){
      filter_self_correlations <- FALSE
    }
    # Do filtering
    if( !is.null(p_adj_threshold) || 
        !is.null(cor_threshold) || 
        filter_self_correlations ){
        # Filter associations
        result <- .association_filter(result, 
                                      p_adj_threshold,
                                      cor_threshold,
                                      assay1, 
                                      assay2, 
                                      filter_self_correlations,
                                      verbose)
    }
    # Matrix or table?
    if (mode == "matrix" && !is.null(result) ) {
        # Create matrices from table
        result <- .association_table_to_matrix(result, verbose)
        
        # If matrix contains rows or columns that have only NAs, error occur in hclust
        if( any(rowSums(is.na(result$cor)) == ncol(result$cor)) || 
            any(colSums(is.na(result$cor)) == nrow(result$cor)) ){
            message("\nCorrelation matrices cannot be sorted, because correlation matrix ",
                    "contains rows and/or columns that contain only NAs.")
            sort <- FALSE
        }
        # If sort was specified and there are more than 1 features
        if(sort && nrow(result$cor) > 1 && ncol(result$cor) > 1 ){
            # Sort associations
            result <- .association_sort(result, verbose)
        }
    }
    # If result includes only one element, return only the element
    if( length(result) == 1 ){
        result <- result[[1]]
    }
    return(result)
}

################################ HELP FUNCTIONS ################################
# This function is for testing if experiment can be found from MAE
#' @importFrom MultiAssayExperiment experiments
.test_experiment_of_mae <- function(x, experiment){
    # If experiment is numeric and bigger than the number of experiments
    if( is.numeric(experiment) && experiment > length(experiments(x)) ){
        stop(paste0("'", deparse(substitute(experiment)), "' is greater than the",
                    " number of experiments in MAE object."), call. = FALSE)
    }
    # Negation of "if value is character and can be found from experiments or
    # if value is numeric and is smaller or equal to the list of experiments.
    if( !( is.character(experiment) && experiment %in% names(experiments(x)) || 
        is.numeric(experiment) && experiment <= length(experiments(x)) ) ){
        stop(paste0("'", deparse(substitute(experiment)), "'", 
                    " must be numeric or character value specifying", 
                    " experiment in experiment(x)."), call. = FALSE)
    }
}
####################### .cross_association_test_data_type ######################
# This function tests if values match with chosen method. With numeric methods, 
# numeric values are expected, and with categorical method factor or character 
# values are expected. Otherwise gives an error.

# Input: assay and method
# Output: assay
.cross_association_test_data_type <- function(assay, method){
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

########################### .check_if_paired_samples ###########################
# Check if samples are paired
.check_if_paired_samples <- function(assay1, assay2){
  if( !all(colnames(assay1) == colnames(assay2)) ){
    stop("Experiments are not paired or samples are in wrong order.", 
         "Check that colnames match between experiments.", call. = FALSE)
  }
}

############################# .calculate_association ###########################
# This function calculates correlations between every feature-pair. For numeric
# values uses cor.test() cor() and for categorical values Goodman and Kruskal's tau test.

# Input: Assays that share samples but that have different features and different parameters.
# Output: Correlation table including correlation values (and p-values and adjusted p-values)
#' @importFrom stats p.adjust
.calculate_association <- function(assay1, 
                                   assay2, 
                                   method = c("spearman", "categorical", "kendall", "pearson"), 
                                   p_adj_method, 
                                   test_significance, 
                                   show_warnings, 
                                   paired,
                                   verbose,
                                   direction,
                                   association_FUN = NULL, 
                                   ...){
    # Check method if association_FUN is not NULL
    if( is.null(association_FUN) ){
        method <- match.arg(method)
        # Get function name for message
        function_name <- ifelse(method == "categorical", "mia:::.calculate_gktau", 
                                ifelse(test_significance, "stats::cor.test", "stats::cor"))
        
        # Test if data is in right format
        .cross_association_test_data_type(assay1, method)
        .cross_association_test_data_type(assay2, method)
    } else{
        # Get name of function
        function_name <- deparse(substitute(association_FUN))
        test_significance <- FALSE
        p_adj_method <- NULL
    }
    
    # Message details of calculation
    if(verbose){
        message( paste0("Calculating correlations...\n",
                    "direction: ", direction, 
                    ", function: ", function_name, 
                    ", method: ", method,
                    ", test_significance: ", test_significance,
                    ", p_adj_method: ",
                    ifelse(!is.null(p_adj_method), p_adj_method, "-"),
                    ", paired: ", paired,
                    ", show_warnings: ", show_warnings) )
    }
  
    # If association_FUN is provided by user, use appropriate function.
    # Otherwise, choose correct method for numeric and categorical data
    if( !is.null(association_FUN) ){
        FUN_ <- .calculate_association_with_own_function
    } else if( method %in% c("kendall", "pearson","spearman") ) {
        FUN_ <- .calculate_association_for_numeric_values
    } else if( method == "categorical" ){
        FUN_ <- .calculate_association_for_categorical_values
    } 
    
    # Get all the sample/feature pairs
    if( paired ){
        .check_if_paired_samples(assay1, assay2)
        variable_pairs <- data.frame(Var1 = colnames(assay1), Var2 = colnames(assay2))
    } else{
        variable_pairs <- expand.grid(colnames(assay1), colnames(assay2))
    }
    
    # If function is stats::cor, then calculate associations directly with matrices
    # Otherwise, loop through variable_pairs
    if( function_name == "stats::cor" ){
        correlations_and_p_values <- .calculate_stats_cor(assay1 = assay1,
                                                          assay2 = assay2,
                                                          method = method,
                                                          variable_pairs = variable_pairs,
                                                          show_warnings = show_warnings)
    } else{
        correlations_and_p_values <- .calculate_association_table(variable_pairs = variable_pairs,
                                                                  FUN_ = FUN_, 
                                                                  test_significance = test_significance, 
                                                                  assay1 = assay1, 
                                                                  assay2 = assay2,
                                                                  method = method,
                                                                  show_warnings = show_warnings, 
                                                                  association_FUN = association_FUN, 
                                                                  ...)
    }
    
    # If there are p_values, adjust them
    if( !is.null(correlations_and_p_values$pval) ){
        correlations_and_p_values$p_adj <- 
            p.adjust(correlations_and_p_values$pval,
                     method = p_adj_method)
    }
    return(correlations_and_p_values)
}

######################### .sort_variable_pairs_row_wise ########################
# This function adds two columns to variable pair data.frame. This function
# sorts each variable pair in alphabetical order. First additional column includes
# first variable, and second second variable from variable pair.

# Input: variable pair data.frame
# Output: variable pair data.frame with two additional columns

#' @importFrom reshape2 melt dcast
.sort_variable_pairs_row_wise <- function(variable_pairs){
  
  variable_pairs_sorted <- variable_pairs
  # Create column that tells in which feature/sample pair variable is
  variable_pairs_sorted$ID <- 1:nrow(variable_pairs_sorted)
  # Convert into longer format so that one column have all variable names
  variable_pairs_sorted <- variable_pairs_sorted %>% reshape2::melt(id.vars = "ID")
  # Order based on ID (variable pair number) and value (alphabetical order of variable name)
  variable_pairs_sorted <- variable_pairs_sorted[order(variable_pairs_sorted$ID , variable_pairs_sorted$value), ]
  # First member of each variable pair is Var1 and second is Var2
  variable_pairs_sorted$variable <- paste0("Var", rep(c(1,2), times = nrow(variable_pairs_sorted)/2) )
  # Convert into wider format. Now table have own columns for each variable name and
  # variable names are in right order.
  variable_pairs_sorted  <- variable_pairs_sorted %>% reshape2::dcast(ID ~ variable, value.var = "value")
  # Take only columns that include variable names
  variable_pairs_sorted  <- variable_pairs_sorted[ , colnames(variable_pairs_sorted) %in%  c("Var1", "Var2") ]
  # Add new colnames
  colnames(variable_pairs_sorted) <- c("Var1_sorted", "Var2_sorted")
  # Add sorted feature pairs to original
  variable_pairs <- cbind(variable_pairs, variable_pairs_sorted)
  return(variable_pairs)
}
######################### .calculate_association_table #########################
# This function calculates association by looping through feature/sample-pairs.

# Input: feature/sample_pairs, and assays
# Output: correlation table with variable names in Var1 and Var2, and correlation 
# values in cor. Additionally, table can also include pval for p-values.
#' @importFrom dplyr filter left_join
.calculate_association_table <- function(variable_pairs,
                                         FUN_, 
                                         test_significance, 
                                         assay1, 
                                         assay2,
                                         method,
                                         show_warnings, 
                                         association_FUN, 
                                         ...){
    # Are assays identical?
    assays_identical <- identical(assay1, assay2)
    # If they are identical, we can make calculation faster
    # row1 vs col2 equals to row2 vs col1
    if( assays_identical ){
        # Unfactor/convert into characters so that filter works correctly
        variable_pairs$Var1 <- as.character(variable_pairs$Var1)
        variable_pairs$Var2 <- as.character(variable_pairs$Var2)
        
        variable_pairs_all <- variable_pairs
        # Sort features row-wise
        variable_pairs_all <- .sort_variable_pairs_row_wise(variable_pairs_all)
        # Get unique feature/sample pairs
        variable_pairs <- 
            variable_pairs_all[ !duplicated(variable_pairs_all[ , c("Var1_sorted", "Var2_sorted")]), ]
    }
    
    # Calculate correlations
    correlations_and_p_values <- apply(variable_pairs, 1, 
                                       FUN = FUN_, 
                                       test_significance = test_significance, 
                                       assay1 = assay1, 
                                       assay2 = assay2,
                                       method = method,
                                       show_warnings = show_warnings, 
                                       association_FUN = association_FUN, 
                                       ...)
  
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
    
    # If assays were identical, and duplicate variable pairs were dropped
    if( assays_identical ){
        # Change names so that they are not equal to colnames of variable_pairs
        colnames(variable_pairs)[1:2] <- c("Var1_", "Var2_")
        # Combine feature-pair names with correlation values and p-values
        correlations_and_p_values <- cbind(variable_pairs, correlations_and_p_values)
        
        # Combine two tables so that values are assigned to larger table with all the 
        # variable pairs
        correlations_and_p_values <- dplyr::left_join(variable_pairs_all, 
                                                      correlations_and_p_values, 
                                                      by = c("Var1_sorted", "Var2_sorted"))
        # Drop off additional columns
        correlations_and_p_values <- correlations_and_p_values[ !colnames(correlations_and_p_values) %in% 
                                                                  c("Var1_sorted", "Var2_sorted", 
                                                                    "Var1_", "Var2_")]
    } else{
        # Otherwise just add variable names
        correlations_and_p_values <- cbind(variable_pairs, correlations_and_p_values)
    }
    
    # Create factors from variable names
    correlations_and_p_values$Var1 <- factor(correlations_and_p_values$Var1)
    correlations_and_p_values$Var2 <- factor(correlations_and_p_values$Var2)
    
    return(correlations_and_p_values)
}

############################# .calculate_stats_cor #############################
# This function calculates correlations with stats::cor. Compared to other functions,
# it can take whole matrices as an input.

# Input: assays
# Output: correlation table
#' @importFrom stats cor 
#' @importFrom reshape2 melt
.calculate_stats_cor <- function(assay1, assay2, method, variable_pairs, show_warnings){
    # If user does not want warnings, 
    # suppress warnings that might occur when calculating correlaitons (NAs...)
    # or p-values (ties, and exact p-values cannot be calculated...)
    if( show_warnings ){
      correlations <- stats::cor(assay1, assay2, 
                                              method = method,
                                              use = "pairwise.complete.obs")
    } else{
      correlations <- suppressWarnings(stats::cor(assay1, assay2, 
                                                               method = method,
                                                               use = "pairwise.complete.obs") )
    }
    # melt matrices into long format, so that it match with output of other functions
    correlations <- reshape2::melt(correlations)
    # Adjust names
    colnames(correlations) <- c("Var1", "Var2", "cor")
    # Adjust order or drop pairs, e.g., if only paired samples were wanted
    correlations <- correlations[ correlations$Var1 == variable_pairs[, 1] & 
                                    correlations$Var2 == variable_pairs[, 2], ]
    return(correlations)
}

################### .calculate_association_for_numeric_values ##################
# This function calculates correlation between feature-pair. If significance test
# is specified, then it uses cor.test(), if not then cor().

# Input: Vector of names that belong to feature pair, and assays.
# Output: Correlation value or list that includes correlation value and p-value.
#' @importFrom stats cor.test
.calculate_association_for_numeric_values <- function(feature_pair, test_significance, 
                                                      assay1, assay2, method, show_warnings,
                                                      ...){
    # Get features
    feature1 <- assay1[ , feature_pair[1]]
    feature2 <- assay2[ , feature_pair[2]]

    # Calculate correlation
    # If user does not want warnings, 
    # suppress warnings that might occur when calculating correlations (NAs...)
    # or p-values (ties, and exact p-values cannot be calculated...)
    if( show_warnings ){
        temp <- cor.test(feature1,
                            feature2, 
                            method = method,
                            use = "pairwise.complete.obs")
    } else {
        temp <- suppressWarnings( cor.test(feature1,
                                            feature2, 
                                            method = method,
                                            use = "pairwise.complete.obs") )
    }
    
    # Take only correlation and p-value
    temp <- c(temp$estimate, temp$p.value)

    return(temp)
}

################# .calculate_association_for_categorical_values ################
# This function calculates correlation between feature-pair. Correlation value is 
# calculated with Goodman and Kruskal's tau test. P-value is calculated with Pearson's
# chi-squared test.

# Input: Vector of names that belong to feature pair, and assays.
# Output: Correlation value or list that includes correlation value and p-value.
.calculate_association_for_categorical_values <- 
    function(feature_pair,
             test_significance, 
             assay1,
             assay2,
             show_warnings,
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
                             test_significance = test_significance,
                             show_warnings)
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

################### .calculate_association_with_own_function  ##################
# This function calculates (dis-)similarity between feature-pair with function 
# that is specified by user.

# Input: Vector of names that belong to feature pair, and assays.
# Output: Correlation value or list that includes correlation value and p-value.
.calculate_association_with_own_function <- function(feature_pair,
                                         assay1, assay2, 
                                         association_FUN, 
                                         show_warnings,
                                         test_significance,
                                         ...){
  # Get features
  feature1 <- assay1[ , feature_pair[1]]
  feature2 <- assay2[ , feature_pair[2]]
  # Create a matrix 
  feature_mat <- rbind(feature1, feature2)
  
  # If user does not want warnings, 
  # suppress warnings that might occur when calculating correlations (NAs...)
  # or p-values (ties, and exact p-values cannot be calculated...)
  # Use try-catch to catch errors that might occur.
  if( show_warnings ){
      temp <- tryCatch({
          do.call(association_FUN, args = c(list(feature_mat), list(...)))
      },
      error = function(cond) {
          stop(paste0("Error occurred during calculation. Check that 
                      'association_FUN' fulfills requirements."), 
              call. = FALSE)
      }
      )
  } else {
      temp <- tryCatch({
          suppressWarnings( do.call(association_FUN, args = c(list(feature_mat), list(...))) )
      },
      error = function(cond) {
          stop(paste0("Error occurred during calculation. Check that ",
                      "'association_FUN' fulfills requirements."), 
              call. = FALSE)
      }
      )
  }
  
  # If temp's length is not 1, then function does not return single numeric value for each pair
  if( length(temp) != 1 ){
      stop(paste0("Error occurred during calculation. Check that ", 
                      "'association_FUN' fulfills requirements."), 
           call. = FALSE)
  } 
  return(temp)
}

############################## .association_filter #############################
# This filters off features that do not have any value under specified thresholds. 

# Input: Correlation table and thresholds
# Output: Filtered correlation table (or NULL if there are no observations after filtering)
.association_filter <- function(result,
                                p_adj_threshold,
                                cor_threshold, 
                                assay1,
                                assay2,
                                filter_self_correlations, 
                                verbose){
    # Give message if verbose == TRUE
    if(verbose){
      message( paste0("\nFiltering results...\np_adj_threshold: ",
                      ifelse(!is.null(p_adj_threshold),  p_adj_threshold, "-"), 
                      ", cor_threshold: ", 
                      ifelse(!is.null(cor_threshold), cor_threshold, "-"), 
                      ", filter_self_correlations: ", 
                      ifelse(filter_self_correlations,
                             filter_self_correlations, "-")) )
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
        message("\nNo significant correlations with the given criteria.\n")
        return(NULL)
    }
    # Adjust levels
    result$Var1 <- factor(result$Var1)
    result$Var2 <- factor(result$Var2)
    
    # Filter self correlations if it's specified
    if ( filter_self_correlations ) {
        # Take only those rows where features differ
        result <- result[result$Var1 != result$Var2, ]
    }
    # Adjust rownames, numbers from 1 to nrow()
    rownames(result) <- seq(nrow(result))
    return(result)
}

############################## .association_sort  ##############################
# This sorts correlation, p-value and adjusted p-values matrices based on 
# hierarchical clustering

# Input: List of matrices (cor, p-values and adjusted p-values / matrix (cor)
# Output: Lst of sorted matrices (cor, p-values and adjusted p-values / matrix (cor)
#' @importFrom stats hclust
.association_sort <- function(result, verbose){
    # Give message if verbose == TRUE
    if(verbose){
      message("\nSorting results...")
    }
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

######################### .association_table_to_matrix #########################
# This function converts correlation table into matrices

# Input: Correlation table
# Output: List of matrices (cor, p-values and adjusted p-values / matrix (cor)
#' @importFrom dplyr select
#' @importFrom tidyr pivot_wider
.association_table_to_matrix <- function(result, verbose){
    # Give message if verbose == TRUE
    if(verbose){
      message("\nConverting table into matrices...")
    }
    # Correlation matrix is done from Var1, Var2, and cor columns
    # Select correct columns
    cor <- result %>% dplyr::select("Var1", "Var2", "cor") %>% 
      # Create a tibble, colum names from Var2, values from cor,
      # first column includes Var1
      tidyr::pivot_wider(names_from = "Var2", values_from = "cor") %>%
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
        pval <- result %>% dplyr::select("Var1", "Var2", "pval") %>% 
          tidyr::pivot_wider(names_from = "Var2", values_from = "pval") %>% 
          as.data.frame()
        rownames(pval) <- pval$Var1
        pval$Var1 <- NULL
        pval <- as.matrix(pval)
        result_list[["pval"]] <- pval
    } 
    # If adjusted p_values exist, then create a matrix and add to the result list
    if( !is.null(result$p_adj) ){
        p_adj <- result %>% dplyr::select("Var1", "Var2", "p_adj") %>% 
          tidyr::pivot_wider(names_from = "Var2", values_from = "p_adj") %>% 
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
.calculate_gktau <- function(x, y, test_significance = FALSE, show_warnings){
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
    # Do the Pearson's chi-squared test.
    # If user does not want warnings,
    # suppress warnings that might occur when there are ties, and exact p-value
    # cant be calculated
    if( show_warnings ){
      temp <- chisq.test(x, y)
    } else {
      temp <- suppressWarnings( chisq.test(x, y) )
    }
    
    # Take the p-value
    p_value <- temp$p.value
    # Result is combination of tau and p-value
    return(list(estimate = tau, p.value = p_value))
}
