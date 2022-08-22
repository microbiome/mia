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
#'    \code{altExp(x)} of \code{TreeSummarizedExperiment} object. Alternatively, 
#'    \code{experiment2} can also be \code{TreeSE} object when \code{x} is \code{TreeSE} object.
#'    (By default: \code{experiment2 = 2} when \code{x} is \code{MAE} and 
#'    \code{experiment2 = x} when \code{x} is \code{TreeSE})
#'    
#' @param assay_name1 A single character value for selecting the
#'   \code{\link[SummarizedExperiment:SummarizedExperiment-class]{assay}} of 
#'   experiment 1 to be transformed. (By default: \code{assay_name1 = "counts"})
#'   
#' @param assay_name2 A single character value for selecting the
#'   \code{\link[SummarizedExperiment:SummarizedExperiment-class]{assay}} of 
#'   experiment 2 to be transformed. (By default: \code{assay_name2 = "counts"})
#'   
#' @param abund_values1 a single \code{character} value for specifying which
#'   assay of experiment 1 to use for calculation.
#'   (Please use \code{assay_name1} instead. At some point \code{abund_values1}
#'   will be disabled.)
#'   
#' @param abund_values2 a single \code{character} value for specifying which
#'   assay of experiment 2 to use for calculation.
#'   (Please use \code{assay_name2} instead. At some point \code{abund_values2}
#'   will be disabled.)
#' 
#' @param altExp1 A single numeric or character value specifying alternative experiment
#'   from the altExp of experiment 1. If NULL, then the experiment is itself 
#'   and altExp option is disabled. 
#'   (By default: \code{altExp1 = NULL})
#'   
#' @param altExp2 A single numeric or character value specifying alternative experiment
#'   from the altExp of experiment 2. If NULL, then the experiment is itself 
#'   and altExp option is disabled. 
#'   (By default: \code{altExp2 = NULL})
#'   
#' @param colData_variable1 A character value specifying column(s) from colData
#'   of experiment 1. If colData_variable1 is used, assay_name1 is disabled.
#'   (By default: \code{colData_variable1 = NULL})
#'   
#' @param colData_variable2 A character value specifying column(s) from colData
#'   of experiment 2. If colData_variable2 is used, assay_name2 is disabled.
#'   (By default: \code{colData_variable2 = NULL})
#' 
#' @param MARGIN A single numeric value for selecting if association are calculated
#'   row-wise / for features (1) or column-wise / for samples (2). Must be \code{1} or
#'   \code{2}. (By default: \code{MARGIN = 1}) 
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
#'    (By default: \code{sort = FALSE})
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
#'    when \code{MARGIN = 1}. (By default: \code{paired = FALSE})
#'
#' @param ... Additional arguments:
#'    \itemize{
#'        \item{\code{symmetric}}{ A single boolean value for specifying if 
#'        measure is symmetric or not. When \code{symmetric = TRUE}, associations
#'        are calculated only for unique variable-pairs, and they are assigned to 
#'        corresponding variable-pair. This decreases the number of calculations in 2-fold 
#'        meaning faster execution. (By default: \code{symmetric = FALSE}) }
#'        \item{\code{association_FUN}}{ A function that is used to calculate (dis-)similarity
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
#' data("HintikkaXOData")
#' mae <- HintikkaXOData
#' 
#' # Subset so that less observations / quicker to run, just for example
#' mae[[1]] <- mae[[1]][1:20, 1:10]
#' mae[[2]] <- mae[[2]][1:20, 1:10]
#' # Transform data
#' mae[[1]] <- transformSamples(mae[[1]], method = "rclr")
#' 
#' # Calculate cross-correlations
#' result <- getExperimentCrossAssociation(mae, method = "pearson", assay_name2 = "nmr")
#' # Show first 5 entries
#' head(result, 5)
#' 
#' # Use altExp option to specify alternative experiment from the experiment
#' altExp(mae[[1]], "Phylum") <- agglomerateByRank(mae[[1]], rank = "Phylum")
#' # Transform data
#' altExp(mae[[1]], "Phylum") <- transformSamples(altExp(mae[[1]], "Phylum"), method = "relabundance")
#' # When mode = matrix, matrix is returned
#' result <- getExperimentCrossAssociation(mae, experiment2 = 2, 
#'                                         assay_name1 = "relabundance", assay_name2 = "nmr",
#'                                         altExp1 = "Phylum", 
#'                                         method = "pearson", mode = "matrix")
#' # Show first 5 entries
#' head(result, 5)
#' 
#' # testExperimentCorrelation returns also significances
#' # filter_self_correlations = TRUE filters self correlations
#' # With p_adj_threshold it is possible to filter those features that do no have
#' # any correlations that have p-value under threshold
#' result <- testExperimentCrossAssociation(mae[[1]], experiment2 = mae[[1]], method = "pearson",
#'                                          filter_self_correlations = TRUE,
#'                                          p_adj_threshold = 0.05)
#' # Show first 5 entries
#' head(result, 5)
#' 
#' # Also getExperimentCrossAssociation returns significances when 
#' # test_signicance = TRUE
#' # Warnings can be suppressed by using show_warnings = FALSE
#' result <- getExperimentCrossAssociation(mae[[1]], experiment2 = mae[[2]], method = "pearson",
#'                                         assay_name2 = "nmr",
#'                                         mode = "matrix", test_significance = TRUE,
#'                                         show_warnings = FALSE)
#'                                         
#' # Returned value is a list of matrices
#' names(result)
#' 
#' # Calculate Bray-Curtis dissimilarity between samples. If dataset includes
#' # paired samples, you can use paired = TRUE.
#' result <- getExperimentCrossAssociation(mae[[1]], mae[[1]], MARGIN = 2, paired = FALSE,
#'                                         association_FUN = vegan::vegdist, method = "bray")
#'                                         
#' 
#' # If experiments are equal and measure is symmetric (e.g., taxa1 vs taxa2 == taxa2 vs taxa1),
#' # it is possible to speed-up calculations by calculating association only for unique
#' # variable-pairs. Use "symmetric" to choose whether to measure association for only
#' # other half of of variable-pairs.
#' result <- getExperimentCrossAssociation(mae, experiment1 = "microbiota", experiment2 = "microbiota", 
#'                                         assay_name1 = "counts", assay_name2 = "counts",
#'                                         symmetric = TRUE)
#' 
#' # For big data sets, calculation might take long. To make calculations quicker, you can take
#' # a random sample from data. In a complex biological problems, random sample
#' # can describe the data enough. Here our random sample is 30 % of whole data.
#' sample_size <- 0.3
#' tse <- mae[[1]]
#' tse_sub <- tse[ sample( seq_len( nrow(tse) ), sample_size * nrow(tse) ), ]
#' result <- testExperimentCrossAssociation(tse_sub)
#' 
#' # It is also possible to choose variables from colData and calculate association
#' # between assay and sample metadata or between variables of sample metadata
#' mae[[1]] <- estimateDiversity(mae[[1]])
#' # colData_variable works similarly to assay_name. Instead of fetching an assay
#' # named assay_name from assay slot, it fetches a column named colData_variable
#' # from colData.
#' result <- getExperimentCrossAssociation(mae[[1]], assay_name1 = "counts", 
#'                                         colData_variable2 = c("shannon", "coverage"))
#'                                         
NULL

#' @rdname getExperimentCrossAssociation
#' @aliases getExperimentCrossCorrelation
#' @export
setGeneric("getExperimentCrossAssociation", signature = c("x"),
           function(x, ...)
               standardGeneric("getExperimentCrossAssociation"))

#' @rdname getExperimentCrossAssociation
#' @aliases getExperimentCrossCorrelation
#' @export
setMethod("getExperimentCrossAssociation", signature = c(x = "MultiAssayExperiment"),
    function(x,
           experiment1 = 1,
           experiment2 = 2,
           assay_name1 = abund_values1, abund_values1 = "counts",
           assay_name2 = abund_values2, abund_values2 = "counts",
           altExp1 = NULL,
           altExp2 = NULL,
           colData_variable1 = NULL,
           colData_variable2 = NULL,
           MARGIN = 1,
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
                                          assay_name1 = assay_name1,
                                          assay_name2 = assay_name2,
                                          altExp1 = altExp1,
                                          altExp2 = altExp2,
                                          colData_variable1 = colData_variable1,
                                          colData_variable2 = colData_variable2,
                                          MARGIN = MARGIN,
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
#' @aliases getExperimentCrossCorrelation
#' @importFrom MultiAssayExperiment MultiAssayExperiment ExperimentList
#' @importFrom SingleCellExperiment altExps
#' @export
setMethod("getExperimentCrossAssociation", signature = "SummarizedExperiment",
    function(x, experiment2 = x, ...){
        ############################## INPUT CHECK #############################
        # If y is  SE or TreeSE object
        if( is(experiment2, "SummarizedExperiment") ){}
        # If y is  character specifying name of altExp, 
        else if( is.character(experiment2) && experiment2 %in% names(altExps(x)) ){}
        # If y is numeric value specifying altExp
        else if( is.numeric(experiment2) && experiment2 <= length(altExps(x)) ){} 
        # If y does not match, then give error
        else{
            stop("'experiment2' must be SE object, or numeric or character",
                 " value specifying experiment in altExps(x) or character value ",
                 " specifying column(s) from 'colData(x)' or it must be NULL.", 
                 call. = FALSE)
        }
        ############################ INPUT CHECK END ###########################
        # Fetch data sets and create a MAE object
        exp1 <- x
        args <- list()
        # If y is character or numeric, it specifies altExp
        if ( is.character(experiment2) || is.numeric(experiment2) ){
            exp2 <- altExps(x)[[experiment2]]
        } else {
            exp2 <- experiment2
        }
        # Create a MAE
        experiments <- ExperimentList(exp1 = exp1, exp2 = exp2)
        exp2_num <- 2
        x <- MultiAssayExperiment(experiments = experiments)
        # Call method with MAE object as an input
        .get_experiment_cross_association(x = x,
                                          experiment1 = 1,
                                          experiment2 = exp2_num,
                                          ...)
    }
)

#' @rdname getExperimentCrossAssociation
#' @aliases testExperimentCrossCorrelation
#' @export
setGeneric("testExperimentCrossAssociation", signature = c("x"),
           function(x, ...)
               standardGeneric("testExperimentCrossAssociation"))

#' @rdname getExperimentCrossAssociation
#' @aliases testExperimentCrossCorrelation
#' @export
setMethod("testExperimentCrossAssociation", signature = c(x = "ANY"),
          function(x, ...){
              getExperimentCrossAssociation(x, test_significance = TRUE, ...)
          }
)

############################# Methods for aliases ##############################
#' @rdname getExperimentCrossAssociation
#' @aliases testExperimentCrossAssociation
#' @export
setGeneric("testExperimentCrossCorrelation", signature = c("x"),
           function(x, ...)
               standardGeneric("testExperimentCrossCorrelation"))

#' @rdname getExperimentCrossAssociation
#' @aliases testExperimentCrossAssociation
#' @export
setMethod("testExperimentCrossCorrelation", signature = c(x = "ANY"),
          function(x, ...){
              getExperimentCrossAssociation(x, test_significance = TRUE, ...)
          }
)

#' @rdname getExperimentCrossAssociation
#' @aliases getExperimentCrossAssociation
#' @export
setGeneric("getExperimentCrossCorrelation", signature = c("x"),
           function(x, ...)
               standardGeneric("getExperimentCrossCorrelation"))

#' @rdname getExperimentCrossAssociation
#' @aliases getExperimentCrossAssociation
#' @export
setMethod("getExperimentCrossCorrelation", signature = c(x = "ANY"),
          function(x, ...){
              getExperimentCrossAssociation(x, ...)
          }
)
############################## MAIN FUNCTIONALITY ##############################
# This function includes all the main functionality.
#' @importFrom S4Vectors unfactor
.get_experiment_cross_association <- function(x,
                                              experiment1 = 1,
                                              experiment2 = 2,
                                              assay_name1 = "counts",
                                              assay_name2 = "counts",
                                              altExp1 = NULL,
                                              altExp2 = NULL,
                                              colData_variable1 = NULL,
                                              colData_variable2 = NULL,
                                              MARGIN = 1,
                                              method = c("spearman", "categorical", "kendall", "pearson"),
                                              mode = c("table", "matrix"),
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
    # Check and fetch tse objects
    tse1 <- .check_and_get_altExp(tse1, altExp1)
    tse2 <- .check_and_get_altExp(tse2, altExp2)
    # Check that experiments have same amount of samples
    if( ncol(tse1) != ncol(tse2) ){
        stop("Samples must match between experiments.",
             call. = FALSE)
    }
    # If variables from coldata are specified check them. Otherwise,
    # check assay_name1
    if( !is.null(colData_variable1) ){
        tse1 <- .check_and_subset_colData_variables(tse1, colData_variable1)
    } else{
        .check_assay_present(assay_name1, tse1)
    }
    if( !is.null(colData_variable2) ){
        tse2 <- .check_and_subset_colData_variables(tse2, colData_variable2)
    } else{
        .check_assay_present(assay_name2, tse2)
    }
    # Check MARGIN
    if( !is.numeric(MARGIN) && !MARGIN %in% c(1, 2) ){
      stop("'MARGIN' must be 1 or 2.", call. = FALSE)
    }
    # Check method
    # method is checked in .calculate_association. Otherwise association_FUN would
    # not work. (It can be "anything", and it might also have method parameter.)
    # Check mode
    mode <- match.arg(mode, c("table", "matrix"))
    p_adj_method <- match.arg(p_adj_method,
                              c("fdr", "BH", "bonferroni", "BY", "hochberg", 
                                "holm", "hommel", "none"))
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
    # Fetch assays to correlate, if variables from coldata are specified, take 
    # coldata, otherwise take assay
    if( !is.null(colData_variable1) ){
        assay1 <- colData(tse1)
        assay1 <- as.matrix(assay1)
        assay1 <- t(assay1)
    } else{
        assay1 <- assay(tse1, assay_name1)
    }
    if( !is.null(colData_variable2) ){
        assay2 <- colData(tse2)
        assay2 <- as.matrix(assay2)
        assay2 <- t(assay2)
    } else{
        assay2 <- assay(tse2, assay_name2)
    }
    
    # Transposes tables to right format, if row is specified
    if( MARGIN == 1 ){
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
    result <- .calculate_association(assay1, assay2, method, 
                                     p_adj_method, 
                                     test_significance, 
                                     show_warnings, paired, 
                                     verbose, MARGIN,
                                     assay_name1, assay_name2,
                                     altExp1, altExp2,
                                     colData_variable1, colData_variable2,
                                     ...)
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
    if( mode == "matrix" && !is.null(result) ) {
        # Create matrices from table
        result <- .association_table_to_matrix(result, verbose)
        # Sort associations if specified
        if( sort ){
            result <- .association_sort(result, verbose)
        }
        # Adjust names of all matrices
        result <- lapply(result, FUN = function(x){
            rownames(x) <- colnames(assay1)[ as.numeric(rownames(x)) ]
            colnames(x) <- colnames(assay2)[ as.numeric(colnames(x)) ]
            return(x)
        })
        # If result list includes only one matrix, return only the matrix
        if( length(result) == 1 ){
            result <- result[[1]]
        }
    } else if( !is.null(result) ){
        # Sort associations if specified
        if( sort ){
            result <- .association_sort(result, verbose)
            # Get levels
            levels1 <- unique( colnames(assay1)[ as.numeric(levels(result$Var1)) ] )
            levels2 <- unique( colnames(assay2)[ as.numeric(levels(result$Var2)) ] )
            # Unfactor so that factors do not affect when neames are adjusted
            result$Var1 <- unfactor(result$Var1)
            result$Var2 <- unfactor(result$Var2)
        }
        # Adjust names
        result$Var1 <- colnames(assay1)[ as.numeric(result$Var1) ]
        result$Var2 <- colnames(assay2)[ as.numeric(result$Var2) ]
        # Adjust factors
        if( sort ){
            # If the order is specified by sorting
            result$Var1 <- factor(result$Var1, levels = levels1)
            result$Var2 <- factor(result$Var2, levels = levels2)
        } else{
            # Otherwise the order is same as in input
            result$Var1 <- factor(result$Var1, levels = unique(result$Var1))
            result$Var2 <- factor(result$Var2, levels = unique(result$Var2))
        }
    }
    
    return(result)
}

################################ HELP FUNCTIONS ################################
# This function is for testing if experiment can be found from MAE
#' @importFrom MultiAssayExperiment experiments
.test_experiment_of_mae <- function(x, experiment){
    # If experiment is numeric and bigger than the number of experiments
    if( is.numeric(experiment) && experiment > length(experiments(x)) ){
        stop("'", deparse(substitute(experiment)), "' is greater than the",
             " number of experiments in MAE object.",
             call. = FALSE)
    }
    # Negation of "if value is character and can be found from experiments or
    # if value is numeric and is smaller or equal to the list of experiments.
    if( !( is.character(experiment) && experiment %in% names(experiments(x)) || 
        is.numeric(experiment) && experiment <= length(experiments(x)) ) ){
        stop("'", deparse(substitute(experiment)), "'", 
             " must be numeric or character value specifying", 
             " experiment in experiment(x).",
             call. = FALSE)
    }
    # Check experiment's class
    obj <- x[[experiment]]
    if( !(class(obj) == "TreeSummarizedExperiment" || 
          class(obj) == "SummarizedExperiment") ){
        stop("The class of experiment specified by ", 
             deparse(substitute(experiment)), " must be 'TreeSummarizedExperiment' ",
             "or 'SummarizedExperiment'.",
             call. = FALSE)
    }
}
############################# .check_and_get_altExp ############################
# This function checks if altExp is specified. If so, then it returns alternative
# experiment from altExp.

# Input: (Tree)SE
# Output: (Tree)SE
.check_and_get_altExp <- function(tse, altExp){
    # Get the variable names
    altExp_name <- deparse(substitute(altExp))
    exp_num <- substr(altExp_name, nchar(altExp_name), nchar(altExp_name))
    tse_name <- paste0("experiment ", exp_num)
    
    # If altExp is specified, check and get it. Otherwise return the original object
    if( !is.null(altExp) ){
        # Check altExp
        .check_altExp_present(altExp, tse, altExp_name, tse_name)
        # Get altExp and return it
        tse <- altExp(tse, altExp)
    }
    return(tse)
}
###################### .check_and_subset_colData_variables #####################
# This function checks if columns can be found from colData. Additionally, 
# integers are converted into numeric and factors to character.

# Input: (Tree)SE and character
# Output: (Tree)SE
.check_and_subset_colData_variables <- function(tse, variables){
    # Get variable name
    variable_name <- deparse(substitute(variables))
    var_num <- substr(variable_name, 
                      start = nchar(variable_name), stop = nchar(variable_name))
    # Check that variables can be found
    if( !(is.character(variables) &&
          all( variables %in% colnames(colData(tse))) ) ){
        stop("'", variable_name, "' must be a character value specifying ",
             "column(s) from colData of experiment ", var_num, ".",
             call. = FALSE)
    }
    # Get coldata
    coldata <- colData(tse)[ , variables, drop = FALSE]
    # Get the classes of variables
    classes <- unlist( lapply(coldata, class) )
    # IF there are factors, convert them into character
    if( any( classes == "factor") ){
        unfactor <- coldata[ , classes == "factor", drop = FALSE]
        unfactor <- lapply(unfactor, unfactor)
        coldata[ , classes == "factor"] <- as.data.frame(unfactor)
        classes[classes == "factor"] <- "character"
    }
    # IF there are integers, convert them into numeric
    int_or_double <- classes == "integer" | classes == "double"
    if( any( int_or_double ) ){
        as_numeric <- coldata[ , int_or_double, drop = FALSE]
        as_numeric <- lapply(as_numeric, as.numeric)
        coldata[ , int_or_double] <- as.data.frame(as_numeric)
        classes[int_or_double] <- "numeric"
    }
    # Check that all the variables have the same class
    if( length(unique(classes)) > 1 ){
        stop(" Variables specified by '", variable_name, "' do not share a same class.", 
             call. = FALSE)
    }
    # Replace the colData with new, subsetted colData
    coldata <- DataFrame(coldata)
    colData(tse) <- coldata
    return(tse)
}

####################### .cross_association_test_data_type ######################
# This function tests if values match with chosen method. With numeric methods, 
# numeric values are expected, and with categorical method factor or character 
# values are expected. Otherwise gives an error.

# Input: assay and method
# Output: assay
.cross_association_test_data_type <- function(assay, method, 
                                              colData_variable){
    # Different available methods
    numeric_methods <- c("kendall", "pearson","spearman")
    categorical_methods <- c("categorical")
    # Get message
    if( !is.null(colData_variable) ){
        message <- "Variables of colData"
    } else{
        message <- "Assay, specified by 'assay_name',"
    }
    # Check if method match with values, otherwise give an error.
    # For numeric methods, expect only numeric values. For categorical methods, 
    # expect only factors/characters.
    if (method %in% numeric_methods && !is.numeric(assay)) {
        # If there are no numeric values, give an error
        stop(message, " of 'experiment' does not ",
             "include numeric values. Choose categorical method for 'method'.",
             call. = FALSE)
    } else if (method %in% categorical_methods && !is.character(assay)) {
        # If there are no factor values, give an error
        stop(message, " specified by 'assay_name', of 'experiment' does not ",
             "include factor or character values. Choose numeric method for 'method'.",
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
                                   MARGIN,
                                   assay_name1, assay_name2,
                                   altExp1, altExp2,
                                   colData_variable1, colData_variable2,
                                   association_FUN = NULL,
                                   ...){
    # Check method if association_FUN is not NULL
    if( is.null(association_FUN) ){
        method <- match.arg(method)
        # Get function name for message
        function_name <- ifelse(method == "categorical", "mia:::.calculate_gktau", 
                                ifelse(test_significance, "stats::cor.test", "stats::cor"))
        
        # Test if data is in right format
        .cross_association_test_data_type(assay1, method, 
                                          colData_variable1)
        .cross_association_test_data_type(assay2, method, 
                                          colData_variable2)
    } else{
        # Get name of function
        function_name <- deparse(substitute(association_FUN))
        test_significance <- FALSE
        p_adj_method <- NULL
    }
    
    # Message details of calculation
    if(verbose){
        message( 
            "Calculating correlations...\n",
            "altExp1: ", ifelse(!is.null(altExp1), altExp1, "-"), 
            ", altExp2: ", ifelse(!is.null(altExp2), altExp2, "-"),
            ifelse(!is.null(colData_variable1), 
                paste0(", assay_name1: -, colData_variable1: ", 
                       paste(colData_variable1, collapse = " + ")), 
                paste0(", assay_name1: ", assay_name1, ", colData_variable1: -")),
            ifelse(!is.null(colData_variable2), 
                paste0(", assay_name2: -, colData_variable2: ", 
                       paste(colData_variable2, collapse = " + ")), 
                paste0(", assay_name2: ", assay_name2, ", colData_variable2: -")),
            "\nMARGIN: ", MARGIN, 
            ", function: ", function_name, 
            ", method: ", method,
            ", test_significance: ", test_significance,
            ", p_adj_method: ",
            ifelse(!is.null(p_adj_method), p_adj_method, "-"),
            ", paired: ", paired,
            ", show_warnings: ", show_warnings, "\n" 
            ) 
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
        variable_pairs <- data.frame( Var1 = seq_len(ncol(assay1)), Var2 = seq_len(ncol(assay2)) )
    } else{
        # Get feature_pairs as indices
        variable_pairs <- expand.grid( seq_len(ncol(assay1)), seq_len(ncol(assay2)) )
    }
    
    # If function is stats::cor, then calculate associations directly with matrices
    # Otherwise, loop through variable_pairs
    if( function_name == "stats::cor" ){
        correlations_and_p_values <- .calculate_stats_cor(assay1 = assay1,
                                                          assay2 = assay2,
                                                          method = method,
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
    
    # Get the order based on original order of variable-pairs
    order <- match( paste0(variable_pairs$Var1, "_", 
                           variable_pairs$Var2),
                    paste0(correlations_and_p_values$Var1, "_", 
                           correlations_and_p_values$Var2) )
    # Order the table
    correlations_and_p_values <- correlations_and_p_values[ order, ]
    # Remove rownames / make them to be in increasing order
    rownames(correlations_and_p_values) <- NULL
    
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
#' @importFrom DelayedMatrixStats rowMins rowMaxs
.sort_variable_pairs_row_wise <- function(variable_pairs){
    # Convert into matrix
    variable_pairs <- as.matrix(variable_pairs)
    # Get rowmins
    rowmins <- rowMins( variable_pairs )
    # Get rowmaxs
    rowmaxs <- rowMaxs( variable_pairs )
    # Convert back into data.frame
    variable_pairs <- as.data.frame( variable_pairs )
    # Add rowMins
    variable_pairs$Var1_sorted <- rowmins
    # Add rowMaxs
    variable_pairs$Var2_sorted <- rowmaxs
    return(variable_pairs)
}
######################### .calculate_association_table #########################
# This function calculates association by looping through feature/sample-pairs.

# Input: feature/sample_pairs, and assays
# Output: correlation table with variable names in Var1 and Var2, and correlation 
# values in cor. Additionally, table can also include pval for p-values.
.calculate_association_table <- function(variable_pairs,
                                         FUN_, 
                                         test_significance, 
                                         assay1, 
                                         assay2,
                                         method,
                                         show_warnings, 
                                         association_FUN, 
                                         symmetric = FALSE,
                                         ...){
    # Check symmetric
    if( !.is_a_bool(symmetric) ){
        stop("'symmetric' must be a boolean value.", 
             call. = FALSE)
    }
    # 
    # If user has specified that the measure is symmetric
    if( symmetric ){
        # Are assays identical? If so, calculate only unique variable-pairs
        assays_identical <- identical(assay1, assay2)
    } else{
        # Calculate all variable-pairs
        assays_identical <- FALSE
    }
    
    # If they are identical, we can make calculation faster
    # row1 vs col2 equals to row2 vs col1
    if( assays_identical ){
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
    # otherwise transpose into the same orientation as feature-pairs and then convert it to df
    if( is.vector(correlations_and_p_values) ){
      correlations_and_p_values <- data.frame(correlations_and_p_values)
    } else{
      correlations_and_p_values  <- t(correlations_and_p_values)
      correlations_and_p_values <- as.data.frame(correlations_and_p_values)
    }
    # Give names
    if( ncol(correlations_and_p_values) == 1 ){
      colnames(correlations_and_p_values) <- c("cor")
    }
    else if( ncol(correlations_and_p_values) == 2 ){
      colnames(correlations_and_p_values) <- c("cor", "pval")
    } else{
        stop("Unexpected error occurred during calculation.",
             call. = FALSE)
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
        correlations_and_p_values <- 
            correlations_and_p_values[ , !colnames(correlations_and_p_values) %in% 
                                           c("Var1_sorted", 
                                             "Var2_sorted", 
                                             "Var1_", 
                                             "Var2_") ]
        
    } else{
        # Otherwise just add variable names
        correlations_and_p_values <- cbind(variable_pairs, correlations_and_p_values)
    }
    
    return(correlations_and_p_values)
}

############################# .calculate_stats_cor #############################
# This function calculates correlations with stats::cor. Compared to other functions,
# it can take whole matrices as an input.

# Input: assays
# Output: correlation table
#' @importFrom stats cor 
#' @importFrom tidyr pivot_longer
.calculate_stats_cor <- function(assay1, assay2, method, show_warnings){
    # If user does not want warnings, 
    # suppress warnings that might occur when calculating correlations (NAs...)
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
    # 
    correlations <- as.data.frame(correlations)
    # Convert row names and column names into indices, so that they equal to 
    # other functions
    rownames(correlations) <- seq_len( nrow(correlations) )
    colnames(correlations) <- seq_len( ncol(correlations) )
    correlations$ID <- rownames(correlations)
    # melt matrix into long format, so that it match with output of other functions
    correlations <- correlations %>% tidyr::pivot_longer(cols = !"ID")
    
    # Adjust names
    colnames(correlations) <- c("Var1", "Var2", "cor")
    # Convert into data,frame
    correlations <- as.data.frame( correlations )
    # Convert indices back to numeric
    correlations$Var1 <- as.numeric(correlations$Var1)
    correlations$Var2 <- as.numeric(correlations$Var2)
    
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
            stop(paste0("Error occurred during calculation. Check, e.g., that ",
                "'association_FUN' fulfills requirements. 'association_FUN' ",
                "threw a following error:\n",  cond),
                call. = FALSE)
        })
    } else {
        temp <- tryCatch({
            suppressWarnings( do.call(association_FUN, args = c(list(feature_mat), list(...))) )
        },
        error = function(cond) {
            stop(paste0("Error occurred during calculation. Check, e.g., that ",
                    "'association_FUN' fulfills requirements. 'association_FUN' ",
                    "threw a following error:\n",  cond),
                 call. = FALSE)
        })
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
      message( "Filtering results...\np_adj_threshold: ",
                      ifelse(!is.null(p_adj_threshold),  p_adj_threshold, "-"), 
                      ", cor_threshold: ", 
                      ifelse(!is.null(cor_threshold), cor_threshold, "-"), 
                      ", filter_self_correlations: ", 
                      ifelse(filter_self_correlations,
                             filter_self_correlations, "-"), "\n" )
    }
    
    # Which features have significant correlations?
    if ( !is.null(result$p_adj) && !is.null(p_adj_threshold) ) {
        # Get those feature-pairs that have significant correlations
        result <- result[result$p_adj < p_adj_threshold & !is.na(result$p_adj), ]
    }
    # Which features have correlation over correlation threshold?
    if ( !is.null(cor_threshold) ) {
        # Get those feature-pairs that have correlations over threshold
        result <- result[abs(result$cor) > cor_threshold & !is.na(result$cor), ]
    }
    
    # If there are no significant correlations
    if ( nrow(result) == 0 ) {
        message("No significant correlations with the given criteria.\n")
        return(NULL)
    }
    
    # Filter self correlations if it's specified
    if ( filter_self_correlations ) {
        # Take only those rows where features differ
        result <- result[result$Var1 != result$Var2, ]
    }
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
      message("Sorting results...\n")
    }
    
    # Is the type of result table or matrix?
    is_dataframe <- is.data.frame(result)
    # If the type is data.frame, convert result first to matrix
    if( is_dataframe ){
        result_mat <- .association_table_to_matrix(result, verbose = FALSE)
        # Fetch data matrices
        correlations <- result_mat$cor
        p_values <- result_mat$pval
        p_values_adjusted <- result_mat$p_adj
    } else{
        # Fetch data matrices
        correlations <- result$cor
        p_values <- result$pval
        p_values_adjusted <- result$p_adj
    }
    
    # If matrix contains rows or columns that have only NAs, error occur in hclust
    if( (any(rowSums(is.na(correlations)) == ncol(correlations)) || 
         any(colSums(is.na(correlations)) == nrow(correlations))) ){
        message("Result cannot be sorted, because it ",
                "contains variable(s) whose correlation was not possible to calculate ",
                " since they all were NAs.\n")
        return(result)
    }
    # If there are only one variable, return as it is
    if( nrow(correlations) == 1 && ncol(result$cor) == 1 ){
        return(result)
    }
    
    # Order in visually appealing order
    tmp <- correlations
    rownames(tmp) <- NULL
    colnames(tmp) <- NULL
    
    # Do hierarchical clustering, use try-catch to catch errors that might occur
    # if data contains NAs
    row_index <- tryCatch({
        hclust(as.dist(1 - cor(t(tmp),
                               use="pairwise.complete.obs")))$order
    },
    error = function(cond) {
        stop(paste0("Error occurred during sorting. Possible reason is that ",
                    "correlation matrix includes NAs. Try with 'sort = FALSE'."), 
             call. = FALSE)
    }
    )
    col_index <- tryCatch({
        hclust(as.dist(1 - cor(tmp,
                               use="pairwise.complete.obs")))$order
    },
    error = function(cond) {
        stop(paste0("Error occurred during sorting. Possible reason is that ",
                    "correlation matrix includes NAs. Try with 'sort = FALSE'."), 
             call. = FALSE)
    }
    )
    # Get the order of features from hierarchical clustering
    rownames <- rownames(correlations)[row_index]
    colnames <- colnames(correlations)[col_index]
    
    # If the format of result was data.frame
    if( is_dataframe ){
        # Add order as factor levels
        result$Var1 <- factor(result$Var1, levels = row_index)
        result$Var2 <- factor(result$Var2, levels = col_index)
        
    } else{
        # Otherwise, the format was matrix
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
        result <- result_list
    }
    
    return(result)
}

######################### .association_table_to_matrix #########################
# This function converts correlation table into matrices

# Input: Correlation table
# Output: List of matrices (cor, p-values and adjusted p-values / matrix (cor)
#' @importFrom tidyr pivot_wider
#' @importFrom DelayedArray rowSums colSums
.association_table_to_matrix <- function(result, verbose){
    # Give message if verbose == TRUE
    if(verbose){
      message("Converting table into matrices...\n")
    }
    
    # Correlation matrix is done from Var1, Var2, and cor columns
    cor <- result %>%
        # Convert into long format
        tidyr::pivot_wider(id_cols = "Var1", names_from = "Var2", values_from = "cor") %>%
        # Convert into data.frame
        as.data.frame()
    # Give rownames and remove additional column
    rownames(cor) <- cor$Var1
    cor$Var1 <- NULL
    # Convert into matrix
    cor <- as.matrix(cor)
    # Remove empty rows and columns
    non_empty_rows <- rowSums(is.na(cor)) < ncol(cor)
    non_empty_cols <- colSums(is.na(cor)) < nrow(cor) 
    cor <- cor[ non_empty_rows, non_empty_cols, drop = FALSE ]
    # Create a result list
    result_list <- list(cor = cor)
    # If p_values exist, then create a matrix and add to the result list
    if( !is.null(result$pval) ){
        pval <- result %>%
            # Convert into long format
            tidyr::pivot_wider(id_cols = "Var1", names_from = "Var2", values_from = "pval") %>%
            # Convert into data.frame
            as.data.frame()
        # Adjust rownames and remove an additional column
        rownames(pval) <- pval$Var1
        pval$Var1 <- NULL
        # Convert into matrix
        pval <- as.matrix(pval)
        # Remove empty rows and columns
        pval <- pval[ non_empty_rows, non_empty_cols, drop = FALSE ]
        # Add it to result list
        result_list[["pval"]] <- pval
    } 
    # If adjusted p_values exist, then create a matrix and add to the result list
    if( !is.null(result$p_adj) ){
        p_adj <- result %>%
            # Convert into long format
            tidyr::pivot_wider(id_cols = "Var1", names_from = "Var2", values_from = "p_adj") %>%
            # Convert into data.frame
            as.data.frame()
        # Adjust rownames and remove an  additional column
        rownames(p_adj) <- p_adj$Var1
        p_adj$Var1 <- NULL
        # Convert into matrix
        p_adj <- as.matrix(p_adj)
        # Remove empty rows and columns
        p_adj <- p_adj[ non_empty_rows, non_empty_cols, drop = FALSE ]
        # Add it to result list
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
