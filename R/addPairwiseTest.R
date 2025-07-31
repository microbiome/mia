#' @name getPairwiseTest
#' 
#' @title Get Pairwise Statistical Test
#' 
#' @description 
#' \code{getPairwiseTest()} performs pairwise comparisons between the levels of a
#' grouping variable using either Wilcoxon rank-sum, Student's \emph{t} test or
#' Kruskal–Wallis \emph{H} test (with optional pairwise follow-ups).  
#' It automatically pulls the response variable from assays, \code{rowData}, or
#' \code{colData} of a \link[TreeSummarizedExperiment]{TreeSummarizedExperiment}
#' (or base \link[SummarizedExperiment]{SummarizedExperiment}) object, returning
#' a tidy data frame of statistics and effect sizes.
#' 
#' @param x
#' \code{\link[TreeSummarizedExperiment:TreeSummarizedExperiment-class]{TreeSummarizedExperiment}}.
#'
#' @param assay.type \code{NULL} or \code{character scalar}. Specifies the 
#' assay to use for abundance values. Must be one of the available assays 
#' in \code{x}. (Default: \code{NULL})
#'
#' @param row.var \code{NULL} or \code{character scalar}. Specifies a variable 
#' from \code{rowData(x)} to test. (Default: \code{NULL})
#'
#' @param col.var \code{NULL} or \code{character scalar}. Specifies a variable 
#' from \code{colData(x)} to test. (Default: \code{NULL})
#'
#' @param group \code{character scalar}. Specifies a grouping variable, 
#' either from \code{rowData(x)} or \code{colData(x)} depending on the input mode.
#'
#' @param significance.method \code{character scalar}. Statistical method to 
#' use: one of \code{"wilcoxon"}, \code{"t.test"}, or \code{"kruskal"}. 
#' (Default: \code{"wilcoxon"})
#'
#' @param p.adjust.method \code{character scalar}. Method for p-value adjustment. 
#' Must be one of \code{p.adjust.methods}. (Default: \code{"fdr"})
#'
#' @param paired \code{logical scalar}. Whether to perform a paired test. 
#' (Default: \code{FALSE})
#'
#' @param return.type \code{character scalar}. Specifies the format of the result: 
#' one of \code{"full"}, \code{"summary"}, or \code{"both"}. (Default: \code{"full"})
#'
#' @param include.effect \code{logical scalar}. Whether to compute effect sizes. 
#' (Default: \code{TRUE})
#'
#' @param ... Additional arguments passed to the test function (e.g. 
#' \code{alternative = "greater"}).
#'
#' @section Effect Sizes:
#' Effect sizes are computed with \code{rstatix} helpers:
#' \itemize{
#'   \item Wilcoxon → \code{rstatix::wilcox_effsize()} (rank-biserial $r$)
#'   \item t-test  → \code{rstatix::cohens_d()} (Hedges-adjusted $g$)
#'   \item Kruskal → \code{rstatix::kruskal_effsize()} (eta-squared $\eta^{2}$)
#' }
#'
#' @return
#' \strong{If \code{method != "kruskal"}} a single
#' \code{tibble} with one row per comparison, containing:
#' \itemize{
#'   \item \code{group1}, \code{group2} – compared levels  
#'   \item \code{n1}, \code{n2} – sample sizes  
#'   \item \code{p}, \code{p.adj} – raw and adjusted P-values  
#'   \item \code{statistic}, \code{df} (for t-test)  
#'   \item \code{mean_group1}, \code{mean_group2}, \code{log2FC}  
#'   \item \code{effsize}, \code{magnitude} (if \code{include.effect = TRUE})
#' }
#'
#' \strong{If \code{method == "kruskal"}} a list with:
#' \itemize{
#'   \item \code{$global} – one-row tibble from \code{kruskal_test()}  
#'   \item \code{$pairwise} – pairwise Wilcoxon tests as above
#' }
#'
#' @examples
#' # Load example data
#' data(GlobalPatterns, package = "mia")
#' tse <- GlobalPatterns
#' 
#' # Transform to relative abundances
#' tse <- transformAssay(tse, method = "relabundance") 
#' 
#' # Pairwise test on assay data grouped by sample type
#' result <- getPairwiseTest(
#'     tse, 
#'     assay.type = "relabundance",
#'     group = "SampleType",
#'     method = "wilcoxon"
#' )
#' 
#' # Add results to metadata
#' tse <- addPairwiseTest(
#'     tse,
#'     assay.type = "relabundance", 
#'     group = "SampleType",
#'     name = "sample_type_comparison"
#' )
#' 
#' @seealso
#' \code{\link[rstatix]{pairwise_wilcox_test}},
#' \code{\link[rstatix]{pairwise_t_test}},
#' \code{\link[rstatix]{kruskal_test}}
#' 
NULL

#' @rdname getPairwiseTest
#' @export
setMethod("getPairwiseTest", signature(x = "SummarizedExperiment"), 
    function( x, assay.type = NULL, features = NULL, row.var = NULL, 
              col.var = NULL, group, facet.by = NULL, comp.by = NULL,
              pair.by = NULL,
              ...) {
    # Add rownames to rowData so that they are available for testing
    rowData(x)[["rownames"]] <- rownames(x)
    args <- c(list(
        tse = x, assay.type = assay.type, features = features,
        row.var = row.var, col.var = col.var, group = group, 
        facet.by = facet.by, comp.by = comp.by), pair.by = pair.by,
        list(...)
    )
    temp <- do.call(.check_input_for_pairwise_test, args)
    # Get the data from the object
    df <- do.call(.get_data_for_pairwise_test, args)

    res <- .calculate_pairwise_test(df, c(assay.type, col.var, row.var), 
                                    group, facet.by, pair.by,
                                    comp.by, features, ...)
    return(res)
})

#' @rdname getPairwiseTest
#' @export
setMethod("addPairwiseTest", signature(x = "SummarizedExperiment"), 
    function( x, name = "pairwiseTest", ...) {
    # get pairwise test
    res <- getPairwiseTest( x, ...)
    
    # Add results to metadata
    x <- .add_values_to_metadata(x, name, res)
    return(x)
})

################################ HELP FUNCTIONS ################################

# Validate user input for pairwise testing
.check_input_for_pairwise_test <- function(
    tse, assay.type, features, row.var, col.var, x, group,
    pair.by = NULL, facet.by = NULL, comp.by = NULL, 
    mark.significance = FALSE, ...){
    
    # Either assay.type. row.var or col.var must be specified
    if( sum(c(is.null(assay.type), is.null(row.var), is.null(col.var))) != 2L ){
        stop("Please specify either 'assay.type', 'row.var', or 'col.var'.",
             call. = FALSE)
    }
    # features cannot be specified if row.var or col.var is specified
    if( is.null(assay.type) && !is.null(features) ){
        stop("'features' can only be specified when 'assay.type is ",
             "specified.", call. = FALSE)
    }
    # As features points to rownames, the TreeSE must have rownames and features
    # must match them
    if( !is.null(features) && is.null(rownames(tse)) ){
        stop("'object' must have rownames.", call. = FALSE)
    }
    if( !(is.null(features) ||
          (is.character(features) && all(features %in% rownames(tse)) )) ){
        stop("'features' must be NULL or single character value specifying ",
             " rownames.", call. = FALSE)
    }
    # If assay was specified, check that it is correct.
    if( !is.null(assay.type) ){
        .check_assay_present(assay.type, tse)
    }
    if( !.is_a_bool(mark.significance) ){
        stop("'mark.significance' must be TRUE or FALSE.", call. = FALSE)
    }
    # Check colData/rowData variables
    temp <- .check_metadata_variable(tse, row.var, row = TRUE)
    temp <- .check_metadata_variable(tse, col.var, col = TRUE)
    temp <- .check_metadata_variable(
        tse, group,
        row = length(c(row.var, assay.type))>0,
        col = length(c(col.var, assay.type))>0
    )
    temp <- .check_metadata_variable(
        tse, comp.by,
        row = length(c(row.var, assay.type))>0,
        col = length(c(col.var, assay.type))>0
    )
    temp <- .check_metadata_variable(
        tse, facet.by,
        row = length(c(row.var, assay.type))>0,
        col = length(c(col.var, assay.type))>0
    )
    temp <- .check_metadata_variable(tse, pair.by, col = TRUE)

    # Check that repeated measures are balanced when 'pair.by' is used
    if (!is.null(pair.by)) {
        coldata <- as.data.frame(colData(tse))
        
        sample_counts <- coldata %>%
            group_by(.data[[pair.by]]) %>%
            summarise(n_samples = n(), .groups = "drop")
        
        if (length(unique(sample_counts$n_samples)) > 1) {
            stop(
                "When 'pair.by' is specified, all subjects must have the ",
                "same number of samples.\nRemove 'pair.by' or filter your ",
                "data to include only subjects with balanced repeats.",
                call. = FALSE
            )
        }
    }
    
    # x must be character or factor in box plots
    if( !is.null(group) && is.numeric(tse[[group]]) ){
        stop("'x' must specify categorical value.", call. = FALSE)
    }
    return(NULL)
}


# This function retrieves the data from TreeSE and outputs a data.frame, ready
# for tests.
.get_data_for_pairwise_test <- function(
    tse, group = NULL, facet.by = NULL, comp.by = NULL, 
    assay.type, row.var, col.var, pair.by = NULL, ...){
    
    # If assay.type is specified, get melted data
    all_vars <- c(group, facet.by, comp.by)
    if( !is.null(assay.type) ){
        # Specify whether to retrieve data from rowData or colData
        row_vars <- vapply(all_vars, function(x){
            x %in% colnames(rowData(tse))
        }, logical(1L))
        col_vars <- all_vars[ !row_vars ]
        row_vars <- all_vars[ row_vars ]
        #
        df <- meltSE(
            tse, assay.type = assay.type,
            col.var = "id",
            add.col = c(col.var, pair.by, col_vars),
            add.row = c(row.var, row_vars)
        )
    }
    # If row.var was specified, get the data from rowData
    if( !is.null(row.var) ){
        df <- rowData(tse)[, c(row.var, all_vars), drop = FALSE]
    }
    # If col.var was specified, get the data from colData
    if( !is.null(col.var) ){
        df <- colData(tse)[, c(col.var, pair.by, all_vars), drop = FALSE]
    }
    # Check that dependent variable is numeric
    if( !is.numeric(df[[c(assay.type, col.var, row.var)]]) ){
        stop("The Y-variable must be numeric.", call. = FALSE)
    }

    return(df)
}

.calculate_pairwise_test <- function(df, y, group, facet.by, pair.by, comp.by, features,
                                     significance.method = "wilcoxon", 
    p.adjust.method = "fdr", paired = !is.null(pair.by),
       include.effect = TRUE, mark.significance = FALSE,
       digits = 3, ...) {
    # Basic validation
    if (!.is_a_bool(paired)) {
        stop("'paired' must be TRUE or FALSE.", call. = FALSE)
    }
    
    supported_methods <- c("wilcoxon", "wilcox.test", "wilcoxon_test","wilcox-test",
                           "ttest", "t.test", "t_test", "t-test",
                           "kruskal", "kruskal.test", "kruskal_test", "kruskal-test")
    if ( !(.is_a_string(significance.method) && significance.method %in% supported_methods) ) {
        stop("'significance.method' must be one of: ",
             paste0("'", supported_methods, "'", collapse = ", "), call. = FALSE)
    }
    
    if ( !.is_a_string(p.adjust.method) ) {
        stop("'p.adjust.method' must be a character string.", call. = FALSE)
    }
    
    if ( !.is_a_bool(mark.significance) ) {
        stop("'mark.significance' must be TRUE or FALSE.", call. = FALSE)
    }
    
    if ( !.is_an_integer(digits) ) {
        stop("'digits' must be a single integer value.", call. = FALSE)
    }
    # Get variables. group will specify the grouping variables, i.e.,
    # we test the significance for these groups separately.
    grouping_vars <- c(facet.by, group)
    # comp.by specify groups for comparison. If they are not
    # specified, we make comparison between group.
    comparison_vars <- c(comp.by) |> unique()
    if( is.null(comparison_vars) ){
        comparison_vars <- group
    }
    grouping_vars <- grouping_vars[ !grouping_vars %in% comparison_vars ]
    
    formula <- as.formula(paste(y, "~", comparison_vars))
    
    if ( significance.method %in% c("wilcoxon", "wilcox.test") ) {
        FUN <- rstatix::pairwise_wilcox_test
        effect_fun <- rstatix::wilcox_effsize
    } else if ( significance.method %in% c("t.test", "ttest") ) {
        FUN <- rstatix::pairwise_t_test
        effect_fun <- rstatix::cohens_d
    } else if ( significance.method %in% c("kruskal", "kruskal.test") ) {
        
        if (paired) {
            stop("Kruskal-Wallis does not support paired comparisons.", call. = FALSE)
        }
        
        global <- rstatix::kruskal_test(formula, data = df, ...)
        pw <- rstatix::dunn_test(
            formula, data = df, p.adjust.method = p.adjust.method, ...
        )
        effect_fun <- rstatix::kruskal_effsize
        pw <- .add_group_means_log2fc(df, group, pw, y)
        return(list(global = global, pairwise = pw))
    }
    
    res <- df |>
        as.data.frame() |>
        # Create grouping to test these groups separately
        dplyr::group_split(across(all_of(grouping_vars))) |>
        purrr::map_df(function(df_group) {
            tryCatch({
                # If we calculate paired analysis, sort data so that correct
                # samples are matched with correct subjects.
                df_group <- df_group %>% arrange(across(all_of(pair.by)))
                # This runs pairwise comparisons automatically if there are more
                # than 2 groups
                res_test <- FUN(
                    formula = as.formula(paste0(y, " ~ ", comparison_vars)),
                    data = df_group,
                    paired = paired,
                    p.adjust.method = "none"
                )
                # Add grouping cols back to pairwise test results
                res_test <- dplyr::bind_cols(
                    res_test,
                    df_group %>% select(all_of(grouping_vars)) %>% distinct()
                )
                if ( include.effect ) {
                    # Calculate effect sizes for this group
                    res_eff <- effect_fun(
                        formula,
                        data = df_group,
                        paired = paired,
                        ...
                    )
                    # Add grouping cols back to effect size results
                    res_eff <- dplyr::bind_cols(
                        res_eff,
                        df_group %>% select(all_of(grouping_vars)) %>% distinct()
                    )
                    # Merge effect sizes into pairwise test results by group1/group2
                    res_merged <- merge(res_test, res_eff, 
                        by = intersect(names(res_test), names(res_eff)))
                    return(res_merged)
                } else {
                    return(res_test)
                }
            }, error = function(e) NULL)
        }) |>
        rstatix::adjust_pvalue(method = p.adjust.method)
    
    res <- .add_group_means_log2fc(df, c(grouping_vars, comparison_vars), res, y)
    
    # Mark or round significance
    if (mark.significance) {
        res <- rstatix::add_significance(res)
    } else {
        res <- rstatix::p_round(res, digits = digits)
    }
    
    # Subset to relevant features if needed
    if (!is.null(features)) {
        if (!is.null(group)) {
            res <- res[res[["rownames"]] %in% features, , drop = FALSE]
            attr_data <- attributes(res)[["args"]][["data"]]
            attr_data <- attr_data[attr_data[["rownames"]] %in% features, , drop = FALSE]
            attributes(res)[["args"]][["data"]] <- attr_data
        } else if ("FeatureID" %in% colnames(res)) {
            res <- res[
                res$group1 %in% features & res$group2 %in% features,
                , drop = FALSE
            ]
            attr_data <- attributes(res)[["args"]][["data"]]
            attr_data <- attr_data[attr_data[["FeatureID"]] %in% features, , drop = FALSE]
            attributes(res)[["args"]][["data"]] <- attr_data
        }
    }
    
    return(res)
}

.add_group_means_log2fc <- function(df, group, res, y) {
    # Compute means for every group (can be 1 or multiple columns)
    means_df <- df |>
        as.data.frame() |>
        dplyr::group_by(across(all_of(group))) |>
        dplyr::summarise(mean = mean(.data[[y]], na.rm = TRUE), .groups = "drop")
    
    comparison_vars <- setdiff(group, colnames(res))
    grouping_vars <- setdiff(group, comparison_vars)
    
    res <- dplyr::left_join(res, means_df, 
        by = c(grouping_vars, setNames(comparison_vars, "group1"))) %>%
        dplyr::rename(mean_group1 = mean)
    
    res <- dplyr::left_join(res, means_df, 
        by = c(grouping_vars, setNames(comparison_vars, "group2"))) %>%
        dplyr::rename(mean_group2 = mean)
    
    # Handle zero values in log2FC calculation
    res$log2FC <- with(res, {
        ifelse(mean_group1 == 0 | mean_group2 == 0, 
               NA_real_,
               log2((mean_group2 + 1e-8) / (mean_group1 + 1e-8)))
    })
    return(res)
}

# This function checks whether variable can be found from colData or rowData.
.check_metadata_variable <- function(
        tse, var, row = FALSE, col = FALSE, multiple = FALSE,
        var.name = .get_name_in_parent(var)){
    if( !.is_a_bool(multiple) ){
        stop("'multiple' must be TRUE or FALSE.", call. = FALSE)
    }
    # If the variable is not NULL
    if( !is.null(var) ){
        # It must be a string and found from colData/rowData
        is_string <- ifelse(multiple, is.character(var), .is_a_string(var))
        check_values <- c()
        check_values <- c(check_values, if(col) colnames(colData(tse)))
        check_values <- c(check_values, if(row) colnames(rowData(tse)))
        var_found <- all( var %in% check_values )
        if( !(is_string && var_found) ){
            stop("'", var.name, "' must be ", ifelse(multiple, "", "a single "),
                 "character value from the following options: '",
                 paste0(check_values, collapse = "', '"), "'", call. = FALSE)
        }
    }
    return(NULL)
}