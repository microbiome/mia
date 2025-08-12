#' @title Get group-wise statistical comparisons
#' 
#' @description
#' Provides modular functions for pairwise comparisons, omnibus tests, and 
#' posthoc analyses with support for both parametric and non-parametric methods.
#'
#' \itemize{
#'   \item \code{getPairwiseDA()}: Two-group comparisons (Wilcoxon, t-test)
#'   \item \code{getOmnibusDA()}: Multi-group omnibus tests (Kruskal-Wallis, Friedman)
#'   \item \code{getPosthocDA()}: Posthoc pairwise comparisons after omnibus testing
#'   \item \code{getDA()}: Unified interface dispatching to appropriate test type
#'   \item \code{add*DA()}: Variants that store results in object metadata
#' }
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
#' either from \code{rowData(x)} or \code{colData(x)}.
#' 
#' @param comp.by \code{NULL} or \code{character scalar}. Specifies a
#' variable from \code{colData(x)} or \code{rowData(x)} which is used to
#' compare observations. (Default: \code{NULL})
#' 
#' @param facet.by \code{NULL} or \code{character scalar}. Specifies a
#' variable from \code{colData(x)} or \code{rowData(x)} which is used to facet
#' or group observations. (Default: \code{NULL})
#' 
#' @param pair.by \code{NULL} or \code{character scalar}.Specifies a
#' variable from \code{colData(x)} which identifies paired observations 
#' (e.g., subject ID). (Default: \code{"NULL"})
#'
#' @param da.method \code{character scalar}. Statistical method to 
#' use: one of \code{"wilcoxon"}, \code{"t.test"}, \code{"kruskal"}, 
#' or \code{"Friedman"}. (Default: \code{"wilcoxon"})
#' 
#' @param da.type Character scalar for \code{getDA()} and \code{addDA()}. 
#' One of "pairwise", "omnibus", or "posthoc". (Default: varies by function)
#'   
#' @param p.adjust.method \code{character scalar}. Method for p-value adjustment. 
#' Must be one of \code{p.adjust.methods}. (Default: \code{"fdr"})
#'
#' @param include.effect \code{logical scalar}. Whether to compute effect sizes. 
#' (Default: \code{TRUE})
#'
#' @param ... Additional arguments passed to the test function (e.g. 
#' \code{alternative = "greater"}).
#'
#' @section Test Methods:
#' \describe{
#'   \item{\strong{Pairwise Tests}}{
#'     \itemize{
#'       \item \strong{Wilcoxon}: Non-parametric test for two groups
#'         \itemize{
#'           \item Unpaired: Mann-Whitney U test
#'           \item Paired: Wilcoxon signed-rank test
#'           \item Effect size: Rank-biserial correlation (r)
#'         }
#'       \item \strong{t-test}: Parametric test for two groups
#'         \itemize{
#'           \item Unpaired: Student's t-test (equal var) or Welch t-test
#'           \item Paired: Paired t-test
#'           \item Effect size: Cohen's d (Hedges-adjusted)
#'         }
#'     }
#'   }
#'   \item{\strong{Omnibus Tests}}{
#'     \itemize{
#'       \item \strong{Kruskal-Wallis}: Non-parametric ANOVA for ≥3 independent groups
#'         \itemize{
#'           \item Effect size: Eta-squared (η²)
#'           \item Posthoc: Dunn's test with rank-based effect sizes
#'         }
#'       \item \strong{Friedman}: Non-parametric repeated-measures ANOVA for ≥3 paired groups
#'         \itemize{
#'           \item Requires \code{pair.by} variable
#'           \item Effect size: Kendall's W (coefficient of concordance)
#'           \item Posthoc: Pairwise Wilcoxon signed-rank tests
#'         }
#'     }
#'   }
#' }
#'
#' @return
#' \describe{
#'   \item{\code{get*DA()}}{A \code{tibble} with test results containing:}
#'   \item{}{
#'     \itemize{
#'       \item \code{group1}, \code{group2}: Comparison groups
#'       \item \code{p}, \code{p.adj}: Raw and adjusted p-values
#'       \item \code{effsize}: Effect size measure
#'       \item \code{mean_group1}, \code{mean_group2}: Group means
#'       \item \code{log2FC}: Log2 fold change (group2 vs group1)
#'       \item Additional columns: \code{n1}, \code{n2}, method-specific statistics
#'     }
#'   }
#'   \item{\code{add*DA()}}{Input object with results added to \code{metadata()}}
#' }
#'
#' For omnibus tests (\code{getOmnibusDA()}), additional attributes:
#' \itemize{
#'   \item \code{attr(result, "global")}: Overall omnibus test results
#'   \item \code{attr(result, "effect")}: Global effect size (if requested)
#' }
#'
#' @seealso 
#' \code{\link[rstatix]{pairwise_wilcox_test}},
#' \code{\link[rstatix]{pairwise_t_test}},
#' \code{\link[rstatix]{kruskal_test}},
#' \code{\link[rstatix]{friedman_test}},
#' \code{\link[rstatix]{cohens_d}},
#' \code{\link[rstatix]{wilcox_effsize}}
#'
#' @name getDA
#' @export
#'
#' @examples
#' library(dplyr)
#' # Load example data
#' data(GlobalPatterns, package = "mia")
#' tse <- GlobalPatterns
#' 
#' # Transform to relative abundances
#' tse <- transformAssay(tse, method = "relabundance")
#' tse <- tse[690:700, ]
#' 
#' result_pairwise <- getPairwiseDA(
#'     tse,
#'     assay.type = "relabundance",
#'     group = "SampleType",
#'     da.method = "wilcoxon"
#' )
#' head(result_pairwise)
#' 
#' # Multi-group omnibus test
#' result_omnibus <- getOmnibusDA(
#'     tse,
#'     assay.type = "relabundance", 
#'     group = "SampleType",
#'     da.method = "kruskal"
#' )
#' print(result_omnibus)
#' 
#' # Unified interface
#' result_unified <- getDA(
#'     tse_two,
#'     da.type = "pairwise",
#'     assay.type = "relabundance",
#'     group = "SampleType"
#' )
#' 
#' # Add results to metadata
#' tse <- addPairwiseDA(
#'     tse_two,
#'     assay.type = "relabundance",
#'     group = "SampleType", 
#'     name = "wilcoxon_results"
#' )
#' 
#' @seealso
#' \code{\link[rstatix]{pairwise_wilcox_test}},
#' \code{\link[rstatix]{pairwise_t_test}},
#' \code{\link[rstatix]{kruskal_test}}
#' 
NULL

#' @rdname getDA
#' @export
setMethod("getDA", signature(x = "SummarizedExperiment"),
    function(x, da.type = c("pairwise", "omnibus", "posthoc"), ...) {
        da.type <- match.arg(da.type)
        FUN <- switch(
            da.type,
            pairwise = getPairwiseDA,
            omnibus  = getOmnibusDA,
            posthoc  = getPosthocDA)
        FUN(x, ...)
    }
)

#' @rdname getDA
#' @export
setMethod("addDA", signature(x = "SummarizedExperiment"),
    function(x, da.type = c("pairwise", "omnibus", "posthoc"),
           name = "DAtest", ...) {
        da.type <- match.arg(da.type)
        res <- getDA(x, da.type, ...)
        x <- .add_values_to_metadata(x, name, res)
        return(x)
    }
)

#' @rdname getDA
#' @export
setMethod("getPairwiseDA", signature(x = "SummarizedExperiment"), 
    function(x, assay.type = NULL, features = NULL, row.var = NULL, 
            col.var = NULL, group, facet.by = NULL, comp.by = NULL,
            pair.by = NULL, da.method = "wilcoxon", ...) {
        .run_DA(x, assay.type, features, row.var, col.var, group, 
               facet.by, comp.by, pair.by, da.method, ...)
    }
)

#' @rdname getDA
#' @export
setMethod("addPairwiseDA", signature(x = "SummarizedExperiment"), 
    function(x, name = "pairwiseDA", ...) {
        res <- getPairwiseDA(x, ...)
        x <- .add_values_to_metadata(x, name, res)
        return(x)
    }
)

#' @rdname getDA
#' @export
setMethod("getOmnibusDA", signature(x = "SummarizedExperiment"), 
    function(x, assay.type = NULL, features = NULL, row.var = NULL, 
            col.var = NULL, group, facet.by = NULL, comp.by = NULL,
            pair.by = NULL, da.method = "kruskal", ...) {
        .run_DA(x, assay.type, features, row.var, col.var, group, 
               facet.by, comp.by, pair.by, da.method, ...)
    }
)

#' @rdname getDA
#' @export
setMethod("addOmnibusDA", signature(x = "SummarizedExperiment"), 
    function(x, name = "omnibusDA", ...) {
        res <- getOmnibusDA(x, ...)
        x <- .add_values_to_metadata(x, name, res)
        return(x)
    }
)

#' @rdname getDA
#' @export
setMethod("getPosthocDA", signature(x = "SummarizedExperiment"), 
    function(x, assay.type = NULL, features = NULL, row.var = NULL, 
            col.var = NULL, group, facet.by = NULL, comp.by = NULL,
            pair.by = NULL, da.method = "wilcoxon", ...) {
        .run_DA(x, assay.type, features, row.var, col.var, group, 
               facet.by, comp.by, pair.by, da.method, ...)
    }
)

#' @rdname getDA
#' @export
setMethod("addPosthocDA", signature(x = "SummarizedExperiment"), 
    function(x, name = "posthocDA", ...) {
        res <- getPosthocDA(x, ...)
        x <- .add_values_to_metadata(x, name, res)
        return(x)
    }
)

################################ HELP FUNCTIONS ################################

# Validate user input for pairwise testing
.check_input_for_DA <- function(
    tse, assay.type, features, row.var, col.var, x, group,
    pair.by = NULL, facet.by = NULL, comp.by = NULL, da.method,
    paired = !is.null(pair.by), ...
) {
    # Either assay.type. row.var or col.var must be specified
    if( sum(c(is.null(assay.type), is.null(row.var), is.null(col.var))) != 2L ){
        stop("Please specify either 'assay.type', 'row.var', or 'col.var'.",
             call. = FALSE)
    }
    # features cannot be specified if row.var or col.var is specified
    if( is.null(assay.type) && !is.null(features) ){
        stop("'features' can only be specified when 'assay.type' is ",
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
    # validate method
    da.method <- .check_method(da.method)
    # check method-specific requirements
    if ( da.method == "friedman" && !paired ) {
        stop("Friedman test requires a 'pair.by' variable.", call. = FALSE)
    }
    if ( da.method == "kruskal" && paired ) {
        stop("Kruskal-Wallis test does not support paired designs.", call. = FALSE)
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

# This function checks validity of method from user
.check_method <- function(method) {
    supported <- c("wilcoxon", "ttest", "kruskal", "dunns", "friedman")
    method <- tolower(method)
    if ( !(method %in% supported) ) {
        stop("Unsupported method: '", method, "'. Must be one of: ",
             paste(supported, collapse = ", "), call. = FALSE)
    }
    method
}

# This function retrieves the data from TreeSE and outputs a data.frame, ready
# for tests.
.get_data_for_DA <- function(
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

# internal wrapper to run DA
.run_DA <- function(x, assay.type, features, row.var, col.var, group,
                    facet.by, comp.by, pair.by, da.method, ...) {
    
    if (!"rownames" %in% colnames(rowData(x))) {
        rowData(x)[["rownames"]] <- rownames(x)
    }
    
    args <- list(tse = x, assay.type = assay.type, features = features,
                 row.var = row.var, col.var = col.var, group = group,
                 facet.by = facet.by, comp.by = comp.by, pair.by = pair.by, 
                 da.method = da.method, ...)
    
    do.call(.check_input_for_DA, args)
    df <- do.call(.get_data_for_DA, args)
    
    .calc_DA(df, c(assay.type, col.var, row.var), group, facet.by, 
             pair.by, comp.by, features, da.method, ...)
}

#' @importFrom rstatix pairwise_wilcox_test wilcox_effsize pairwise_t_test cohens_d
#' @importFrom rstatix friedman_test friedman_effsize kruskal_test kruskal_effsize
#' @importFrom rstatix dunn_test adjust_pvalue add_significance p_round
#' @importFrom dplyr across all_of distinct group_by summarise left_join 
#' @importFrom dplyr rename arrange select
# Actual function calling various statistical tests
.calc_DA <- function(
    df, y, group, facet.by, pair.by, comp.by, features, 
    da.method, p.adjust.method = "fdr", 
    paired = !is.null(pair.by), include.effect = TRUE, 
    digits = 3, ...
) {
    # Basic validation
    if (!.is_a_bool(paired)) {
        stop("'paired' must be TRUE or FALSE.", call. = FALSE)
    }
    
    if ( !.is_a_string(p.adjust.method) ) {
        stop("'p.adjust.method' must be a character string.", call. = FALSE)
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
    # grouping_vars: used to split data into groups for repeated tests
    # comparison_vars: the variable(s) we are testing across
    if( is.null(comparison_vars) ){
        comparison_vars <- group
    }
    grouping_vars <- grouping_vars[ !grouping_vars %in% comparison_vars ]
    
    formula <- as.formula(paste(y, "~", paste(comparison_vars, collapse = "+")))
    
    res <-  global <- effect <- NULL
    method_funs <- list(
        wilcoxon = .calc_pairwise,
        ttest    = .calc_pairwise,
        kruskal  = .calc_omnibus,
        dunns    = .calc_posthoc,
        friedman = .calc_omnibus
    )
    
    FUN <- method_funs[[da.method]]
    
    res <- FUN(
        df, formula, grouping_vars, comparison_vars, pair.by, p.adjust.method,
        paired, include.effect, y, 
        da.method, ...
    )
    # Safety check: handle empty results
    if (is.null(res) || nrow(res) == 0) {
        warning("Statistical test returned empty results. This may be due to:",
                "\n- Insufficient sample size",
                "\n- All values being identical", 
                "\n- Missing data",
                "\n- Inappropriate test for the data structure",
                call. = FALSE)
        return(res)  # Return empty result early
    }
    
    # Subset to relevant features if needed
    if (!is.null(features)) {
        res <- .subset_features(res, features)
    }
    
    res <- .clean_DA_results(res, da.method)
    return(res)
}

# Calculates pairwise statistics
.calc_pairwise <- function(
    df, formula, grouping_vars, comparison_vars, pair.by, 
    p.adjust.method, paired, include.effect, y, 
    da.method, ...
) {
    # 1. Determine appropriate statistical functions
    if (da.method %in% c("wilcoxon")) {
        test_fun <- pairwise_wilcox_test
        effect_fun <- wilcox_effsize
    } else {
        test_fun <- pairwise_t_test
        effect_fun <- cohens_d
    }
    
    pw <- df |>
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
                res_test <- test_fun(
                    formula = formula,
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
                        df_group %>% 
                            select(all_of(grouping_vars)) %>% 
                            distinct()
                    )
                    # Merge effect sizes into pairwise test results by 
                    # group1/group2
                    res_merged <- merge(res_test, res_eff, 
                                        by = intersect(names(res_test), 
                                                       names(res_eff)))
                    return(res_merged)
                } else {
                    return(res_test)
                }
            }, error = function(e) NULL)
        }) |>
        rstatix::adjust_pvalue(method = p.adjust.method)
    pw <- .add_group_means_log2fc(df, c(grouping_vars, comparison_vars), pw, y)
    return(pw)
}

# Calculates omnibus statistics
.calc_omnibus <- function(
    df, formula, grouping_vars, comparison_vars, pair.by, p.adjust.method,
    paired, include.effect, y, 
    da.method, ...
) {
    effect <- NULL
    
    # ----------------------------- Friedman Test ------------------------------
    if (da.method %in% c("friedman")) {
        friedman_formula <- as.formula(
            paste(y, "~", paste(comparison_vars, collapse = "+"), "|", pair.by)
        )
        global <- friedman_test(
            formula = friedman_formula, 
            data = df, 
            ...
        )
        if (include.effect) {
            effect <- friedman_effsize(friedman_formula, data = df, ...)
        }
        
    # ------------------------------ Kruskal-Wallis Test ----------------------
    } else {
        global <- kruskal_test(formula, data = df, ...)
        
        if (include.effect) {
            effect <- kruskal_effsize(formula, data = df, ...)
        }
    }
    
    if (!is.null(effect)) {
        attr(global, "effect") <- effect
        global <- dplyr::left_join(global, effect, by = c("n", ".y."))
    }
    class(global) <- c("stat_test_result", class(global))
    
    return(global)
}

# Calculates dunn's posthoc
.calc_dunns <- function(df, formula, grouping_vars, comparison_vars,
    p.adjust.method, include.effect,
    y, ...
) {
    pw <- dunn_test(formula = formula, data = df, 
                             p.adjust.method = p.adjust.method)
    
    if (include.effect) {
        eff_pw <- wilcox_effsize(formula = formula, data = df)
        pw <- merge(pw, eff_pw, by = intersect(names(pw), names(eff_pw)))
    }
    
    pw <- .add_group_means_log2fc(df, c(grouping_vars, comparison_vars), pw, y)
    
    class(pw) <- c("stat_test_result", class(pw))
    return(pw)
}

# Calculates general posthoc
.calc_posthoc <- function(
    df, formula, grouping_vars, comparison_vars, pair.by, p.adjust.method,
    paired, include.effect, y, 
    da.method, ...
) {
    # ----------------------------- Friedman Test with Posthoc -----------------
    if (da.method %in% c("friedman", "wilcoxon")) {
        friedman_formula <- as.formula(
            paste(y, "~", paste(comparison_vars, collapse = "+"), "|", pair.by)
        )
        global <- .calc_omnibus(
            df, friedman_formula, grouping_vars, comparison_vars, pair.by, 
            p.adjust.method, paired, include.effect, y, 
            da.method = "friedman", ...
        )

        pw <- .calc_pairwise(
            df, formula, grouping_vars, comparison_vars, pair.by, 
            p.adjust.method, paired, include.effect, 
            y, da.method, ...
        )

    # ------------------------ Kruskal-Wallis with Dunn's Test ----------------
    } else {
        global <- .calc_omnibus(
            df, formula, grouping_vars, comparison_vars, pair.by = NULL, 
            p.adjust.method, paired, include.effect, y, 
            da.method = "kruskal", ...
        )

        pw <- .calc_dunns(
            df, formula, grouping_vars, comparison_vars, p.adjust.method, 
            include.effect, y, 
            da.method, ...
        )
    }

    attr(pw, "global_test") <- global
    return(pw)
}

# Calculates group means and fold change
.add_group_means_log2fc <- function(df, group, res, y) {
    # Safety check for empty results
    if (is.null(res) || nrow(res) == 0) {
        return(res)
    }
    
    # Safety check for required columns
    if (!all(c("group1", "group2") %in% colnames(res))) {
        warning("Results missing required 'group1' and 'group2' columns. ",
                "Skipping mean and log2FC calculations.", call. = FALSE)
        return(res)
    }
    
    # Compute means for every group (can be 1 or multiple columns)
    means_df <- df |>
        as.data.frame() |>
        dplyr::group_by(across(all_of(group))) |>
        dplyr::summarise(mean = mean(.data[[y]], na.rm = TRUE), .groups = "drop")
    
    # Safety check for means calculation
    if (nrow(means_df) == 0) {
        warning("Could not calculate group means. Check data structure.", 
                call. = FALSE)
        return(res)
    }
    
    comparison_vars <- setdiff(group, colnames(res))
    grouping_vars <- setdiff(group, comparison_vars)
    
    # Safely join means
    tryCatch({
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
    }, error = function(e) {
        warning("Error calculating group means and log2FC: ", e$message, 
                ". Results returned without these columns.", call. = FALSE)
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

# Useful to subset features entered by user
.subset_features <- function(res, features) {
    if (!is.null(features) && nrow(res) > 0) {
        if ("rownames" %in% colnames(res)) {
            res <- res[res[["rownames"]] %in% features, ,
                       drop = FALSE]
            # Also filter attr data if present
            attr_data <- attr(res, "args")$data
            if (!is.null(attr_data)) {
                attr_data <- attr_data[
                    attr_data[["rownames"]] %in% features, ,
                    drop = FALSE
                ]
                attr(res, "args")$data <- attr_data
            }
        } else if ("FeatureID" %in% colnames(res)) {
            res <- res[res[["FeatureID"]] %in% features, ,
                       drop = FALSE]
            attr_data <- attr(res, "args")$data
            if (!is.null(attr_data)) {
                attr_data <- attr_data[
                    attr_data[["FeatureID"]] %in% features, ,
                    drop = FALSE
                ]
                attr(res, "args")$data <- attr_data
            }
        }
    }
    return(res)
}

# Clean DA results before returning
.clean_DA_results <- function(res, da.method) {
    # Safety check for empty or null results
    if (is.null(res) || nrow(res) == 0) {
        return(res)
    }
    
    # Save attributes before modifying
    attr_list <- attributes(res)
    
    if (da.method %in% c("wilcoxon", "ttest", "dunns")) {
        res <- res[, !names(res) %in% c("statistic", "magnitude"), drop = FALSE]
        desired_order <- c(".y.", "group1", "group2", "p", "p.adj", 
                           "effsize", "log2FC", "mean_group1", "mean_group2", 
                           "n1", "n2", "df", "method")
        
    } else if (da.method %in% c("kruskal", "friedman")) {
        res <- res[, !names(res) %in% c("magnitude"), drop = FALSE]
        desired_order <- c(".y.", "n", "statistic", 
                           "df", "p", "method", "effsize")
        
    } else {
        return(res)  # no cleaning
    }
    
    existing_cols <- intersect(desired_order, names(res))
    remaining_cols <- setdiff(names(res), existing_cols)
    final_order <- c(existing_cols, remaining_cols)
    res <- res[, final_order, drop = FALSE]
    
    # Restore original attributes (except dimensions)
    attributes(res) <- modifyList(attr_list, attributes(res))
    
    return(res)
}

