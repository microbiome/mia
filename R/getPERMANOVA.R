#' @name
#' getPERMANOVA
#'
#' @title
#' Calculate PERMANOVA (Permutational Multivariate Analysis of Variance)
#'
#' @description
#' These functions perform PERMANOVA to assess the significance of group
#' differences based on a specified dissimilarity matrix. The results can be
#' returned directly or added to metadata in an object of class
#' \code{TreeSummarizedExperiment}.
#'
#' @details
#' PERMANOVA is a non-parametric method used to test whether the centroids of
#' different groups (as defined by the formula or covariates) are significantly
#' different in terms of multivariate space.
#'
#' PERMANOVA relies on the assumption of group homogeneity, meaning the groups
#' should be distinct and have similar variances within each group. This
#' assumption is essential as PERMANOVA is sensitive to differences in
#' within-group dispersion, which can otherwise confound results.
#' This is why the functions return homogeneity test results by default.
#'
#' The functions utilize \code{\link[vegan:adonis2]{vegan::adonis2}} to compute
#' PERMANOVA. For homogeneity testing,
#' \code{\link[vegan:betadisper]{vegan::betadisper}} along
#' with \code{\link[vegan:permutest]{vegan::permutest}} are utilized by default,
#' which allow testing for equal dispersion across groups and validate the
#' homogeneity assumption.
#'
#' PERMANOVA and distance-based redundancy analysis (dbRDA) are closely related
#' methods for analyzing multivariate data. PERMANOVA is non-parametric, making
#' fewer assumptions about the data. In contrast, dbRDA assumes a linear
#' relationship when constructing the ordination space, although it also
#' employs a PERMANOVA-like approach to test the significance of predictors
#' within this space. dbRDA offers a broader scope overall, as it provides
#' visualization of the constrained ordination, which can reveal patterns and
#' relationships. However, when the underlying data structure is non-linear,
#' the results from the two methods can differ significantly due to dbRDA's
#' reliance on linear assumptions.
#'
#' With \code{pairwise = TRUE}, PERMANOVA is run between each pair of groups of
#' the first variable in \code{formula} (or the first \code{col.var}), with the
#' other variables as covariates. Pairwise tests are post-hoc comparisons that
#' follow a significant global test. Their p-values are adjusted across the
#' pairs with \code{p.adj.method}.
#'
#' @return
#' \code{getPERMANOVA} returns the PERMANOVA results or a list containing the
#' PERMANOVA results and homogeneity test results if
#' \code{test.homogeneity = TRUE}. With \code{pairwise = TRUE}, each result is
#' a \code{data.frame} with one row per pair of groups: the groups
#' (\code{group1}, \code{group2}), the number of samples (\code{n}), the
#' statistics of the grouping variable and the adjusted p-value
#' (\code{p.adj}). \code{addPERMANOVA} adds these results to metadata of
#' \code{x}.
#'
#' @inheritParams runCCA
#' @inheritParams addAlpha

#' @param name \code{Character scalar}. A name for the results that will be
#' stored to metadata. (Default: \code{"permanova"})
#'
#' @param method \code{Character scalar}. A dissimilarity metric used in
#' PERMANOVA and group dispersion calculation. (Default: \code{"bray"})
#'
#' @param test.homogeneity \code{Logical scalar}. Should the homogeneity of
#' group dispersions be evaluated? (Default: \code{TRUE})
#'
#' @param pairwise \code{Logical scalar}. Should PERMANOVA be run between each
#' pair of groups of the first variable? (Default: \code{FALSE})
#'
#' @param p.adj.method \code{Character scalar}. Specifies adjustment method of
#' p-values when \code{pairwise = TRUE}. Passed to \code{p.adjust} function.
#' (Default: \code{"fdr"})
#'
#' @param ... additional arguments passed to \code{getDissimilarity} and
#' \code{vegan::adonis2}.
#' \itemize{
#'   \item \code{by}: \code{Character scalar}. Specifies how significance is
#'   calculated. (Default: \code{"margin"})
#'
#'   \item \code{na.action}: \code{function}. Action to take when missing
#'   values for any of the variables in \code{formula} are encountered.
#'   (Default: \code{na.fail})
#'
#'   \item \code{full} \code{Logical scalar}. should all the results from the
#'   homogeneity calculations be returned. When \code{FALSE}, only
#'   summary tables are returned. (Default: \code{FALSE})
#'
#'   \item \code{homogeneity.test}: \code{Character scalar}. Specifies
#'   the significance test used to analyse
#'   \code{\link[vegan:betadisper]{vegan::betadisper}} results.
#'   Options include 'permanova'
#'   (\code{\link[vegan:permutest]{vegan::permutest}}), 'anova'
#'   (\code{\link[stats:anova]{stats::anova}}) and 'tukeyhsd'
#'   (\code{\link[stats:TukeyHSD]{stats::TukeyHSD}}).
#'   (Default: \code{"permanova"})
#' }
#'
#' @examples
#'
#' data(GlobalPatterns)
#' tse <- GlobalPatterns
#'
#' # Apply relative transformation
#' tse <- transformAssay(tse, method = "relabundance")
#' # Perform PERMANOVA
#' tse <- addPERMANOVA(
#'     tse,
#'     assay.type = "relabundance",
#'     method = "bray",
#'     formula = x ~ SampleType,
#'     permutations = 99
#'     )
#' # The results are stored to metadata
#' metadata(tse)[["permanova"]]
#'
#' # Calculate dbRDA
#' rda_res <- getRDA(
#'     tse, assay.type = "relabundance", method = "bray",
#'     formula = x ~ SampleType, permutations = 99)
#' # Significance results are similar to PERMANOVA
#' attr(rda_res, "significance")
#'
#' # Run PERMANOVA using a pre-calculated dissimilarity matrix from metadata.
#' # By doing like this, you can use all dissimilarity methods from
#' # getDissimilarity() in PERMANOVA.
#' tse <- addDissimilarity(tse, method = "unifrac", name = "unifrac_dist")
#' res_precalc <- getPERMANOVA(
#'     tse,
#'     dis.name = "unifrac_dist",
#'     formula = x ~ SampleType,
#'     permutations = 99
#' )
#'
#' # Compare each pair of geographical locations, with diet as a covariate.
#' # The p-values are adjusted across the three pairs.
#' data(peerj13075)
#' tse <- transformAssay(peerj13075, method = "relabundance")
#' res_pairwise <- getPERMANOVA(
#'     tse,
#'     assay.type = "relabundance",
#'     formula = x ~ Geographical_location + Diet,
#'     pairwise = TRUE,
#'     permutations = 99
#' )
#' res_pairwise[["permanova"]]
#'
#' @seealso
#' For more details on the actual implementation see
#' \code{\link[vegan:adonis2]{vegan::adonis2}},
#' \code{\link[vegan:betadisper]{vegan::betadisper}}, and
#' \code{\link[vegan:permutest]{vegan::permutest}}. See also
#' \code{\link[=runCCA]{addCCA}} and \code{\link[=runCCA]{addRDA}}
#'
#' @references
#' \insertRef{Anderson2001}{mia}
#'
NULL

#' @export
#' @rdname getPERMANOVA
setMethod("getPERMANOVA", "SingleCellExperiment", function(x,  ...){
    # Get altexp if specified
    x <- .check_and_get_altExp(x, ...)
    res <- callNextMethod(x, ...)
    return(res)
    }
)

#' @export
#' @rdname getPERMANOVA
setMethod("getPERMANOVA", "SummarizedExperiment",
    function(
        x, assay.type = NULL, dis.name = NULL, formula = NULL, col.var = NULL,
         ...){
        ############################# Input check ##############################
        # Formula must be either correctly specified or not specified
        if( !(is.null(formula) || is(formula, "formula")) ){
            stop("'formula' must be formula or NULL.", call. = FALSE)
        }
        # User can also specify covariates with col.var
        if( !(is.null(col.var) || (is.character(col.var) &&
                all(col.var %in% colnames(colData(x))))) ){
            stop("'col.var' must specify column from colData(x) or be NULL.",
                call. = FALSE)
        }
        if( !is.null(formula) && !is.null(col.var) ){
            stop("Specify either 'formula' or 'col.var'.", call. = FALSE)
        }
        # User can specify either abundance matrix or dissimilarity matrix from
        # metadata.
        if( sum(c(is.null(assay.type), is.null(dis.name))) != 1L ){
            stop("Either 'assay.type' or 'dis.name' must be specified.",
                 call. = FALSE)
        }
        if( !is.null(assay.type) ){
            # Get assay
            .check_assay_present(assay.type, x)
            mat <- assay(x, assay.type)
            # UniFrac needs the tree, which is taken from the TreeSE
            if( is(x, "TreeSummarizedExperiment") &&
                    identical(list(...)[["method"]], "unifrac") ){
                mat <- getDissimilarity(x, assay.type = assay.type, ...)
            }
        } else{
            # Get dissimilarity matrix
            .check_metadata_present(dis.name, x)
            mat <- metadata(x)[[dis.name]]
        }
        ########################### Input check end ############################
        # Get formula and related sample metadata as DF
        temp <- .get_formula_and_covariates(x, formula, col.var)
        formula <- temp[["formula"]]
        covariates <- temp[["variables"]]
        # Calculate PERMANOVA with matrix/dist method
        res <- getPERMANOVA(
            mat, formula = formula, data = covariates, ...)
        return(res)
    }
)

#' @export
#' @rdname getPERMANOVA
#' @importFrom stats p.adjust.methods
setMethod("getPERMANOVA", "ANY", function(
        x, formula, data, method = "bray", test.homogeneity = TRUE,
        pairwise = FALSE, p.adj.method = "fdr", ...){
    is_dist <- inherits(x, "dist")
    if( !(is.matrix(x) || is_dist) ){
        stop("'x' must be a matrix or a 'dist' object.", call. = FALSE)
    }
    if( !is(formula, "formula") ){
        stop("'formula' must be formula or NULL.", call. = FALSE)
    }
    if( !(is.data.frame(data) || is.matrix(data) || is(data, "DFrame")) ){
        stop("'data' must be data.frame or coercible to one.", call. = FALSE)
    }
    n_samples <- if( is_dist ) attr(x, "Size") else ncol(x)
    if( n_samples != nrow(data) ){
        stop("Number of samples in 'x' should match with number of rows in ",
            "'data'.", call. = FALSE)
    }
    if( !.is_a_string(method) ){
        stop("'method' must be a single character value.", call. = FALSE)
    }
    if( !.is_a_bool(test.homogeneity) ){
        stop("'test.homogeneity' must be TRUE or FALSE.", call. = FALSE)
    }
    if( !.is_a_bool(pairwise) ){
        stop("'pairwise' must be TRUE or FALSE.", call. = FALSE)
    }
    p.adj.method <- match.arg(p.adj.method, p.adjust.methods)
    # Run the analysis below for each pair of groups
    if( pairwise ){
        res <- .calculate_pairwise_permanova(
            x, formula, data, method, test.homogeneity, p.adj.method, ...)
        return(res)
    }
    # Calculate dissimilarities
    if( !is_dist ){
        x <- getDissimilarity(x |> t(), method = method, ...)
    }
    # Calculate PERMANOVA
    res <- .calculate_permanova(x, formula = formula, data = data, ...)
    # Test homogeneity
    if( test.homogeneity ){
        homogeneity <- .calculate_homogeneity(x, data, ...)
        res <- list(permanova = res, homogeneity = homogeneity)
    }
    return(res)
})

#' @export
#' @rdname getPERMANOVA
setMethod("addPERMANOVA", "SummarizedExperiment",
    function(x, name = "permanova", ...){
    if( !.is_a_string(name) ){
        stop("'name' must be a single character value.", call. = FALSE)
    }
    # Calculate permanova
    res <- getPERMANOVA(x, ...)
    # Add results to metadata
    x <- .add_values_to_metadata(x, name, res)
    return(x)
    }
)

################################ HELP FUNCTIONS ################################

# This function is internal function to perform PERMANOVA from abundance
# matrix or distance matrix, formula and sample metadata table.
#' @importFrom vegan adonis2
.calculate_permanova <- function(
        x, formula, data, by = "margin", na.action = na.fail, ...){
    #
    if( !.is_a_string(by) ){
        stop("'by' must be a single character value.", call. = FALSE)
    }
    #
    # Get abundance data into correct orientation. Samples must be in rows and
    # features in columns.
    if( is.matrix(x) ){
        x <- x |> t()
    }
    # Covariate table must be data.frame. Unused factor levels (e.g. after
    # subsetting samples) are dropped: vegan::adonis2() can crash R with them.
    data <- droplevels(data.frame(data, check.names = FALSE))
    # A grouping variable with one group cannot explain any variation
    vars <- intersect(all.vars(formula), colnames(data))
    one_group <- vapply(data[vars], function(v)
        (is.factor(v) || is.character(v)) && length(unique(na.omit(v))) < 2L,
        logical(1L))
    if( any(one_group) ){
        stop("Variables with only one group cannot be tested: ",
            paste0("'", vars[one_group], "'", collapse = ", "), ".",
            call. = FALSE)
    }
    # Instead of letting na.action pass through, give informative error
    # about missing values.
    if( any(is.na(data)) && isTRUE(all.equal(na.action, na.fail)) ){
        stop("Variables contain missing values. Set na.action to na.exclude ",
            "to remove samples with missing values.", call. = FALSE)
    }
    # This next step ensures that the formula points correctly to abundance
    # table or distance matrix
    formula[[2]] <- as.name("x")
    environment(formula) <- list2env(
        list(x = x), parent = environment(formula))
    # Calculate PERMANOVA
    res <- adonis2(
        formula = formula, data = data, by = by, na.action = na.action, ...)
    return(res)
}

# This function runs getPERMANOVA() for each pair of groups in the first
# variable of the formula and adjusts the p-values across the pairs.
#' @importFrom stats as.dist
#' @importFrom utils combn
.calculate_pairwise_permanova <- function(
        x, formula, data, method, test.homogeneity, p.adj.method, ...){
    data <- data.frame(data, check.names = FALSE)
    var <- all.vars(formula[[3]])[1L]
    group <- data[[var]]
    if( !(is.factor(group) || is.character(group) || is.logical(group)) ||
            nlevels(factor(group)) < 2L ){
        stop("With 'pairwise = TRUE', the first variable on the right-hand ",
            "side of 'formula' must be categorical with at least two groups.",
            call. = FALSE)
    }
    group <- factor(group)
    # The grouping variable comes first so that its homogeneity result is the
    # first row
    data <- data[, c(var, setdiff(colnames(data), var)), drop = FALSE]
    pairs <- combn(levels(group), 2L, simplify = FALSE)
    res <- lapply(pairs, function(pair){
        # Samples with a missing group are kept so that 'na.action' handles
        # them as in the global test
        keep <- is.na(group) | group %in% pair
        x_sub <- if( is.matrix(x) ) x[, keep, drop = FALSE] else
            as.dist(as.matrix(x)[keep, keep])
        getPERMANOVA(
            x_sub, formula = formula, data = data[keep, , drop = FALSE],
            method = method, test.homogeneity = test.homogeneity, ...)
    })
    info <- data.frame(
        group1 = vapply(pairs, "[", character(1L), 1L),
        group2 = vapply(pairs, "[", character(1L), 2L),
        n = vapply(pairs, function(pair) sum(group %in% pair), integer(1L)))
    if( !test.homogeneity ){
        res <- .combine_pairwise_tables(res, info, p.adj.method, var)
        return(res)
    }
    homogeneity <- lapply(res, function(r) r[["homogeneity"]])
    if( isTRUE(list(...)[["full"]]) ){
        homogeneity <- lapply(homogeneity, function(h) h[["summary"]])
    }
    res <- list(
        permanova = .combine_pairwise_tables(
            lapply(res, function(r) r[["permanova"]]), info, p.adj.method, var),
        homogeneity = .combine_pairwise_tables(
            homogeneity, info, p.adj.method))
    return(res)
}

# This function takes the row of the grouping variable from each pair's result
# table (the first row for homogeneity results) and adds the pair information
# and adjusted p-values.
#' @importFrom stats p.adjust
.combine_pairwise_tables <- function(tabs, info, p.adj.method, row = 1L){
    if( is.character(row) &&
            !all(vapply(tabs, function(t) row %in% rownames(t), logical(1L))) ){
        stop("With 'pairwise = TRUE', the PERMANOVA table must have a row ",
            "for '", row, "'. Use by = \"terms\" or by = \"margin\" ",
            "without interactions of the grouping variable.", call. = FALSE)
    }
    tab <- do.call(rbind, lapply(tabs, function(t) as.data.frame(t)[row, ]))
    rownames(tab) <- NULL
    # The p-value column is "p adj" for TukeyHSD and "Pr(>F)" otherwise
    p <- tab[[intersect(c("Pr(>F)", "p adj"), colnames(tab))]]
    res <- cbind(info, tab, p.adj = p.adjust(p, method = p.adj.method))
    return(res)
}
