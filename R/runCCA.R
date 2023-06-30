#' Canonical Correspondence Analysis
#'
#' These functions perform Canonical Correspondence Analysis on data stored
#' in a \code{SummarizedExperiment}.
#'
#' @param x For \code{calculate*} a 
#'   \code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}} 
#'   or a numeric matrix with columns as samples 
#'
#'   For \code{run*} a
#'   \code{\link[SingleCellExperiment:SingleCellExperiment]{SingleCellExperiment}}
#'   or a derived object.
#'
#' @param formula If \code{x} is a
#'   \code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#'   a formula can be supplied. Based on the right-hand side of the given formula
#'   \code{colData} is subset to \code{variables}.
#'   
#'   \code{variables} and \code{formula} can be missing, which turns the CCA analysis 
#'   into a CA analysis and dbRDA into PCoA/MDS.
#'
#' @param variables When \code{x} is a \code{SummarizedExperiment},
#'   \code{variables} can be used to specify variables from \code{colData}. 
#'   
#'   When \code{x} is a matrix, \code{variables} is a \code{data.frame} or 
#'   an object coercible to one containing the variables to use. 
#'    
#'   All variables are used. Please subset, if you want to consider only some of them. 
#'   \code{variables} and \code{formula} can be missing, which turns the CCA analysis 
#'   into a CA analysis and dbRDA into PCoA/MDS.
#' 
#' @param test_signif a logical scalar, should the PERMANOVA and analysis of
#'   multivariate homogeneity of group dispersions be performed.
#'   (By default: \code{test_signif = TRUE})
#'   
#' @param scale a logical scalar, should the expression values be standardized?
#'   \code{scale} is disabled when using \code{*RDA} functions. Please scale before
#'   performing RDA (Check examples.)
#'
#' @param assay.type a single \code{character} value for specifying which
#'   assay to use for calculation.
#'
#' @param exprs_values a single \code{character} value for specifying which
#'   assay to use for calculation.
#'   (Please use \code{assay.type} instead.)
#'   
#' @param assay_name a single \code{character} value for specifying which
#'   assay to use for calculation.
#'   (Please use \code{assay.type} instead. At some point \code{assay_name}
#'   will be disabled.)
#'   
#' @param altexp String or integer scalar specifying an alternative experiment
#'   containing the input data.
#'
#' @param name String specifying the name to be used to store the result in the
#'   reducedDims of the output.
#'
#' @param ... additional arguments passed to vegan::cca or vegan::dbrda
#' 
#' @details
#'   *CCA functions utilize \code{vegan:cca} and *RDA functions \code{vegan:dbRDA}.
#'   By default dbRDA is done with euclidean distances which equals to RDA.
#'   
#'   Significance tests are done with \code{vegan:anova.cca} (PERMANOVA), and
#'   \code{vegan:betadisper} and \code{stats:anova}
#'   (multivariate homogeneity of groups dispersions (variances)).
#'
#' @return
#' For \code{calculateCCA} a matrix with samples as rows and CCA dimensions as
#' columns
#'
#' For \code{runCCA} a modified \code{x} with the results stored in
#' \code{reducedDim} as the given \code{name}
#'
#' @name runCCA
#' @seealso
#' For more details on the actual implementation see \code{\link[vegan:cca]{cca}}
#' and \code{\link[vegan:dbrda]{dbrda}}
#'
#' @examples
#' library(scater)
#' data(GlobalPatterns)
#' GlobalPatterns <- runCCA(GlobalPatterns, data ~ SampleType)
#' plotReducedDim(GlobalPatterns,"CCA", colour_by = "SampleType")
#' 
#' # Fetch significance results
#' attr(reducedDim(GlobalPatterns, "CCA"), "significance")$summary
#'
#' GlobalPatterns <- runRDA(GlobalPatterns, data ~ SampleType)
#' plotReducedDim(GlobalPatterns,"CCA", colour_by = "SampleType")
#' 
#' # To scale values when using *RDA functions, use transformCounts(MARGIN = "features", 
#' tse <- GlobalPatterns
#' tse <- transformCounts(tse, MARGIN = "features", method = "z")
#' # Data might include taxa that do not vary. Remove those because after z-transform
#' # their value is NA
#' tse <- tse[ rowSums( is.na( assay(tse, "z") ) ) == 0, ]
#' # Calculate RDA
#' tse <- runRDA(tse, formula = data ~ SampleType, 
#'               assay.type = "z", name = "rda_scaled", na.action = na.omit)
#' # Plot
#' plotReducedDim(tse,"rda_scaled", colour_by = "SampleType")
NULL

#' @rdname runCCA
setGeneric("calculateCCA", signature = c("x"),
           function(x, ...)
               standardGeneric("calculateCCA"))

#' @rdname runCCA
setGeneric("runCCA", signature = c("x"),
           function(x, ...)
               standardGeneric("runCCA"))

#' @rdname runCCA
setGeneric("calculateRDA", signature = c("x"),
           function(x, ...)
               standardGeneric("calculateRDA"))

#' @rdname runCCA
setGeneric("runRDA", signature = c("x"),
           function(x, ...)
               standardGeneric("runRDA"))

.remove_special_functions_from_terms <- function(terms){
    names(terms) <- terms
    m <- regexec("^Condition\\(([^\\(\\)]*)\\)$|^([^\\(\\)]*)$", terms)
    m <- regmatches(terms, m)
    terms <- vapply(m,
                    function(n){
                        n <- n[seq.int(2L,length(n))]
                        n[n != ""]
                    },
                    character(1))
    terms
}

.get_dependent_var_name <- function(formula){
    # Get dependent variable from factors
    dep_var <- rownames(attr(terms(formula), "factors"))[1]
    # If it is NULL, get it from the formula
    if( is.null(dep_var) ){
        dep_var <- as.character(formula)[2]
    }
    return(dep_var)
}

#' @importFrom stats as.formula
.calculate_cca <- function(x, formula, variables, scale = TRUE, ...){
    .require_package("vegan")
    # input check
    if(!.is_a_bool(scale)){
        stop("'scale' must be TRUE or FALSE.", call. = FALSE)
    }
    #
    x <- as.matrix(t(x))
    variables <- as.data.frame(variables)
    if(ncol(variables) > 0L && !missing(formula)){
        dep_var_name <- .get_dependent_var_name(formula)
        assign(dep_var_name, x)
        # recast formula in current environment
        form <- as.formula(paste(as.character(formula)[c(2,1,3)],
                                 collapse = " "))
        cca <- vegan::cca(form, data = variables, scale = scale, ...)
        X <- cca$CCA
    } else if(ncol(variables) > 0L) {
        cca <- vegan::cca(X = x, Y = variables, scale = scale, ...)
        X <- cca$CCA
    } else {
        cca <- vegan::cca(X = x, scale = scale, ...)
        X <- cca$CA
    }
    ans <- X$u
    attr(ans, "rotation") <- X$v
    attr(ans, "eigen") <- X$eig
    attr(ans, "cca") <- cca
    ans
}

#' @export
#' @rdname runCCA
setMethod("calculateCCA", "ANY", .calculate_cca)

#' @importFrom stats terms
#' @importFrom SummarizedExperiment colData
.get_variables_from_data_and_formula <- function(x, formula){
    if(missing(formula)){
        return(NULL)
    }
    terms <- rownames(attr(terms(formula),"factors"))
    terms <- terms[terms != as.character(formula)[2L]]
    terms <- .remove_special_functions_from_terms(terms)
    if(!all(terms %in% colnames(colData(x)))){
        stop("All variables on the right hand side of 'formula' must be ",
             "present in colData(x).", call. = FALSE)
    }
    colData(x)[,terms,drop=FALSE]
}

.get_formula_from_data_and_variables <- function(x, variables){
    # Check that variables specify columns from colData
    if( all(!is.character(variables)) ){
        stop("'variables' should be a character specifying variables from ",
             "colData(x).", 
             call. = FALSE)
    }
    if(!all(variables %in% colnames(colData(x)))){
        stop("All variables must be ",
             "present in colData(x).", call. = FALSE)
    }
    # Create a formula based on variables
    formula <- as.formula(paste0("x ~ ", 
                                 paste(variables, collapse = " + ")))
    formula
}

#' @export
#' @rdname runCCA
setMethod("calculateCCA", "SummarizedExperiment",
    function(x, formula, variables, ..., test_signif = TRUE,
             assay.type = assay_name, assay_name = exprs_values, exprs_values = "counts")
    {
        # Check assay.type and get assay
        .check_assay_present(assay.type, x)
        mat <- assay(x, assay.type)
        # Check test_signif
        if( !.is_a_bool(test_signif) ){
            stop("'test_signif' must be TRUE or FALSE.", call. = FALSE)
        }

        # If formula is missing but variables are not
        if( !missing(variables) && missing(formula) ){
            # Create a formula based on variables
            formula <- .get_formula_from_data_and_variables(x, variables)
            # Get the data from colData
            variables <- colData(x)[ , variables, drop = FALSE]
        } else{
            # Otherwise if formula is provided, get variables based on formula
            # (If formula is not provided variables is just empty data.frame)
            variables <- .get_variables_from_data_and_formula(x, formula)
        } 
        cca <- .calculate_cca(mat, formula, variables, ...)
        
        # Test significance if specified
        if( test_signif ){
            res <- .test_rda(mat, attr(cca, "cca"), variables, ...)
            attr(cca, "significance") <- res
        }
        return(cca)
    }
)

#' @export
#' @rdname runCCA
#' @importFrom SingleCellExperiment reducedDim<-
setMethod("runCCA", "SingleCellExperiment",
    function(x, ..., altexp = NULL, name = "CCA")
    {
        if (!is.null(altexp)) {
          y <- altExp(x, altexp)
        } else {
          y <- x
        }
        # Calculate CCA
        cca <- calculateCCA(y, ...)
        # If samples do not match / there were samples without appropriate metadata
        # and they are now removed
        if( !all( colnames(x) %in% rownames(cca) ) ){
            # Take a subset
            x_sub <- x[ , rownames(cca) ]
            # Add CCA
            reducedDim(x_sub, name) <- cca
            # Add subset and original data  to MAE
            exp_list <- ExperimentList(original = x, subset = x_sub)
            x <- MultiAssayExperiment(exp_list)
            # Give a message
            message("After CCA, certain samples are removed. The result object ",
                    "is MAE which includes the original and subsetted data.")
        } else{
            # Otherwose put the result to reducedDim if all the samples are found
            reducedDim(x, name) <- cca
        }
        x
    }
)

.calculate_rda <- function(x, formula, variables, ...){
    .require_package("vegan")
    #
    # Transpose and ensure that the table is in matrix format
    x <- as.matrix(t(x))
    # If formula is missing (vegan:dbrda requires formula)
    if( missing(formula) ){
        formula <- x ~ 1  
    }
    # Get the dependent variable
    dep_var_name <- .get_dependent_var_name(formula)
    # Dependent variable is the assay x. It does not matter what is the left-side
    # of the formula; it is always the assay
    assign(dep_var_name, x)
    # recast formula in current environment
    form <- as.formula(paste(as.character(formula)[c(2,1,3)],
                             collapse = " "))
    # If variables are provided
    if( !missing(variables) ){
        # Convert into data.frame
        variables <- as.data.frame(variables)
        # Calculate RDA with variables
        rda <- vegan::dbrda(formula = form, data = variables, ...)
    } else{
        # Otherwise calculate RDA without variables
        rda <- vegan::dbrda(formula = form, ...)
    }
    # Get CCA
    X <- rda$CCA
    # If variable(s) do not explain inertia at all, CCA is NULL. Then take CA
    if( is.null(X) ){
        X <- rda$CA
    }
    ans <- X$u
    attr(ans, "rotation") <- X$v
    attr(ans, "eigen") <- X$eig
    attr(ans, "rda") <- rda
    ans
}

# Perform PERMANOVA and homogeneity analysis to RDA object
#' @importFrom vegan anova.cca betadisper
#' @importFrom stats anova
.test_rda <- function(mat, rda, variables, ...){
    # Perform permanova for whole model and for variables
    permanova_model <- anova.cca(rda, by = NULL, ...)
    if( !is.null(variables) ){
        res <- .test_rda_vars(
            mat, rda, variables, permanova_model, by = "margin", ...)
    } else{
        res <- list(permanova = permanova_model)
    }
    return(res)
}

.test_rda_vars <- function(mat, rda, variables, permanova_model, by = "margin", ...){
    permanova <- anova.cca(rda, by = by, ...)
    # Perform homogeneity analysis
    mat <- t(mat)
    # Get the dissimilarity matrix based on original dissimilarity index
    # provided by user.
    dist_mat <- vegdist(mat, ...)
    # For all variables run the analysis
    homogeneity <- lapply(colnames(variables), function(x){
        # Get variable values
        var <- variables[[x]]
        # Run betadisper.
        # Suppress possible warnings: "some squared distances are negative and changed to zero"
        # Suppress possible messages: "missing observations due to 'group' removed"
        suppressWarnings(
        suppressMessages(
        betadisper_res <- betadisper(dist_mat, group = var)
        )
        )
        # Run permanova
        anova_res <- anova( betadisper_res, ... )
        # Return the results as a list
        res <- list(betadisper = betadisper_res, anova = anova_res)
        return(res)
    })
    names(homogeneity) <- colnames(variables)
    # Create a table from the results
    # PERMANOVAs
    table_model <- as.data.frame(permanova_model)
    table <- as.data.frame(permanova)
    # Take only model results
    table <- rbind(table_model[1, ], table)
    # Take homogeneity results of the different groups
    homogeneity_tab <- lapply(homogeneity, function(x) as.data.frame(x$anova)[1, ])
    homogeneity_tab <- do.call(rbind, homogeneity_tab)
    # Adjust colnames of tables, add info
    colnames(table) <- paste0(colnames(table), " (PERMANOVA)")
    colnames(homogeneity_tab) <- paste0(colnames(homogeneity_tab), " (homogeneity)")
    # Combine and adjust rownames
    table <- merge(table, homogeneity_tab, by=0, all.x = TRUE)
    rownames(table) <- table[["Row.names"]]
    table[["Row.names"]] <- NULL
    
    # Add the results to original RDA object
    res <- list(
        summary = table,
        permanova = list(model = permanova_model, variables = permanova),
        homogeneity = homogeneity)
    return(res)
}
#' @export
#' @rdname runCCA
setMethod("calculateRDA", "ANY", .calculate_rda)

#' @export
#' @rdname runCCA
setMethod("calculateRDA", "SummarizedExperiment",
    function(x, formula, variables, ..., test_signif = TRUE,
             assay.type = assay_name, assay_name = exprs_values, exprs_values = "counts")
    {
        # Check assay.type and get assay
        .check_assay_present(assay.type, x)
        mat <- assay(x, assay.type)
        # Check test_signif
        if( !.is_a_bool(test_signif) ){
            stop("'test_signif' must be TRUE or FALSE.", call. = FALSE)
        }
        
        # If formula is missing but variables are not
        if( !missing(variables) && missing(formula) ){
            # Create a formula based on variables
            formula <- .get_formula_from_data_and_variables(x, variables)
            # Get the data from colData
            variables <- colData(x)[ , variables, drop = FALSE]
        } else{
            # Otherwise if formula is provided, get variables based on formula
            # (If formula is not provided variables is just empty data.frame)
            variables <- .get_variables_from_data_and_formula(x, formula)
        } 
        # Calculate RDA
        rda <- .calculate_rda(mat, formula, variables, ...)
        
        # Test significance if specified
        if( test_signif ){
            res <- .test_rda(mat, attr(rda, "rda"), variables, ...)
            attr(rda, "significance") <- res
        }
        return(rda)
    }
)

#' @export
#' @rdname runCCA
#' @importFrom SingleCellExperiment reducedDim<-
setMethod("runRDA", "SingleCellExperiment",
    function(x, ..., altexp = NULL, name = "RDA")
    {
        if (!is.null(altexp)) {
          y <- altExp(x, altexp)
        } else {
          y <- x
        }
        # Calculate RDA
        rda <- calculateRDA(y, ...)
        # If samples do not match / there were samples without appropriate metadata
        # and they are now removed
        if( !all( colnames(x) %in% rownames(rda) ) ){
            # Take a subset
            x_sub <- x[ , rownames(rda) ]
            # Add RDA
            reducedDim(x_sub, name) <- rda
            # Add subset and original data  to MAE
            exp_list <- ExperimentList(original = x, subset = x_sub)
            x <- MultiAssayExperiment(exp_list)
            # Give a message
            message("After RDA, certain samples are removed. The result object ",
                    "is MAE which includes the original and subsetted data.")
        } else{
            # Otherwise add the RDA to original data's reducedDim
            reducedDim(x, name) <- rda
        }
        return(x)
    }
)
