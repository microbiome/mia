#' Canonical Correspondence Analysis and Redundancy Analysis
#'
#' These functions perform Canonical Correspondence Analysis on data stored
#' in a \code{SummarizedExperiment}.
#'
#' @inheritParams getDominant
#' @inheritParams getDissimilarity
#'
#' @details
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
#' @param variables \code{Character scalar}. When \code{x} is a \code{SummarizedExperiment},
#'   \code{variables} can be used to specify variables from \code{colData}. 
#'   
#'   When \code{x} is a matrix, \code{variables} is a \code{data.frame} or 
#'   an object coercible to one containing the variables to use. 
#'    
#'   All variables are used. Please subset, if you want to consider only some of them. 
#'   \code{variables} and \code{formula} can be missing, which turns the CCA analysis 
#'   into a CA analysis and dbRDA into PCoA/MDS.
#' 
#' @param test.signif \code{Logical scalar}. Should the PERMANOVA and analysis of
#'   multivariate homogeneity of group dispersions be performed.
#'   (Default: \code{TRUE})
#'   
#' @param altexp \code{Character scalar} or \code{integer scalar}. Specifies an alternative experiment
#'   containing the input data.
#' 
#' @param name \code{Character scalar}. A name for the column of the 
#'   \code{colData} where results will be stored. (Default: \code{"CCA"})
#' 
#' @param scores \code{Character scalar}. Specifies scores to be returned. Must be
#' 'wa' (site scores found as weighted averages (cca) or weighted sums (rda) of
#' v with weights Xbar, but the multiplying effect of eigenvalues removed) or
#' 'u' ((weighted) orthonormal site scores). (Default: \code{'wa'})
#' 
#' @param exprs_values Deprecated. Use \code{assay.type} instead.
#'
#' @param ... additional arguments passed to vegan::cca or vegan::dbrda and
#' other internal functions.
#' \itemize{
#'   \item{\code{method} a dissimilarity measure to be applied in dbRDA and
#'   possible following homogeneity test. (By default: \code{method="euclidean"})}
#'   \item{\code{scale}: \code{Logical scalar}. Should the expression values be
#'   standardized? \code{scale} is disabled when using \code{*RDA} functions.
#'   Please scale before performing RDA. (Default: \code{TRUE})}
#'   \item{\code{full} \code{Logical scalar}. should all the results from the
#'   significance calculations be returned. When \code{full=FALSE}, only
#'   summary tables are returned. (Default: \code{FALSE})}
#'   \item{\code{homogeneity.test}: \code{Character scalar}. Specifies
#'   the significance test used to analyse \code{vegan::betadisper} results.
#'   Options include 'permanova' (\code{vegan::permutest}), 'anova'
#'   (\code{stats::anova}) and 'tukeyhsd' (\code{stats::TukeyHSD}).
#'   (By default: \code{homogeneity.test="permanova"})}
#'   \item{\code{permutations} a numeric value specifying the number of permutations 
#'   for significance testing in \code{vegan::anova.cca}. (By default: \code{permutations=999})}
#' }
#' 
#' @details
#'   *CCA functions utilize \code{vegan:cca} and *RDA functions \code{vegan:dbRDA}.
#'   By default dbRDA is done with euclidean distances which equals to RDA.
#'   
#'   Significance tests are done with \code{vegan:anova.cca} (PERMANOVA). Group
#'   dispersion, i.e., homogeneity within groups is analyzed with 
#'   \code{vegan:betadisper} (multivariate homogeneity of groups dispersions (variances))
#'   and statistical significance of homogeneity is tested with a test
#'   specified by \code{homogeneity.test} parameter.
#'
#' @return
#' For \code{getCCA} a matrix with samples as rows and CCA dimensions as
#' columns. Attributes include calculated \code{cca}/\code{rda} object and
#' significance analysis results.
#'
#' For \code{addCCA} a modified \code{x} with the results stored in
#' \code{reducedDim} as the given \code{name}.
#'
#' @name runCCA
#' @seealso
#' For more details on the actual implementation see \code{\link[vegan:cca]{cca}}
#' and \code{\link[vegan:dbrda]{dbrda}}
#'
#' @examples
#' library(scater)
#' data(GlobalPatterns)
#' GlobalPatterns <- addCCA(GlobalPatterns, data ~ SampleType)
#' plotReducedDim(GlobalPatterns,"CCA", colour_by = "SampleType")
#' 
#' # Fetch significance results
#' attr(reducedDim(GlobalPatterns, "CCA"), "significance")
#'
#' GlobalPatterns <- addRDA(GlobalPatterns, data ~ SampleType)
#' plotReducedDim(GlobalPatterns,"CCA", colour_by = "SampleType")
#' 
#' # Specify dissimilarity measure
#' GlobalPatterns <- transformAssay(GlobalPatterns, method = "relabundance")
#' GlobalPatterns <- addRDA(
#'     GlobalPatterns, data ~ SampleType, assay.type = "relabundance",
#'     method = "bray")
#' 
#' # To scale values when using *RDA functions, use
#' # transformAssay(MARGIN = "features", ...) 
#' tse <- GlobalPatterns
#' tse <- transformAssay(tse, MARGIN = "features", method = "standardize")
#' # Data might include taxa that do not vary. Remove those because after
#' # z-transform their value is NA
#' tse <- tse[ rowSums( is.na( assay(tse, "standardize") ) ) == 0, ]
#' # Calculate RDA
#' tse <- addRDA(
#'     tse, formula = data ~ SampleType, assay.type = "standardize",
#'     name = "rda_scaled", na.action = na.omit)
#' # Plot
#' plotReducedDim(tse,"rda_scaled", colour_by = "SampleType")
#' # A common choice along with PERMANOVA is ANOVA when statistical significance
#' # of homogeneity of groups is analysed. Moreover, full significance test results
#' # can be returned.
#'  tse <- addRDA(
#'      tse, data ~ SampleType, homogeneity.test = "anova", full = TRUE)
#' # Example showing how to pass extra parameters, such as 'permutations', to anova.cca
#' tse <- addRDA(tse, data ~ SampleType, permutations = 500)
#' 
NULL

#' @rdname runCCA
setGeneric("getCCA", signature = c("x"),
           function(x, ...)
               standardGeneric("getCCA"))

#' @rdname runCCA
setGeneric("addCCA", signature = c("x"),
           function(x, ...)
               standardGeneric("addCCA"))

#' @rdname runCCA
setGeneric("getRDA", signature = c("x"),
           function(x, ...)
               standardGeneric("getRDA"))

#' @rdname runCCA
setGeneric("addRDA", signature = c("x"),
           function(x, ...)
               standardGeneric("addRDA"))

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
.calculate_cca <- function(x, formula, variables, scores,  scale = TRUE, ...){
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
    # Create the matrix to be returned
    ans <- X[[scores]]
    attr(ans, "rotation") <- X$v
    attr(ans, "eigen") <- X$eig
    attr(ans, "cca") <- cca
    ans
}

#' @export
#' @rdname runCCA
setMethod("getCCA", "ANY",
      function(x, ...){
          .calculate_cca(x, ...)
      })

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
setMethod("getCCA", "SummarizedExperiment",
    function(x, formula, variables, test.signif = TRUE,
             assay.type = assay_name, assay_name = exprs_values, exprs_values = "counts",
             scores = "wa", ...)
    {
        # Check assay.type and get assay
        .check_assay_present(assay.type, x)
        mat <- assay(x, assay.type)
        # Check test.signif
        if( !.is_a_bool(test.signif) ){
            stop("'test.signif' must be TRUE or FALSE.", call. = FALSE)
        }
        if( !(.is_a_string(scores) && scores %in% c("wa", "u")) ){
            stop("'scores' must be 'wa' or 'u'.",
                 call. = FALSE)
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
        cca <- .calculate_cca(mat, formula, variables, scores, ...)
        
        # Test significance if specified
        if( test.signif ){
            res <- .test_rda(mat, attr(cca, "cca"), variables, ...)
            attr(cca, "significance") <- res
        }
        return(cca)
    }
)

#' @export
#' @rdname runCCA
#' @aliases getCCA
calculateCCA <- function(x,...){
    getCCA(x,...)
}

#' @export
#' @rdname runCCA
#' @importFrom SingleCellExperiment reducedDim<-
setMethod("addCCA", "SingleCellExperiment",
    function(x, formula, variables, altexp = NULL, name = "CCA", ...)
    {
        if (!is.null(altexp)) {
          y <- altExp(x, altexp)
        } else {
          y <- x
        }
        # Calculate CCA
        cca <- getCCA(y, formula, variables, ...)
        # Add object to reducedDim
        x <- .add_object_to_reduceddim(x, cca, name = name, ...)
        return(x)
    }
)

#' @export
#' @rdname runCCA
#' @aliases addCCA
runCCA <- function(x,...){
    addCCA(x,...)
}

#' @importFrom vegan sppscores<-
.calculate_rda <- function(
        x, formula, variables, scores, method = distance, distance = "euclidean", ...){
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
        rda <- vegan::dbrda(formula = form, data = variables, distance = method, ...)
    } else{
        # Otherwise calculate RDA without variables
        rda <- vegan::dbrda(formula = form, distance = method, ...)
    }
    # Get CCA
    if( !is.null(rda$CCA) ){
        X <- rda$CCA
        # Get species scores. Get only those samples that are included in rda
        # object (some might missing due missing metadata)
        species_scores <- x[ rownames(X[[scores]]),  ]
    } else{
        # If variable(s) do not explain inertia at all, CCA is NULL. Then take CA
        X <- rda$CA
        # Get species scores (whole data since metadata was not in input)
        species_scores <- x
        # If scores is "wa", but they are not available
        if( scores == "wa" ){
            warning("'wa' scores are not available. Defaults to 'u'.", call. = FALSE)
            scores <- "u"
        }
    }
    # Add species scores since they are missing from dbrda object (in cca they are included)
    sppscores(rda) <- species_scores
    # Create the matrix to be returned
    ans <- X[[scores]]
    attr(ans, "rotation") <- X$v
    attr(ans, "eigen") <- X$eig
    attr(ans, "rda") <- rda
    ans
}

# Add RDA/CCA to reducedDim
.add_object_to_reduceddim <- function(
        tse, rda, name, subset_result = FALSE, ...){
    # Test subset
    if( !.is_a_bool(subset_result) ){
        stop("'subset_result' must be TRUE or FALSE.", call. = FALSE)
    }
    #
    # If samples do not match / there were samples without appropriate metadata
    # and they are now removed
    if( !all(colnames(tse) %in% rownames(rda)) && subset_result ){
        # Take a subset
        tse <- tse[ , rownames(rda) ]
        # Give a message
        message("Certain samples are removed from the result because they did ",
                "not include sufficient metadata.")
    } else if( !all(colnames(tse) %in% rownames(rda)) && !subset_result ){
        # If user do not want to subset the data
        # Save attributes from the object
        attr <- attributes(rda)
        attr <- attr[ !names(attr) %in% c("dim", "dimnames")]
        # Find samples that are removed
        samples_not_found <- setdiff(colnames(tse), rownames(rda))
        # Create an empty matrix
        mat <- matrix(nrow = length(samples_not_found), ncol=ncol(rda))
        rownames(mat) <- samples_not_found
        # Combine the data and order it in correct order
        rda <- rbind(rda, mat)
        rda <- rda[colnames(tse), ]
        # Add attributes
        attr <- c(attributes(rda), attr)
        attributes(rda) <- attr
    }
    # Add object to reducedDIm
    reducedDim(tse, name) <- rda
    return(tse)
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

# Test association of variables to ordination
.test_rda_vars <- function(
        mat, rda, variables, permanova_model, by = "margin", full = FALSE,
        homogeneity.test = "permanova", method = distance, distance = "euclidean",
        ...){
    # Check full parameter
    if( !.is_a_bool(full) ){
        stop("'full' must be TRUE or FALSE.", call. = FALSE)
    }
    # Check homogeneity.test
    if( !(is.character(homogeneity.test) &&
          length(homogeneity.test) == 1 &&
          homogeneity.test %in% c("permanova", "anova", "tukeyhsd")) ){
        stop("'homogeneity.test' must be 'permanova', 'anova', or 'tukeyhsd'.",
             call. = FALSE)
    }
    #
    # Perform PERMANOVA
    permanova <- anova.cca(rda, by = by, ...)
    # Create a table from the results
    # PERMANOVAs
    table_model <- as.data.frame(permanova_model)
    permanova_tab <- as.data.frame(permanova)
    # Take only model results
    permanova_tab <- rbind(table_model[1, ], permanova_tab)
    # Add info about total variance
    permanova_tab[ , "Total variance"] <- permanova_tab["Model", 2] + permanova_tab["Residual", 2]
    # Add info about explained variance
    permanova_tab[ , "Explained variance"] <- permanova_tab[ , 2] / permanova_tab[ , "Total variance"]

    # Perform homogeneity analysis
    mat <- t(mat)
    # Get the dissimilarity matrix based on original dissimilarity index
    # provided by user. If the analysis is CCA, disable method; calculate
    # always euclidean distances because CCA is based on Euclidean distances.
    if( length(class(rda)) == 1 && is(rda, 'cca') ){
        dist_mat <- vegdist(mat, method = "euclidean")
    } else{
        dist_mat <- vegdist(mat, method = method, ...)
    }
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
        # Run significance test
        significance <- .homogeneity_significance(betadisper_res, homogeneity.test, ...)
        # Return the results as a list
        models <- list(betadisper_res, significance[["obj"]])
        names(models) <- c("betadisper", homogeneity.test)
        res <- list(
            models = models,
            table = significance[["table"]]
        )
        return(res)
    })
    names(homogeneity) <- colnames(variables)
    # Get models
    homogeneity_model <- lapply(homogeneity, function(x) x[["models"]])
    names(homogeneity_model) <- colnames(variables)
    # Get tables
    homogeneity_tab <- lapply(homogeneity, function(x) x[["table"]])
    # Combine tables
    homogeneity_tab <- do.call(rbind, homogeneity_tab)
    
    # Return whole data or just a tables
    if( full ){
        res <- list(
            permanova = list(
                summary = permanova_tab,
                model = permanova_model,
                variables = permanova),
            homogeneity = list(
                summary = homogeneity_tab,
                variables = homogeneity_model)
        )
    } else{
        res <- list(
            permanova = permanova_tab,
            homogeneity = homogeneity_tab)
    }
    return(res)
}

# Perform statistical test for group homogeneity results
.homogeneity_significance <- function(betadisper_res, homogeneity.test, ...){
    # Run specified significance test
    if( homogeneity.test == "anova" ){
        res <- stats::anova(betadisper_res, ...)
    } else if ( homogeneity.test == "tukeyhsd" ){
        res <- stats::TukeyHSD(betadisper_res, ...)
    } else{
        res <- vegan::permutest(betadisper_res, ...)
    }

    # Get summary table from the results
    if( homogeneity.test == "anova" ){
        tab <- as.data.frame(res)
    } else if( homogeneity.test == "tukeyhsd" ){
        tab <- res[["group"]]
    } else{
        tab <- res[["tab"]]
    }
    
    # Modify permanova/anova table
    if( homogeneity.test != "tukeyhsd" ){
        # Add info about total variance
        tab[ , "Total variance"] <- tab["Groups", "Sum Sq"] + tab["Residuals", "Sum Sq"]
        # Add info about explained variance
        tab[ , "Explained variance"] <- tab[ , "Sum Sq"] / tab[ , "Total variance"]
        # Get only groups row (drop residuals row)
        tab <- tab[1, ]
    }
    
    # Create list from the object and results
    res <- list(
        obj = res,
        table = tab
    )
    return(res)
}

#' @export
#' @rdname runCCA
setMethod("getRDA", "ANY",
      function(x, ...){
          .calculate_rda(x, ...)
      })

#' @export
#' @rdname runCCA
setMethod("getRDA", "SummarizedExperiment",
    function(x, formula, variables, test.signif = TRUE,
             assay.type = assay_name, assay_name = exprs_values, exprs_values = "counts",
             scores = "wa", ...)
    {
        # Check assay.type and get assay
        .check_assay_present(assay.type, x)
        mat <- assay(x, assay.type)
        # Check test.signif
        if( !.is_a_bool(test.signif) ){
            stop("'test.signif' must be TRUE or FALSE.", call. = FALSE)
        }
        if( !(.is_a_string(scores) && scores %in% c("wa", "u", "v")) ){
            stop("'scores' must be 'wa', 'u', or 'v'.",
                 call. = FALSE)
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
        rda <- .calculate_rda(mat, formula, variables, scores, ...)
        
        # Test significance if specified
        if( test.signif ){
            res <- .test_rda(mat, attr(rda, "rda"), variables, ...)
            attr(rda, "significance") <- res
        }
        return(rda)
    }
)

#' @export
#' @rdname runCCA
#' @aliases getRDA
calculateRDA <- function(x,...){
    getRDA(x,...)
}

#' @export
#' @rdname runCCA
#' @importFrom SingleCellExperiment reducedDim<-
setMethod("addRDA", "SingleCellExperiment",
    function(x, formula, variables, altexp = NULL, name = "RDA", ...)
    {
        if (!is.null(altexp)) {
          y <- altExp(x, altexp)
        } else {
          y <- x
        }
        # Calculate RDA
        rda <- getRDA(y, formula, variables, ...)
        # Add object to reducedDim
        x <- .add_object_to_reduceddim(x, rda, name = name, ...)
        return(x)
    }
)

#' @export
#' @rdname runCCA
#' @aliases addRDA
runRDA <- function(x,...){
    addRDA(x,...)
}
