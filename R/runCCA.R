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
#' @param variables When \code{x} is a \code{SummarizedExperiment},
#'   \code{variables} can be used to specify variables from \code{colData}. 
#'   
#'   When \code{x} is a matrix, \code{variables} is a \code{data.frame} or 
#'   an object coercible to one containing the variables to use. 
#'   
#'   Can be missing, which turns the CCA analysis into a CA analysis. 
#'   All variables are used. Please subset, if you want to consider only some of them. 
#'   
#' @param scale a logical scalar, should the expression values be standardized?
#'   \code{scale} is disabled when using \code{*RDA} functions. Please scale before
#'   performing RDA (Check examples.) 
#' 
#' @param assay_name a single \code{character} value for specifying which
#'   assay to use for calculation.
#'
#' @param exprs_values a single \code{character} value for specifying which
#'   assay to use for calculation.
#'   (Please use \code{assay_name} instead.)
#'   
#' @param abund_values a single \code{character} value for specifying which
#'   assay to use for calculation.
#'   (Please use \code{assay_name} instead. At some point \code{abund_values}
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
#' GlobalPatterns <- runRDA(GlobalPatterns, data ~ SampleType)
#' plotReducedDim(GlobalPatterns,"CCA", colour_by = "SampleType")
#' 
#' # To scale values when using *RDA functions, use transformFeatures
#' tse <- GlobalPatterns
#' tse <- transformFeatures(tse, method = "z")
#' # Data might include taxa that do not vary. Remove those because after z-transform
#' # their value is NA
#' tse <- tse[ rowSums( is.na( assay(tse, "z") ) ) == 0, ]
#' # Calculate RDA
#' tse <- runRDA(tse, formula = data ~ SampleType, 
#'               assay_name = "z", name = "rda_scaled", na.action = na.omit)
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
    formula <- as.formula(paste0("data ~ ", 
                                 paste(variables, collapse = " + ")))
    formula
}

#' @export
#' @rdname runCCA
setMethod("calculateCCA", "SummarizedExperiment",
    function(x, formula, variables, ..., 
             assay_name = abund_values, abund_values = exprs_values, exprs_values = "counts")
    {
        mat <- assay(x, assay_name)
        # If formula is missing but variables are not
        if( !missing(formula) ){
            variables <- .get_variables_from_data_and_formula(x, formula)
        } else if( !missing(variables) ){
            # If variables are not missing
            formula <- .get_formula_from_data_and_variables(x, variables)
            variables <- colData(x)[ , variables, drop = FALSE]
        }
        .calculate_cca(mat, formula, variables, ...)
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
        if( all(rownames(cca) != colnames(x)) ){
            # Take a subset
            x_sub <- x[ , rownames(cca) ]
            # Add CCA
            reducedDim(x_sub, name) <- cca
            # Add subset to altExp
            altExp(x, name) <- x_sub
            # Give a message
            message("After CCA, certain samples are removed. Subsetted object with ",
                    "results of CCA analysis is stored in altExp.")
        } else{
            reducedDim(x, name) <- cca
        }
        x
    }
)

.calculate_rda <- function(x, formula, variables, ...){
    .require_package("vegan")
    #
    x <- as.matrix(t(x))
    # If formula is missing
    if( missing(formula) ){
        formula <- x ~ 1  
    }
    dep_var_name <- .get_dependent_var_name(formula)
    assign(dep_var_name, x)
    # recast formula in current environment
    form <- as.formula(paste(as.character(formula)[c(2,1,3)],
                             collapse = " "))
    # If variables are not provided, formula is 'x ~ 1'. 
    if( !missing(variables) ){
        variables <- as.data.frame(variables)
        rda <- vegan::dbrda(formula = form, data = variables, ...)
    } else{
        rda <- vegan::dbrda(formula = form, ...)
    }
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

#' @export
#' @rdname runCCA
setMethod("calculateRDA", "ANY", .calculate_rda)

#' @export
#' @rdname runCCA
setMethod("calculateRDA", "SummarizedExperiment",
    function(x, formula, variables, ..., 
             assay_name = abund_values, abund_values = exprs_values, exprs_values = "counts")
    {
        .check_assay_present(assay_name, x)
        mat <- assay(x, assay_name)
        # If formula is provided, it takes the precedence
        if( !missing(formula) ){
            variables <- .get_variables_from_data_and_formula(x, formula)
        } else if( !missing(variables) ){
            # If variables are provided, create a formula based on them
            formula <- .get_formula_from_data_and_variables(x, variables)
            variables <- colData(x)[ , variables, drop = FALSE]
        }
        .calculate_rda(mat, formula, variables, ...)
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
        if( all(rownames(rda) != colnames(x)) ){
            # Take a subset
            x_sub <- x[ , rownames(rda) ]
            # Add RDA
            reducedDim(x_sub, name) <- rda
            # Add subset to altExp
            altExp(x, name) <- x_sub
            # Give a message
            message("After RDA, certain samples are removed. Subsetted object with ",
                    "results of RDA analysis is stored in altExp.")
        } else{
            reducedDim(x, name) <- rda
        }
        x
    }
)
