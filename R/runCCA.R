#' Canonical Correspondance Analysis
#'
#' These functions perform Canonical Correspondance Analysis on data stored
#' in a \code{SummarizedExperiment}.
#'
#' @param x For \code{calculateCCA} a numeric matrix with columns as samples or
#'   a
#'   \code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}.
#'
#'   For \code{runCCA} a
#'   \code{\link[SingleCellExperiment:SingleCellExperiment]{SingleCellExperiment}}
#'   or a derived object.
#'
#' @param formula If \code{x} is a
#'   \code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#'   a formula can be supplied. Based on the right-hand side of the given formula
#'   \code{colData} is subset to \code{variables}.
#'
#' @param variables a \code{data.frame} or an object coercible to one containing
#'   the variables to use. Can be missing, which turns the CCA analysis into
#'   a CA analysis. All variables are used. Please subset, if you want to
#'   consider only some of them.
#'
#' @param scale Logical scalar, should the expression values be standardized?
#' 
#' @param abund_values a single \code{character} value for specifying which
#'   assay to use for calculation.
#'
#' @param exprs_values a single \code{character} value for specifying which
#'   assay to use for calculation.
#'   (Please use \code{abund_values} instead. At some point \code{exprs_values}
#'   will be disabled.)
#'
#' @param altexp String or integer scalar specifying an alternative experiment
#'   containing the input data.
#'
#' @param name String specifying the name to be used to store the result in the
#'   reducedDims of the output.
#'
#' @param ... additional arguments not used.
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
#' and \code{\link[vegan:cca]{rda}}
#'
#' @examples
#' library(scater)
#' data(GlobalPatterns)
#' GlobalPatterns <- runCCA(GlobalPatterns, data ~ SampleType)
#' plotReducedDim(GlobalPatterns,"CCA", colour_by = "SampleType")
#'
#' GlobalPatterns <- runRDA(GlobalPatterns, data ~ SampleType)
#' plotReducedDim(GlobalPatterns,"CCA", colour_by = "SampleType")
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
    rownames(attr(terms(formula), "factors"))[1]
}

#' @importFrom stats as.formula
.calculate_cca <- function(x, formula, variables, scale = TRUE){
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
        cca <- vegan::cca(form, data = variables, scale = scale)
        X <- cca$CCA
    } else if(ncol(variables) > 0L) {
        cca <- vegan::cca(X = x, Y = variables, scale = scale)
        X <- cca$CCA
    } else {
        cca <- vegan::cca(X = x, scale = scale)
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

#' @export
#' @rdname runCCA
setMethod("calculateCCA", "SummarizedExperiment",
    function(x, formula, ..., abund_values = exprs_values, exprs_values = "counts")
    {
        mat <- assay(x, abund_values)
        variables <- .get_variables_from_data_and_formula(x, formula)
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
        reducedDim(x, name) <- calculateCCA(y, ...)
        x
    }
)

.calculate_rda <- function(x, formula, variables, scale = TRUE, ...){
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
        rda <- vegan::rda(form, data = variables, scale = scale, ...)
        X <- rda$CCA
    } else if(ncol(variables) > 0L) {
        rda <- vegan::rda(X = x, Y = variables, scale = scale, ...)
        X <- rda$CCA
    } else {
        rda <- vegan::rda(X = x, scale = scale, ...)
        X <- rda $CA
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
    function(x, formula, ..., abund_values = exprs_values, exprs_values = "counts")
    {
        mat <- assay(x, abund_values)
        variables <- .get_variables_from_data_and_formula(x, formula)
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
        reducedDim(x, name) <- calculateRDA(y, ...)
        x
    }
)
