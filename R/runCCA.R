#' Canonical Correspondance Analysis
#' 
#' These functions perform Canonical Correspondance Analysis on data stored
#' in a \code{SummarizedExperiment}.
#' 
#' @param x a numeric matrix or a
#'   \code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}.
#' 
#' @param variables a \code{data.frame} or an object coercible to one containing
#'   the variables to use. Can be missing, which turns the CCA analysis into
#'   a CA analysis. All variables are used. Please subset, if you want to 
#'   consider only some of them.
#'
#' @param ntop Numeric scalar specifying the number of features with the highest
#'   variances to use for dimensionality reduction.
#'
#' @param subset_row Vector specifying the subset of features to use for
#'   dimensionality reduction. This can be a character vector of row names, an
#'   integer vector of row indices or a logical vector.
#'
#' @param scale Logical scalar, should the expression values be standardized?
#'
#' @param transposed Logical scalar, is x transposed with cells in rows?
#' 
#' @param dimred String or integer scalar specifying the existing dimensionality
#'   reduction results to use.
#'
#' @param n_dimred Integer scalar or vector specifying the dimensions to use if
#'   dimred is specified.
#'
#' @param altexp String or integer scalar specifying an alternative experiment
#'   containing the input data.
#'
#' @param name String specifying the name to be used to store the result in the
#'   reducedDims of the output.
#' 
#' @param formula If \code{x} is a \code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#'   a formula can be supplied. Based on the right-hand side of the given formula
#'   \code{colData} is subset to \code{variables}.
#'
#' @return 
#'
#' @name runCCA
#' @seealso
#' For more details on the actual implementation see \code{\link[vegan:cca]{cca}}
#' and \code{\link[vegan:cca]{rda}}
#' 
#' 
#' @examples
#' data("GlobalPatterns",package="MicrobiomeExperiment")
#' GlobalPatterns <- runCCA(GlobalPatterns, data ~ SampleType)
#' plotReducedDim(GlobalPatterns,"CCA", colour_by = "SampleType")
#' 
#' GlobalPatterns <- runRDA(GlobalPatterns, data ~ SampleType)
#' plotReducedDim(GlobalPatterns,"CCA", colour_by = "SampleType")
NULL

setGeneric("calculateCCA", signature = c("x"),
           function(x, ...)
               standardGeneric("calculateCCA"))

setGeneric("runCCA", signature = c("x"),
           function(x, ...)
               standardGeneric("runCCA"))

setGeneric("calculateRDA", signature = c("x"),
           function(x, ...)
               standardGeneric("calculateRDA"))

setGeneric("runRDA", signature = c("x"),
           function(x, ...)
               standardGeneric("runRDA"))

.check_dependend_var_name <- function(names, dep_var_name = "data"){
    if(dep_var_name %in% cn){
        stop("'",dep_var_name,"' cannot be the name of a independent variable.",
             " Please rename the variable.", call. = FALSE)
    }
}

#' @importFrom stats terms
#' @importFrom SummarizedExperiment colData
.get_variables_from_data_and_formula <- function(x, formula){
    if(missing(formula)){
        return(NULL)
    }
    cn <- rownames(attr(terms(formula),"factors"))
    cn <- cn[cn != as.character(formula)[2L]]
    if(!all(cn %in% colnames(colData(x)))){
        stop("All variables on the right hand side of 'formula' must be ",
             "present in colData(x).", call. = FALSE)
    }
    .check_dependend_var_name(cn, "data")
    colData(x)[,cn,drop=FALSE]
}

#' @importFrom stats as.formula
#' @importFrom vegan cca
.calculate_cca <- function(x, variables, ntop = nrow(x),
                           subset_row = NULL, scale = TRUE, transposed = FALSE)
{
    # input check
    if(!.is_a_bool(scale)){
        stop("'scale' must be TRUE or FALSE.", call. = FALSE)
    }
    #
    if(!transposed) {
        x <- .get_mat_for_reddim(x, subset_row = subset_row, ntop = ntop,
                                 scale = scale)
    }
    data <- as.matrix(x)
    variables <- as.data.frame(variables)
    .check_dependend_var_name(colnames(variables), "data")
    if(ncol(variables) > 0L){
        form <- as.formula(paste0("data ~ ",
                                  paste(colnames(variables), collapse = " + ")))
        cca <- vegan::cca(form, data = variables, X = x, scale = scale)
        X <- cca$CCA
    } else {
        cca <- vegan::cca(X = data, scale = scale)
        X <- cca$CA
    }
    ans <- X$u
    attr(ans, "rotation") <- X$v
    attr(ans, "eigen") <- X$eig
    attr(ans, "cca") <- cca
    ans
}

#' @export
#' @rdname runDPCoA
setMethod("calculateCCA", "ANY", .calculate_cca)

#' @export
#' @rdname runCCA
setMethod("calculateCCA", "SummarizedExperiment",
    function(x, formula, ..., exprs_values = "counts")
    {
        mat <- assay(x,exprs_values)
        variables <- .get_variables_from_data_and_formula(x, formula)
        .calculate_cca(mat, variables, transposed = !is.null(dimred), ...)
    }
)

#' @export
#' @rdname runCCA
setMethod("calculateCCA", "SingleCellExperiment",
    function(x, formula, ..., exprs_values = "counts", dimred = NULL,
             n_dimred = NULL)
    {
        mat <- .get_mat_from_sce(x, exprs_values = exprs_values,
                               dimred = dimred, n_dimred = n_dimred)
        variables <- .get_variables_from_data_and_formula(x, formula)
        .calculate_cca(mat, variables, transposed = !is.null(dimred), ...)
    }
)

#' @export
#' @rdname runCCA
#' @importFrom SingleCellExperiment reducedDim<-
setMethod("runCCA", "SingleCellExperiment",
    function(x, ..., altexp=NULL, name="CCA")
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

#' @importFrom vegan rda
.calculate_rda <- function(x, variables, ntop = 500,
                           subset_row = NULL, scale = TRUE, transposed = FALSE){
    # input check
    if(!.is_a_bool(scale)){
        stop("'scale' must be TRUE or FALSE.", call. = FALSE)
    }
    #
    if(!transposed) {
        x <- .get_mat_for_reddim(x, subset_row = subset_row, ntop = ntop,
                                 scale = scale)
    }
    data <- as.matrix(x)
    variables <- as.data.frame(variables)
    .check_dependend_var_name(colnames(variables), "data")
    if(ncol(variables) > 0L){
        form <- as.formula(paste0("data ~ ",
                                  paste(colnames(variables), collapse = " + ")))
        rda <- vegan::rda(form, data = variables, X = x, scale = scale)
        X <- rda$CCA
    } else {
        rda <- vegan::rda(X = data, scale = scale)
        X <- rda$CA
    }
    ans <- X$u
    attr(ans, "rotation") <- X$v
    attr(ans, "eigen") <- X$eig
    attr(ans, "rda") <- rda
    ans
}

#' @export
#' @rdname runDPCoA
setMethod("calculateRDA", "ANY", .calculate_rda)



#' @export
#' @rdname runCCA
setMethod("calculateRDA", "SummarizedExperiment",
    function(x, formula, ..., exprs_values = "counts")
    {
        mat <- assay(x, exprs_values)
        variables <- .get_variables_from_data_and_formula(x, formula)
        .calculate_rda(mat, variables, transposed = !is.null(dimred), ...)
    }
)

#' @export
#' @rdname runCCA
setMethod("calculateRDA", "SingleCellExperiment",
    function(x, formula, ..., exprs_values = "counts", dimred = NULL,
             n_dimred = NULL)
    {
        mat <- .get_mat_from_sce(x, exprs_values = exprs_values,
                               dimred = dimred, n_dimred = n_dimred)
        variables <- .get_variables_from_data_and_formula(x, formula)
        .calculate_rda(mat, variables, transposed = !is.null(dimred), ...)
    }
)

#' @export
#' @rdname runCCA
#' @importFrom SingleCellExperiment reducedDim<-
setMethod("runRDA", "SingleCellExperiment",
    function(x, ..., altexp=NULL, name="RDA")
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