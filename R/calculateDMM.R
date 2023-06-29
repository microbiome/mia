#' Dirichlet-Multinomial Mixture Model: Machine Learning for Microbiome Data
#'
#' These functions are accessors for functions implemented in the
#' \code{\link[DirichletMultinomial:DirichletMultinomial-package]{DirichletMultinomial}}
#' package
#'
#' @param x a numeric matrix with samples as rows or a
#'   \code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#'   object.
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
#' @param k the number of Dirichlet components to fit. See
#'   \code{\link[DirichletMultinomial:dmn]{dmn}}
#'
#' @param BPPARAM A
#'   \code{\link[BiocParallel:BiocParallelParam-class]{BiocParallelParam}}
#'   object specifying whether the UniFrac calculation should be parallelized.
#'
#' @param transposed Logical scalar, is x transposed with samples in rows?
#'
#' @param type the type of measure used for the goodness of fit. One of
#'   \sQuote{laplace}, \sQuote{AIC} or \sQuote{BIC}.
#'
#' @param name the name to store the result in
#'   \code{\link[SummarizedExperiment:RangedSummarizedExperiment-class]{metadata}}
#'
#' @param variable a variable from \code{colData} to use as a grouping variable.
#'   Must be a character of factor.
#'
#' @param seed random number seed. See
#'   \code{\link[DirichletMultinomial:dmn]{dmn}}
#'
#' @param ... optional arguments not used.
#'
#' @return
#' \code{calculateDMN} and \code{getDMN} return a list of \code{DMN} objects,
#' one element for each value of k provided.
#'
#' \code{bestDMNFit} returns the index for the best fit and \code{getBestDMNFit}
#' returns a single \code{DMN} object.
#'
#' \code{calculateDMNgroup} returns a
#'   \code{\link[DirichletMultinomial:DMNGroup-class]{DMNGroup}} object
#'
#' \code{performDMNgroupCV} returns a \code{data.frame}
#'
#' @seealso
#' \code{\link[DirichletMultinomial:DMN-class]{DMN-class}},
#' \code{\link[DirichletMultinomial:DMNGroup-class]{DMNGroup-class}},
#' \code{\link[DirichletMultinomial:dmn]{dmn}},
#' \code{\link[DirichletMultinomial:dmngroup]{dmngroup}},
#' \code{\link[DirichletMultinomial:cvdmngroup]{cvdmngroup }},
#' \code{\link[DirichletMultinomial:fitted]{accessors for DMN objects}}
#'
#' @name calculateDMN
#'
#' @examples
#' fl <- system.file(package="DirichletMultinomial", "extdata", "Twins.csv")
#' counts <- as.matrix(read.csv(fl, row.names=1))
#' fl <- system.file(package="DirichletMultinomial", "extdata", "TwinStudy.t")
#' pheno0 <- scan(fl)
#' lvls <- c("Lean", "Obese", "Overwt")
#' pheno <- factor(lvls[pheno0 + 1], levels=lvls)
#' colData <- DataFrame(pheno = pheno)
#'
#' tse <- TreeSummarizedExperiment(assays = list(counts = counts),
#'                                 colData = colData)
#'
#'
#' # Compute DMM algorithm and store result in metadata
#' tse <- cluster(tse, name = "DMM", DmmParam(),
#'                MARGIN = "samples", full = TRUE)
#' 
#' # Get the list of DMN objects
#' metadata(tse)$DMM$dmm
#' # Get which objects fits best
#' metadata(tse)$DMM$best
#' # Get the model that fits best
#' metadata(tse)$DMM$dmm[[metadata(tse)$DMM$best]]
#'
NULL

#' @rdname calculateDMN
#' @export
setGeneric("calculateDMN", signature = c("x"),
           function(x, ...)
               standardGeneric("calculateDMN"))

#' @importFrom DirichletMultinomial dmn
#' @importFrom stats runif
.calculate_DMN <- function(x, k = 1, BPPARAM = SerialParam(),
                           seed = runif(1, 0, .Machine$integer.max), ...){
    if(!is.numeric(k) ||
       length(k) == 0 ||
       anyNA(k) ||
       any(k <= 0) ||
       any(k != as.integer(k))){
        stop("'k' must be an integer vector with positive values only.",
             call. = FALSE)
    }
    #
    old <- getAutoBPPARAM()
    setAutoBPPARAM(BPPARAM)
    on.exit(setAutoBPPARAM(old))
    if (!(bpisup(BPPARAM) || is(BPPARAM, "MulticoreParam"))) {
        bpstart(BPPARAM)
        on.exit(bpstop(BPPARAM), add = TRUE)
    }

    ans <- BiocParallel::bplapply(k, DirichletMultinomial::dmn, count = x,
                                  seed = seed, ...,
                                  BPPARAM = BPPARAM)
    ans
}

#' @rdname calculateDMN
#' @export
setMethod("calculateDMN", signature = c(x = "ANY"), .calculate_DMN)

#' @rdname calculateDMN
#' @export
setMethod("calculateDMN", signature = c(x = "SummarizedExperiment"),
    function(x, assay.type = assay_name, assay_name = exprs_values, exprs_values = "counts", 
             transposed = FALSE, ...){
        .Deprecated(old="calculateDMN", new="cluster", 
                    "Now calculateDMN is deprecated. Use cluster with DMMParam parameter instead.")
        mat <- assay(x, assay.type)
        if(!transposed){
            mat <- t(mat)
        }
        calculateDMN(mat, ...)
    }
)

#' @rdname calculateDMN
#' @importFrom S4Vectors metadata<-
#' @export
runDMN <- function(x, name = "DMN", ...){
    .Deprecated(old="runDMN", new="cluster", 
                "Now runDMN is deprecated. Use cluster with DMMParam parameter instead.")
    if(!is(x,"SummarizedExperiment")){
        stop("'x' must be a SummarizedExperiment")
    }
    metadata(x)[[name]] <- calculateDMN(x, ...)
    x
}

################################################################################
# accessors

#' @importFrom S4Vectors metadata
.get_dmn <- function(x, name){
    dmn <- metadata(x)[[name]]
    if(is.null(dmn)){
        stop("No data found for '",name,"'.", call. = FALSE)
    }
    all_dmn <- vapply(dmn, is, logical(1), "DMN")
    if(!all(all_dmn)){
        stop("'name' does not a list of DMN objects.", call. = FALSE)
    }
    dmn
}

.get_dmn_fit_FUN <- function(type){
    type <- match.arg(type, c("laplace","AIC","BIC"))
    fit_FUN <- switch(type,
                      laplace = DirichletMultinomial::laplace,
                      AIC = DirichletMultinomial::AIC,
                      BIC = DirichletMultinomial::BIC)
    fit_FUN
}

#' @rdname calculateDMN
#' @export
setGeneric("getDMN", signature = "x",
           function(x, name = "DMN", ...)
               standardGeneric("getDMN"))

#' @rdname calculateDMN
#' @importFrom DirichletMultinomial laplace AIC BIC
#' @export
setMethod("getDMN", signature = c(x = "SummarizedExperiment"),
    function(x, name = "DMN"){
        .Deprecated(old="getDMN", new="cluster", 
                    "Now getDMN is deprecated. Use cluster with DMMParam parameter and full parameter set as true instead.")
        .get_dmn(x, name)
    }
)


.get_best_dmn_fit <- function(dmn, fit_FUN){
    fit <- vapply(dmn, fit_FUN, numeric(1))
    which.min(fit)
}

#' @rdname calculateDMN
#' @export
setGeneric("bestDMNFit", signature = "x",
           function(x, name = "DMN", type = c("laplace","AIC","BIC"), ...)
               standardGeneric("bestDMNFit"))

#' @rdname calculateDMN
#' @importFrom DirichletMultinomial laplace AIC BIC
#' @export
setMethod("bestDMNFit", signature = c(x = "SummarizedExperiment"),
    function(x, name = "DMN", type = c("laplace","AIC","BIC")){
        .Deprecated(old="bestDMNFit", new="cluster", 
                    "Now bestDMNFit is deprecated. Use cluster with DMMParam parameter and full parameter set as true instead.")
        #
        dmn <- getDMN(x, name)
        fit_FUN <- .get_dmn_fit_FUN(type)
        #
        .get_best_dmn_fit(dmn, fit_FUN)
    }
)

#' @rdname calculateDMN
#' @export
setGeneric("getBestDMNFit", signature = "x",
           function(x, name = "DMN", type = c("laplace","AIC","BIC"), ...)
               standardGeneric("getBestDMNFit"))

#' @rdname calculateDMN
#' @importFrom DirichletMultinomial laplace AIC BIC
#' @export
setMethod("getBestDMNFit", signature = c(x = "SummarizedExperiment"),
    function(x, name = "DMN", type = c("laplace","AIC","BIC")){
        .Deprecated(old="getBestDMNFit", new="cluster", 
                    "Now getBestDMNFit is deprecated. Use cluster with DMMParam parameter and full parameter set as true instead.")
        dmn <- getDMN(x, name)
        fit_FUN <- .get_dmn_fit_FUN(type)
        dmn[[.get_best_dmn_fit(dmn, fit_FUN)]]
    }
)

################################################################################
# DMN group

#' @rdname calculateDMN
#' @export
setGeneric("calculateDMNgroup", signature = c("x"),
           function(x, ...)
               standardGeneric("calculateDMNgroup"))

#' @importFrom DirichletMultinomial dmngroup
#' @importFrom stats runif
.calculate_DMNgroup <- function(x, variable, k = 1,
                                seed = runif(1, 0, .Machine$integer.max), ...){
    # input check
    if(!is.factor(variable) && is.character(variable)){
        variable <- factor(variable, unique(variable))
    } else if(!is.factor(variable)) {
        stop("'variable' must be a factor or a character value.", call. = FALSE)
    }
    #
    variable <- droplevels(variable)
    dmngroup(x, variable, k = k, seed = seed, ...)
}

#' @rdname calculateDMN
#' @export
setMethod("calculateDMNgroup", signature = c(x = "ANY"), .calculate_DMNgroup)

#' @rdname calculateDMN
#' @export
setMethod("calculateDMNgroup", signature = c(x = "SummarizedExperiment"),
    function(x, variable, 
             assay.type = assay_name, assay_name = exprs_values, exprs_values = "counts", 
             transposed = FALSE, ...){
        mat <- assay(x, assay.type)
        if(!transposed){
            mat <- t(mat)
        }
        variable <- colData(x)[,variable]
        if(is.null(variable)){
            stop("No data found in '",variable,"' column of colData(x).",
                 call. = FALSE)
        }
        calculateDMNgroup(x = mat, variable = variable, ...)
    }
)

################################################################################
# DMN group cross validations

#' @rdname calculateDMN
#' @export
setGeneric("performDMNgroupCV", signature = c("x"),
           function(x, ...)
               standardGeneric("performDMNgroupCV"))

#' @importFrom DirichletMultinomial cvdmngroup
#' @importFrom stats runif
.perform_DMNgroup_cv <- function(x, variable, k = 1,
                                 seed = runif(1, 0, .Machine$integer.max), ...){
    # input check
    if(!is.factor(variable) && is.character(variable)){
        variable <- factor(variable, unique(variable))
    } else if(!is.factor(variable)) {
        stop("'variable' must be a factor or a character value.", call. = FALSE)
    }
    variable <- droplevels(variable)
    if(is.null(names(k)) || !all(names(k) %in% levels(variable))){
        stop("'k' must be named. Names must fit the levels of 'variable'.",
             call. = FALSE)
    }
    #
    cvdmngroup(nrow(x), x, variable, k = k, seed = seed, ...)
}

#' @rdname calculateDMN
#' @export
setMethod("performDMNgroupCV", signature = c(x = "ANY"), .perform_DMNgroup_cv)

#' @rdname calculateDMN
#' @export
setMethod("performDMNgroupCV", signature = c(x = "SummarizedExperiment"),
    function(x, variable, 
             assay.type = assay_name, assay_name = exprs_values, exprs_values = "counts", 
             transposed = FALSE, ...){
        mat <- assay(x, assay.type)
        if(!transposed){
            mat <- t(mat)
        }
        variable <- colData(x)[,variable]
        if(is.null(variable)){
            stop("No data found in '",variable,"' column of colData(x).",
                 call. = FALSE)
        }
        performDMNgroupCV(x = mat, variable = variable, ...)
    }
)
