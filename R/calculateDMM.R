#' Dirichlet-Multinomial Mixture Model: Machine Learning for Microbiome Data
#'
#' These functions are accessors for functions implemented in the
#' \code{\link[DirichletMultinomial:DirichletMultinomial-package]{DirichletMultinomial}}
#' package
#'
#' @param x a numeric matrix with samples as rows or a
#' \code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#' object.
#'
#' @param assay.type \code{Character scalar}. Specifies the name of the
#'   assay used in calculation. (Default: \code{"counts"})
#'
#' @param k \code{Numeric scalar}. The number of Dirichlet components to fit.
#' See \code{\link[DirichletMultinomial:dmn]{dmn}}. (Default: \code{1})
#'
#' @param transposed \code{Logical scalar}. Is \code{x} transposed with samples
#' in rows? (Default: \code{FALSE})
#'
#' @param variable \code{Character scalar}. A variable from \code{colData} to
#' use as a grouping variable. Must be a character of factor.
#'
#' @param seed \code{Numeric scalar}. Random number seed. See
#'   \code{\link[DirichletMultinomial:dmn]{dmn}}
#'
#' @param ... optional arguments not used.
#'
#' @return
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
#' \code{\link[DirichletMultinomial:mixture]{accessors for DMN objects}}
#'
#' @references
#' \insertRef{Holmes2012}{mia}
#'
#' @name calculateDMNgroup
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
#' tse <- TreeSummarizedExperiment(
#'     assays = list(counts = counts),
#'     colData = colData)
#'
#' library(bluster)
#'
#' # Compute DMM algorithm and store result in metadata
#' tse <- addCluster(
#'     tse,
#'     assay.type = "counts",
#'     name = "DMM",
#'     DmmParam(k = 1:3, type = "laplace"),
#'     by = "samples",
#'     full = TRUE
#' )
#'
#' # Get the list of DMN objects
#' metadata(tse)$DMM$dmm
#'
#' # Get and display which objects fits best
#' bestFit <- metadata(tse)$DMM$best
#' bestFit
#'
#' # Get the model that generated the best fit
#' bestModel <- metadata(tse)$DMM$dmm[[bestFit]]
#' bestModel
#'
#' # Get the sample-cluster assignment probability matrix
#' head(metadata(tse)$DMM$prob)
#'
#' # Get the weight of each component for the best model
#' bestModel@mixture$Weight
NULL

# Used by miaViz::plotDMNFit()
.get_dmn_fit_FUN <- function(type){
    type <- match.arg(type, c("laplace","AIC","BIC"))
    fit_FUN <- switch(type,
                        laplace = DirichletMultinomial::laplace,
                        AIC = DirichletMultinomial::AIC,
                        BIC = DirichletMultinomial::BIC)
    fit_FUN
}

# Fits DMN models with k components to a samples x features matrix
#' @importFrom BiocParallel SerialParam bplapply bpisup bpstart bpstop
#'     getAutoBPPARAM setAutoBPPARAM
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

################################################################################
# DMN group

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

#' @rdname calculateDMNgroup
#' @export
setMethod("calculateDMNgroup", signature = c(x = "ANY"), .calculate_DMNgroup)

#' @rdname calculateDMNgroup
#' @export
setMethod("calculateDMNgroup", signature = c(x = "SummarizedExperiment"),
            function(x, variable,
                    assay.type = "counts", transposed = FALSE, ...){
                mat <- assay(x, assay.type)
                if(!transposed){
                    mat <- t(mat)
                }
                variable <- colData(x)[,variable]
                if(is.null(variable)){
                    stop("No data found in '",variable,
                        "' column of colData(x).", call. = FALSE)
                }
                calculateDMNgroup(x = mat, variable = variable, ...)
            }
)

################################################################################
# DMN group cross validations

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

#' @rdname calculateDMNgroup
#' @export
setMethod("performDMNgroupCV", signature = c(x = "ANY"), .perform_DMNgroup_cv)

#' @rdname calculateDMNgroup
#' @export
setMethod("performDMNgroupCV", signature = c(x = "SummarizedExperiment"),
            function(x, variable,
                    assay.type = "counts", transposed = FALSE, ...){
                mat <- assay(x, assay.type)
                if(!transposed){
                    mat <- t(mat)
                }
                variable <- colData(x)[,variable]
                if(is.null(variable)){
                    stop("No data found in '",variable,
                        "' column of colData(x).", call. = FALSE)
                }
                performDMNgroupCV(x = mat, variable = variable, ...)
            }
)
