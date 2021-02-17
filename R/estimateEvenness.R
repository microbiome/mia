#' Estimate Evenness
#'
#' Several functions for calculation of evenness estimates are available.
#'
#' These includes the \sQuote{Camargo}, \sQuote{Pielou}, \sQuote{Simpson},
#' \sQuote{Evar} and \sQuote{Bulla} evenness measures.
#'
#' @param x a \code{\link{SummarizedExperiment}} object
#'
#' @param abund_values the name of the assay used for calculation of the
#'   sample-wise estimates
#'
#' @param index a evenness measurement
#'
#' @param name a name for the column of the colData the results should be stored
#'   in.
#'
#' @param threshold a numeric threshold. assay values below or equal to this
#'   threshold will be set to zero.
#'
#' @param BPPARAM A
#'   \code{\link[BiocParallel:BiocParallelParam-class]{BiocParallelParam}}
#'   object specifying whether calculation of estimates should be parallelized.
#'
#' @param ... additional parameters passed to \code{estimateEvenness}
#'
#' @return \code{x} with additional \code{\link{colData}} named
#'   \code{*name*}
#'
#' @seealso
#' \code{\link[scater:plotColData]{plotColData}}
#'
#' @name estimateEvenness
#'
#' @examples
#' data(GlobalPatterns)
#' se <- GlobalPatterns
#'
#' se <- estimateEvenness(se)
#' head(colData(se))
NULL

#' @rdname estimateEvenness
#' @export
setGeneric("estimateEvenness",signature = c("x"),
           function(x, abund_values = "counts",
                    index = c("camargo", "pielou", "simpson", "evar", "bulla"),
                    name = index, ...)
               standardGeneric("estimateEvenness"))

#' @rdname estimateEvenness
#' @export
setMethod("estimateEvenness", signature = c(x = "SummarizedExperiment"),
    function(x, abund_values = "counts",
             index = c("camargo", "pielou", "simpson", "evar", "bulla"),
             name = index, ..., BPPARAM = SerialParam()){
        # input check
        index<- match.arg(index, several.ok = TRUE)
        if(!.is_non_empty_character(name) || length(name) != length(index)){
            stop("'name' must be a non-empty character value and have the ",
                 "same length than 'index'.",
                 call. = FALSE)
        }
        .check_abund_values(abund_values, x)
        #
        vnss <- BiocParallel::bplapply(index,
                                       .run_evenness,
                                       mat = assay(x, abund_values),
                                       BPPARAM = BPPARAM, ...)
        .add_values_to_colData(x, vnss, name)
    }
)

.calc_bulla_evenness <- function(mat) {
    # Species richness (number of species)
    S <- colSums2(mat > 0, na.rm = TRUE)

    # Relative abundances
    p <- t(mat)/colSums2(mat, na.rm = TRUE)

    i <- seq_len(nrow(p))
    O <- vapply(i,function(i){sum(pmin(p[i,], 1/S[i]))},numeric(1))

    # Bulla's Evenness
    (O - 1/S)/(1 - 1/S)
}

# Camargo's eveness x: species counts zeroes: include zeros Inspired
# by code from Pepijn de Vries and Zhou Xiang at
# researchgate.net/post/How_can_we_calculate_the_Camargo_evenness_index_in_R
# but rewritten here
.calc_camargo_evenness <- function(mat) {
    N <- colSums2(mat > 0, na.rm = TRUE)

    seq <- IntegerList(lapply(N - 1,seq_len))

    x <- mapply(
        function(i, n, s){
            xx <- 0
            for (j in s) {
                xx <- xx + sum(abs(mat[(j + 1):n,i] - mat[j,i]))
            }
            xx
        },
        seq_along(N),
        N,
        seq)
    # Return
    1 - x/(colSums2(mat, na.rm = TRUE) * N)
}

# x: Species count vector
.calc_simpson_evenness <- function(mat) {
    # Species richness (number of species)
    S <- colSums2(mat > 0, na.rm = TRUE)

    # Simpson index
    lambda <-  1 - .get_simpson(mat)

    # Simpson evenness (Simpson diversity per richness)
    (1/lambda)/S
}

# x: Species count vector
.calc_pielou_evenness <- function(mat) {
    # Remove zeroes
    mat[mat == 0] <- NA

    # Species richness (number of detected species)
    S <- colSums2(mat > 0, na.rm = TRUE)

    # Relative abundances
    p <- t(mat)/colSums2(mat, na.rm = TRUE)

    # Shannon index
    H <- (-rowSums2(p * log(p), na.rm = TRUE))

    # Simpson evenness
    H/log(S)
}

# Smith and Wilson’s Evar index
.calc_evar_evenness <- function(mat) {
    N <- colSums2(mat, na.rm = TRUE)

    # Log abundance
    a <- log(mat)
    a[is.na(a) | is.infinite(a)] <- 0

    # Richness
    S <- colSums2(mat > 0, na.rm = TRUE)

    b <- t(a)/S
    c <- rowSums2(b, na.rm = TRUE)
    d <- t((t(a) - c)^2/S)
    d[mat == 0] <- 0

    f <- colSums2(d, na.rm = TRUE)

    (1 - 2/pi * atan(f))
}

.run_evenness <- function(index, mat, threshold = 0, ...){
    if(!is.numeric(threshold) || length(threshold) != 1L){
        stop("'threshold' must be a single numeric value.", call. = FALSE)
    }
    if(threshold > 0){
        mat[mat <= threshold] <- 0
    }
    vnss_FUN <- switch(index,
                       camargo = .calc_camargo_evenness,
                       pielou = .calc_pielou_evenness,
                       simpson = .calc_simpson_evenness,
                       evar = .calc_evar_evenness,
                       bulla = .calc_bulla_evenness)
    vnss_FUN(mat)
}
