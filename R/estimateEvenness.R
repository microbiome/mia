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
#' @param BPPARAM A
#'   \code{\link[BiocParallel:BiocParallelParam-class]{BiocParallelParam}}
#'   object specifying whether calculation of estimates should be parallelized.
#'   (Currently not used)
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
                                       x = x,
                                       mat = assay(x, abund_values),
                                       BPPARAM = BPPARAM, ...)
        .add_values_to_colData(x, vnss, name)
    }
)

.calc_bulla_evenness <- function(x, zeroes = TRUE) {
    if(!.is_a_bool(zeroes)){
        stop("'zeroes' must be TRUE or FALSE.", call. = FALSE)
    }

    if (!zeroes) {
        x[x > 0]
    }

    # Species richness (number of species)
    S <- sum(x > 0, na.rm = TRUE)

    # Relative abundances
    p <- x/sum(x)

    O <- sum(pmin(p, 1/S))

    # Bulla's Evenness
    (O - 1/S)/(1 - 1/S)
}

# Camargo's eveness x: species counts zeroes: include zeros Inspired
# by code from Pepijn de Vries and Zhou Xiang at
# researchgate.net/post/How_can_we_calculate_the_Camargo_evenness_index_in_R
# but rewritten here
.calc_camargo_evenness <- function(x, zeroes = TRUE) {
    if(!.is_a_bool(zeroes)){
        stop("'zeroes' must be TRUE or FALSE.", call. = FALSE)
    }

    if (!zeroes) {
        x[x > 0]
    }

    N <- sum(x > 0, na.rm = TRUE)

    xx <- 0
    for (i in seq_len(N - 1)) {
        xx <- xx + sum(abs(x[(i + 1):N] - x[i]))
    }

    # Return
    1 - xx/(sum(x) * N)
}

# x: Species count vector
.calc_simpson_evenness <- function(x) {

    # Species richness (number of species)
    S <- sum(x > 0, na.rm = TRUE)

    # Simpson index
    lambda <- .get_simpson(x)

    # Simpson evenness (Simpson diversity per richness)
    (1/lambda)/S
}

# x: Species count vector
.calc_pielou_evenness <- function(x) {

    # Remove zeroes
    x <- x[x > 0]

    # Species richness (number of detected species)
    S <- sum(x > 0, na.rm = TRUE)

    # Relative abundances
    p <- x/sum(x)

    # Shannon index
    H <- (-sum(p * log(p)))

    # Simpson evenness
    H/log(S)
}

# Smith and Wilson’s Evar index
.calc_evar_evenness <- function(x, zeroes = TRUE) {
    if(!.is_a_bool(zeroes)){
        stop("'zeroes' must be TRUE or FALSE.", call. = FALSE)
    }

    if (!zeroes) {
        x[x > 0]
    }

    n <- sum(x, na.rm = TRUE)
    d <- rep(NA, n)

    # Log abundance
    a <- log(x)
    a[is.na(a) | is.infinite(a)] <- 0

    # Richness
    S <- sum(x > 0)

    b <- a/S
    c <- sum(b)
    d <- (a - c)^2/S
    d[x != 0] <- 0

    f <- sum(d)

    (1 - 2/pi * atan(f))
}

.run_evenness <- function(x, i, mat, detection = 0, ...){
    if(!is.numeric(detection) || length(detection) == 1L){
        stop("'detection' must be a single numeric value.", call. = FALSE)
    }
    if(detection > 0){
        x[x <= detection] <- 0
    }
    vnss_FUN <- switch(i,
                       camargo = .calc_camargo_evenness,
                       pielou = .calc_pielou_evenness,
                       simpson = .calc_simpson_evenness,
                       evar = .calc_evar_evenness,
                       bulla = .calc_bulla_evenness)
    vnss <- vnss_FUN(mat, ...)
    vnss
}
