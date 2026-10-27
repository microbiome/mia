
.calc_bulla_evenness <- function(mat, ...) {
    # Apply threshold, i.e., set values under and equal this threshold to 0
    mat <- .apply_threshold(mat, ...)

    # Species richness (number of species)
    S <- colSums2(mat > 0, na.rm = TRUE)

    # Relative abundances
    p <- t(mat)/colSums2(mat, na.rm = TRUE)

    i <- seq_len(nrow(p))
    O <- vapply(i,function(i){sum(pmin(p[i,], 1/S[i]))},numeric(1))

    # Bulla's Evenness
    (O - 1/S)/(1 - 1/S)
}

# Camargo's evenness. Inspired by code from Pepijn de Vries and Zhou Xiang
# at researchgate.net/post/How_can_we_calculate_the_Camargo_evenness_index_in_R
# but rewritten here. The index is defined over the S observed species, so
# zeroes are removed before computing it (microbiome/microbiome#231).
.calc_camargo_evenness <- function(mat, ...) {
    # Apply threshold, i.e., set values under and equal this threshold to 0
    mat <- .apply_threshold(mat, ...)
    res <- vapply(seq_len(ncol(mat)), function(i){
        # Observed species only
        x <- as.numeric(mat[, i])
        x <- x[!is.na(x) & x > 0]
        N <- length(x)
        if( N == 0 ){
            return(NA_real_)
        }
        # Sum of pairwise absolute differences, sum_{i<j} |x_i - x_j|,
        # computed from the sorted abundances
        xs <- sort(x)
        xx <- sum(xs * (2 * seq_len(N) - N - 1))
        1 - xx/(sum(x) * N)
    }, numeric(1))
    names(res) <- colnames(mat)
    res
}

# x: Species count vector
.calc_simpson_evenness <- function(mat, ...) {
    # Apply threshold, i.e., set values under and equal this threshold to 0
    mat <- .apply_threshold(mat, ...)

    # Species richness (number of detected species)
    S <- colSums2(mat > 0, na.rm = TRUE)

    # Simpson evenness (Simpson diversity per richness)
    .calc_inverse_simpson(mat)/S
}

# x: Species count vector
.calc_pielou_evenness <- function(mat, ...) {
    # Apply threshold, i.e., set values under and equal this threshold to 0
    mat <- .apply_threshold(mat, ...)

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
.calc_evar_evenness <- function(mat, ...) {
    # Apply threshold, i.e., set values under and equal this threshold to 0
    mat <- .apply_threshold(mat, ...)
    #
    N <- colSums2(mat, na.rm = TRUE)

    # Log abundance
    a <- log(mat)
    a[is.na(a) | is.infinite(a)] <- 0

    # Richness
    S <- colSums2(mat > 0, na.rm = TRUE)

    c <- colSums2(a, na.rm = TRUE)/S
    d <- t((t(a) - c)^2/S)
    d[mat == 0] <- 0

    f <- colSums2(d, na.rm = TRUE)

    (1 - 2/pi * atan(f))
}

# This function keeps only values over threshold and sets other to 0.
.apply_threshold <- function(mat, threshold = 0, ...){
    #
    if( !.is_a_numeric(threshold) ){
        stop("'threshold' must be a single numeric value.", call. = FALSE)
    }
    if( threshold > 0 ){
        mat[mat <= threshold] <- 0
    }
    return(mat)
}

.get_evenness_values <- function(index, mat, threshold = 0, ...){
    if(!is.numeric(threshold) || length(threshold) != 1L){
        stop("'threshold' must be a single numeric value.", call. = FALSE)
    }
    if(threshold > 0){
        mat[mat <= threshold] <- 0
    }
    FUN <- switch(index,
        camargo = .calc_camargo_evenness,
        pielou = .calc_pielou_evenness,
        simpson_evenness = .calc_simpson_evenness,
        evar = .calc_evar_evenness,
        bulla = .calc_bulla_evenness
        )
    res <- FUN(mat = mat, ...)
    res <- unname(res)
    return(res)
}

