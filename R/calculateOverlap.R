
.get_overlap <- function(x, detection = 0, ...){
    ############################# INPUT CHECK ##############################
    # Check detection
    if (!.is_numeric_string(detection)) {
        stop("'detection' must be a single numeric value or coercible to ",
            "one.",
            call. = FALSE)
    }
    detection <- as.numeric(detection)
    ########################### INPUT CHECK END ############################
    x <- t(x)
    # All the sample pairs
    sample_pairs <- as.matrix(expand.grid(colnames(x), colnames(x)))
        
    # Loop through all sample pairs
    result <- apply(sample_pairs, 1, FUN = function(sample_pair){
        # Get samples
        sample1 <- x[ , sample_pair[1]]
        sample2 <- x[ , sample_pair[2]]
        # Calculate overlap
        temp_result <- .calculate_overlap(sample1, sample2, detection)
        })
    # Create a matrix from result vector and give name to rownames and colnames
    result <- matrix(result, ncol = ncol(x))
    colnames(result) <- colnames(x)
    rownames(result) <- colnames(x)
        
    # Convert into distances
    result <- stats::as.dist(result)
    return(result)
}


################################ HELP FUNCTIONS ################################

.calculate_overlap <- function (x, y, detection) {
    # Take those taxa that have abundance over threshold
    inds <- which(x > detection & y > detection)
    x <- x[inds]
    y <- y[inds]
    # Overlap is the average of the sums of the values in each sample
    overlap <- (sum(x) + sum(y))/2
    return(overlap)
}
