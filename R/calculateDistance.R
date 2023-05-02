# calculateDistance function is removed. This function is left for other functions 
# to use
.calculate_distance <- function(mat, FUN = stats::dist, ...){
    # Distance between all samples against all samples
    do.call(FUN, c(list(mat),list(...)))
}
