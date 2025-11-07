#' Generate a MaskedMatrix from Numeric Input
#'
#' Internal helper that ensures input is a 2D numeric matrix and returns a masked version,
#' replacing non-finite values (e.g., NA, NaN, Inf) with `NA` and recording their positions in a logical mask.
#'
#' @param mat A numeric matrix or vector. If a vector, it will be converted to a 1-row matrix.
#'
#' @return A list of class \code{"MaskedMatrix"} with two elements:
#' \describe{
#'   \item{data}{The original matrix with non-finite values replaced by \code{NA}.}
#'   \item{mask}{Logical matrix indicating non-finite entries (TRUE if missing).}
#' }
#'
#' @keywords internal

.mask_value_only <- function(mat) {
    #ensure matrix is at least 2D
    if (is.vector(mat)) {
        mat <- matrix(mat, nrow = 1)
    }
  
    #ensure matrix is not more than 2D
    if (length(dim(mat)) > 2) {
        stop("Input matrix can only have two dimensions or less")
    }
  
    #generate logical mask: TRUE where values are missing
    mask <- !is.finite(mat)  
  
    #create masked matrix
    masked.mat <- mat
    masked.mat[!is.finite(mat)] <- NA
  
    #return as a masked matrix
    return(structure(list(
        data = masked.mat,
        mask = mask
        ), class = "MaskedMatrix"))
}