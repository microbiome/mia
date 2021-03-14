
################################################################################
# internal methods loaded from other packages

.get_mat_from_sce <- scater:::.get_mat_from_sce
.get_mat_for_reddim <- scater:::.get_mat_for_reddim

################################################################################
# integration with other packages

.require_package <- function(pkg){
    if(!requireNamespace(pkg, quietly = TRUE)){
    stop("'",pkg,"' package not found. Please install the '",pkg,"' package ",
         "to use this function.", call. = FALSE)
    }
}

################################################################################
# testing

.is_a_bool <- function(x){
  is.logical(x) && length(x) == 1L && !is.na(x)
}

.is_non_empty_character <- function(x){
  is.character(x) && all(nzchar(x))
}

.is_non_empty_string <- function(x){
  .is_non_empty_character(x) && length(x) == 1L
}

.is_a_string <- function(x){
  is.character(x) && length(x) == 1L
}

.are_whole_numbers <- function(x){
  tol <- 100 * .Machine$double.eps
  abs(x - round(x)) <= tol && !is.infinite(x)
}

.is_numeric_string <- function(x){
  x <- as.character(x)
  suppressWarnings({x <- as.numeric(x)})
  !is.na(x)
}

.is_function <- function(x){
  typeof(x) == "closure" && is(x, "function")
}

.all_are_existing_files <- function(x){
  all(file.exists(x))
}

.get_name_in_parent <- function(x) {
  .safe_deparse(do.call(substitute, list(substitute(x), parent.frame())))
}

.safe_deparse <- function (expr, ...) {
  paste0(deparse(expr, width.cutoff = 500L, ...), collapse = "")
}

################################################################################
# checks

#' @importFrom SummarizedExperiment assays
.check_abund_values <- function(abund_values, x,
                                name = .get_name_in_parent(abund_values)){
    if(!.is_non_empty_string(abund_values)){
        stop("'",name,"' must be a single non-empty character value.",
             call. = FALSE)
    }
    if(!(abund_values %in% names(assays(x)))){
        stop("'",name,"' must be a valid name of assays(x)", call. = FALSE)
    }
}

################################################################################
# internal wrappers for getter/setter

#' @importFrom SummarizedExperiment colData colData<-
#' @importFrom S4Vectors DataFrame
.add_values_to_colData <- function(x, values, name){
    values <- mapply(
        function(value, n){
            value <- DataFrame(value)
            colnames(value)[1L] <- n
            if(ncol(value) > 1L){
                i <- seq.int(2,ncol(value))
                colnames(value)[i] <- paste0(n,"_",colnames(value)[i])
            }
            value
        },
        values,
        name)

    
    newdf <- DataFrame(values)

    # Check if the new column data contains variables with
    # identical names compared to the original column data
    if any(colnames(newdf) %in% colData(x)) {
        shared <- intersect(colnames(newdf), colData(x))
        warning(paste0("The original and newly added colData(...) contains the following 
            identical field names: ", paste(shared, collapse = ", "), ". The new colData(...) 
            will thus contain duplicated field names. Consider using the 'name' argument in 
            the function call to specify alternative names."))
    }

    colData(x) <- cbind(colData(x),newdf)
    
    x
}
