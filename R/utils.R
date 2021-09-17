
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
.check_assay_present <- function(assay_name, x,
                                 name = .get_name_in_parent(assay_name)){
    if(!.is_non_empty_string(assay_name)){
        stop("'",name,"' must be a single non-empty character value.",
             call. = FALSE)
    }
    if(!(assay_name %in% names(assays(x)))){
        stop("'",name,"' must be a valid name of assays(x)", call. = FALSE)
    }
}

################################################################################
# internal wrappers for getter/setter
#' @importFrom S4Vectors isSingleNumber
#' @importFrom S4Vectors new2
.make_zero_col_DataFrame <- function (nrow = 0L){
    # Sourced by S4vectors
    stopifnot(isSingleNumber(nrow))
    if (!is.integer(nrow)) 
        nrow <- as.integer(nrow)
    stopifnot(nrow >= 0L)
    new2("DFrame", nrows = nrow, check = FALSE)
}




#' @importFrom SummarizedExperiment colData colData<-
#' @importFrom S4Vectors DataFrame
.add_values_to_colData <- function(x, values, name){
    # converts each value:name pair into a DataFrame
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
    values <- do.call(cbind, values)
    # check for duplicated values
    f <- colnames(colData(x)) %in% colnames(values)
    if(any(f)) {
        warning("The following values are already present in `colData` and ",
                "will be overwritten: '",
                paste(colnames(colData(x))[f], collapse = "', '"),
                "'. Consider using the 'name' argument to specify alternative ",
                "names.",
                call. = FALSE)
    }
    # keep only unique values
    colData(x) <- cbind(colData(x)[!f], values)
    x
}

# keep dimnames of feature table (assay) consistent with the meta data 
# of sample (colData) and feature (rowData)
.set_feature_tab_dimnames <- function(feature_tab, 
                                      sample_meta, 
                                      feature_meta) {
    if (nrow(sample_meta) > 0 || ncol(sample_meta) > 0) {
        if (ncol(feature_tab) != nrow(sample_meta) 
            || !setequal(colnames(feature_tab), rownames(sample_meta))) {
            stop(
                "The sample ids in feature table are not incompatible ",
                "with those in sample meta",
                call. = FALSE
            )
        }
        if (!identical(colnames(feature_tab), rownames(sample_meta))) {
            feature_tab <- feature_tab[, rownames(sample_meta), drop = FALSE]
        }
    }
    
    if (nrow(feature_meta) > 0 || ncol(feature_meta) > 0) {
        if (nrow(feature_tab) != nrow(feature_meta)
            || !setequal(rownames(feature_tab), rownames(feature_meta))) {
            stop(
                "The feature names in feature table are not incompatible ",
                "with those in feature meta",
                call. = FALSE
            )
        }
        if (!identical(rownames(feature_tab), rownames(feature_meta))) {
            feature_tab <- feature_tab[rownames(feature_meta), , drop = FALSE]
        }
    }
  
    feature_tab
}
