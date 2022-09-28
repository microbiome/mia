
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

.is_an_integer <- function(x){
    is.numeric(x) && length(x) == 1L && x%%1==0
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

.check_altExp_present <- function(altexp, tse, 
                                  altExpName = .get_name_in_parent(altexp),
                                  tse_name = paste0("'", .get_name_in_parent(tse), "'") ){
    # Get class of object
    class <- as.character( class(tse) )
    # If the object does not have altExp slot
    if( !(is(tse, "TreeSummarizedExperiment") ||
          is(tse, "SingleCellExperiment")) ){
        stop("The class of ", tse_name, " is '", class, "' which does not have ",
             "an altExp slot. Please try '", altExpName, " = NULL'.", 
             call. = FALSE)
    }
    # If the object does not contain any altExps
    if( length(altExps(tse)) == 0 ){
        stop("altExps() of ", tse_name, " is empty. ",
             "Please try '", altExpName, " = NULL'.",
             call. = FALSE)
    }
    # altexp must specify altExp
    if( !( ( .is_an_integer(altexp) && altexp<length(altExps(tse)) && altexp>0) ||
           (.is_a_string(altexp) && altexp %in% altExpNames(tse)) ) ){
        stop("'", altExpName, "' must be integer or character specifying an ",
             "alternative experiment from ", tse_name, ".", call. = FALSE)
    }
}

.check_rowTree_present <- function(tree_name, x,
                                   name = .get_name_in_parent(tree_name) ){
    if( !.is_non_empty_string(tree_name) ){
        stop("'", name, "' must be a single non-empty character value.",
             call. = FALSE)
    }
    if( !(tree_name %in% names(x@rowTree)) ){
        stop("'", name, "' must specify a tree from 'x@rowTree'.",
             call. = FALSE)
    }
}

.check_colTree_present <- function(tree_name, x,
                                   name = .get_name_in_parent(tree_name) ){
    if( !.is_non_empty_string(tree_name) ){
        stop("'", name, "' must be a single non-empty character value.",
             call. = FALSE)
    }
    if( !(tree_name %in% names(x@colTree)) ){
        stop("'", name, "' must specify a tree from 'x@colTree'.",
             call. = FALSE)
    }
}

################################################################################
# internal wrappers for getter/setter

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


#' @importFrom S4Vectors metadata metadata<-
.add_values_to_metadata <- function(x, names, values){
    add_metadata <- as.list(values)
    names(add_metadata) <- names
    metadata(x) <- c(metadata(x), add_metadata)
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

#' Parse taxa in different taxonomic levels
#' @param taxa_tab `data.frame` object.
#' 
#' @param sep character string containing a regular expression, separator
#'  between different taxonomic levels, defaults to one compatible with both
#'  GreenGenes and SILVA `; |;"`.
#'  
#' @param column_name a single \code{character} value defining the column of taxa_tab
#'  that includes taxonomical information.
#'  
#' @param removeTaxaPrefixes {\code{TRUE} or \code{FALSE}: Should 
#'  taxonomic prefixes be removed? (default:
#'  \code{removeTaxaPrefixes = FALSE})}
#'  
#' @return  a `data.frame`.
#' @keywords internal
#' @importFrom IRanges CharacterList IntegerList
#' @importFrom S4Vectors DataFrame
#' @noRd
.parse_taxonomy <- function(taxa_tab, sep = "; |;", column_name = "Taxon",
                            removeTaxaPrefixes = FALSE, ...) {
    ############################### Input check ################################
    # Check sep
    if(!.is_non_empty_string(sep)){
      stop("'sep' must be a single character value.",
           call. = FALSE)
    }
    # Check column_name
    if( !(.is_non_empty_string(column_name) && column_name %in% colnames(taxa_tab)) ){
      stop("'column_name' must be a single character value defining column that includes",
           " information about taxonomic levels.",
           call. = FALSE)
    }
    # Check removeTaxaPrefixes
    if(!.is_a_bool(removeTaxaPrefixes)){
      stop("'removeTaxaPrefixes' must be TRUE or FALSE.", call. = FALSE)
    }
    ############################## Input check end #############################
    
    #  work with any combination of taxonomic ranks available
    all_ranks <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
    all_prefixes <- c("k__", "p__", "c__", "o__", "f__", "g__", "s__")
    
    # split the taxa strings
    taxa_split <- CharacterList(strsplit(taxa_tab[, column_name],sep))
    # extract present prefixes
    taxa_prefixes <- lapply(taxa_split, substr, 1L, 3L)
    # match them to the order given by present_prefixes
    taxa_prefixes_match <- lapply(taxa_prefixes, match, x = all_prefixes)
    taxa_prefixes_match <- IntegerList(taxa_prefixes_match)
    # get the taxa values
    if(removeTaxaPrefixes){
      taxa_split <- lapply(taxa_split,
                           gsub,
                           pattern = "([kpcofgs]+)__",
                           replacement = "")
      taxa_split <- CharacterList(taxa_split)
    }
    # extract by order matches
    taxa_split <- taxa_split[taxa_prefixes_match]
    #
    if(length(unique(lengths(taxa_split))) != 1L){
      stop("Internal error. Something went wrong while splitting taxonomic levels.",
           "Please check that 'sep' is correct.", call. = FALSE)
    }
    taxa_tab <- DataFrame(as.matrix(taxa_split))
    colnames(taxa_tab) <- all_ranks
    
    taxa_tab
}
