################################################################################
# This function gives user a message when the package is loaded into the session
.onAttach <- function(libname, pkgname) {
    pkg_version <- utils::packageDescription(pkgname, fields = "Version")
    msg <- paste0(
        "This is ", pkgname, " version ", pkg_version, "\n",
        "- Online documentation and vignettes: https://microbiome.github.io/",
        pkgname, "/",
        "\n",
        "- Online book 'Orchestrating Microbiome Analysis (OMA)': ",
        "https://microbiome.github.io/OMA/docs/devel/"
    )
    packageStartupMessage(msg)
}

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

.is_integer <- function(x){
    is.numeric(x) && all(x%%1==0)
}

.is_an_integer <- function(x){
    .is_integer(x) && length(x) == 1L
}

.are_whole_numbers <- function(x){
    tol <- 100 * .Machine$double.eps
    abs(x - round(x)) <= tol && !is.infinite(x)
}

.is_a_numeric <- function(x){
    is.numeric(x) && length(x) == 1L
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
.check_assay_present <- function(
        assay.type, x, name = .get_name_in_parent(assay.type)){
    if(!.is_non_empty_string(assay.type)){
        stop("'",name,"' must be a single non-empty character value.",
             call. = FALSE)
    }
    if(!(assay.type %in% names(assays(x)))){
        stop("'",name,"' must be a valid name of assays(x)", call. = FALSE)
    }
}

.check_rowTree_present <- function(
        tree.name, x, name = .get_name_in_parent(tree.name) ){
    if( !.is_non_empty_string(tree.name) ){
        stop("'", name, "' must be a single non-empty character value.",
             call. = FALSE)
    }
    if( !(tree.name %in% rowTreeNames(x)) ){
        stop("'", name, "' must specify a tree from 'rowTreeNames(x)'.",
             call. = FALSE)
    }
}

.check_colTree_present <- function(
        tree.name, x, name = .get_name_in_parent(tree.name) ){
    if( !.is_non_empty_string(tree.name) ){
        stop("'", name, "' must be a single non-empty character value.",
             call. = FALSE)
    }
    if( !(tree.name %in% colTreeNames(x)) ){
        stop("'", name, "' must specify a tree from 'colTreeNames(x)'.",
             call. = FALSE)
    }
}

# Check if alternative experiment can be found from altExp slot.
.check_altExp_present <- function(
        altexp, tse, altExpName = .get_name_in_parent(altexp),
        tse_name = .get_name_in_parent(tse), .disable.altexp = FALSE, ...){
    # Disable altExp if specified
    if( !.is_a_bool(.disable.altexp) ){
        stop("'.disable.altexp' must be TRUE or FALSE.", call. = FALSE)
    }
    if( .disable.altexp ){
        altexp <- NULL
    }
    # Check that altexp.name must be an integer or name
    if( !(.is_a_string(altexp) || .is_an_integer(altexp) || is.null(altexp)) ){
        stop(
            "'", altExpName, "' must be a string or an integer.", call. = FALSE)
    }
    # If is not NULL, but the object does not have altExp slot
    if( !is.null(altexp) && !is(tse, "SingleCellExperiment") ){
        stop(
            "'", altExpName, "', is specified but '", tse_name, "' does not ",
            "have altExp slot.", call. = FALSE)
    }
    # Then check that altExp can be found; name or index.
    if( !is.null(altexp) && !altexp %in% c(
        altExpNames(tse), seq_len(length(altExps(tse)))) ){
        stop(
            "'", altExpName, "', does not specify an experiment from altExp ",
            "slot of '", tse_name, "'.", call. = FALSE)
    }
}

# Check whether dimred is present in tse
.check_dimred_present <- function(dimred, x){
    specifies_index <- .is_integer(dimred) && dimred > 0 &&
        dimred <= length(reducedDims(x))
    specifies_name <- .is_a_string(dimred) && dimred %in% reducedDimNames(x)
    if( !specifies_index && !specifies_name ){
        stop("'dimred' must specify name or index from reducedDims(x).",
             call. = FALSE)
    }
    return(NULL)
}

# Check MARGIN parameters. Should be defining rows or columns.
.check_MARGIN <- function(MARGIN, name = .get_name_in_parent(MARGIN)) {
    # MARGIN must be one of the following options
    if( !(length(MARGIN) == 1L && tolower(MARGIN) %in% c(
        1, 2, "1", "2", "features", "samples", "columns", "col", "row",
        "rows", "cols")) ) {
        stop("'", name,"' must be 'rows' or 'cols'.", call. = FALSE)
    }
    # Convert MARGIN to numeric if it is not.
    MARGIN <- ifelse(tolower(MARGIN) %in% c(
        "samples", "columns", "col", 2, "cols"), 2, 1)
    return(MARGIN)
}

# Check if metadata has the specified data.
.check_metadata_present <- function(
        data.type, x, name = .get_name_in_parent(data.type)){
    if( !.is_non_empty_string(data.type) ){
        stop("'" ,name, "' must be a single non-empty character value.",
             call. = FALSE)
    }
    if( !(data.type %in% names(metadata(x))) ){
        stop("'",name,"' must be a valid name of metadata(x)", call. = FALSE)
    }
    return(data.type)
}
################################################################################
# Internal wrappers for getters

# Input: (Tree)SE
# Output: (Tree)SE
.check_and_get_altExp <- function(
        x, altexp = NULL, ...){
    # If altexp is specified, check and get it.
    # Otherwise return the original object
    if( !is.null(altexp) ){
        # Check altexp
        .check_altExp_present(altexp, x, ...)
        # Get altExp and return it
        x <- altExp(x, altexp)
    }
    return(x)
}

################################################################################
# Internal wrappers for setters

# This function adds values to colData (or rowData). The data must be in a list.
# Each element of list represent a column to be added to col/rowData.
#' @importFrom SummarizedExperiment colData colData<- rowData rowData<-
#' @importFrom S4Vectors DataFrame
.add_values_to_colData <- function(
        x, values, name, altexp = NULL, MARGIN = default.MARGIN,
        default.MARGIN = 2, transpose.MARGIN = FALSE, colname = "name",
        ...){
    #
    if( !.is_a_string(colname) ){
        stop("'colname' must be a string.", call. = FALSE)
    }
    #
    # Check if altExp can be found
    .check_altExp_present(altexp, x)
    # Check that MARGIN is correct
    MARGIN <- .check_MARGIN(MARGIN)
    #
    # If transpose.MARGIN is TRUE, transpose MARGIN, i.e. 1 --> 2, and 2 --> 1.
    # In certain functions, values calculated by rows (MARGIN=1) are stored to
    # colData (MARGIN=2) and vice versa.
    if( transpose.MARGIN ){
        MARGIN <- ifelse(MARGIN == 1, 2, 1)
    }
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
    
    # Based on MARGIN, get rowDatra or colData
    FUN <- switch(MARGIN, rowData, colData)
    # If altexp.name was not NULL, then we know that it specifies correctly
    # altExp from the slot. Take the colData/rowData from experiment..
    if( !is.null(altexp) ){
        cd <- FUN( altExp(x, altexp) )
    } else{
        cd <- FUN(x)
    }
    
    # check for duplicated values
    f <- colnames(cd) %in% colnames(values)
    FUN_name <- switch(MARGIN, "rowData", "colData")
    if(any(f)) {
        warning(
            "The following values are already present in `", FUN_name,
            "` and will be overwritten: '",
            paste(colnames(cd)[f], collapse = "', '"),
            "'. Consider using the '", colname,
            "' argument to specify alternative names.",
            call. = FALSE)
    }
    # Keep only unique values
    cd <- cbind( (cd)[!f], values )
    
    # Replace colData with new one
    x <- .add_to_coldata(x, cd, altexp = altexp, MARGIN = MARGIN)
    return(x)
}

# Get feature or sample metadata. Allow hidden usage of MARGIN and altExp.
#' @importFrom SummarizedExperiment rowData colData
.add_to_coldata <- function(
        x, cd, altexp = NULL, .disable.altexp = FALSE,
        MARGIN = default.MARGIN, default.MARGIN = 1, ...){
    #
    if( !.is_a_bool(.disable.altexp) ){
        stop("'.disable.altexp' must be TRUE or FALSE.", call. = FALSE)
    }
    # Check if altExp can be found
    .check_altExp_present(altexp, x, ...)
    # Check that MARGIN is correct
    MARGIN <- .check_MARGIN(MARGIN)
    # Based on MARGIN, add result to rowData or colData
    FUN <- switch(MARGIN, `rowData<-`, `colData<-`)
    # If altexp was specified, add result to altExp. Otherwise add it directly
    # to x.
    if( !is.null(altexp) && !.disable.altexp ){
        altExp(x, altexp) <- FUN( altExp(x, altexp), value = cd )
    } else{
        x <- FUN(x, value = cd)
    }
    return(x)
}

#' @importFrom S4Vectors metadata metadata<-
.add_values_to_metadata <- function(
        x, names, values, altexp = NULL, metadata.name = "name", ...){
    #
    if( !.is_a_string(metadata.name) ){
        stop("'metadata.name' must be a string.", call. = FALSE)
    }
    # Check if altExp can be found
    .check_altExp_present(altexp, x)
    #
    # Create a list and name elements
    add_metadata <- list(values)
    names(add_metadata) <- names
    # Get old metadata
    if( !is.null(altexp) ){
        old_metadata <- metadata( altExp(x, altexp) )
    } else{
        old_metadata <- metadata(x)
    }
    # Check if names match with elements that are already present
    f <- names(old_metadata) %in% names(add_metadata)
    if( any(f) ){
        warning(
            "The following values are already present in `metadata` and will ",
            "be overwritten: '",
            paste(names(old_metadata)[f], collapse = "', '"),
            "'. Consider using the '", metadata.name,
            "' argument to specify alternative ", "names.", call. = FALSE)
    }
    # keep only unique values
    add_metadata <- c( old_metadata[!f], add_metadata )
    # Add metadata to altExp or directly to x
    if( !is.null(altexp) ){
        metadata( altExp(x, altexp) ) <- add_metadata
    } else{
        metadata(x) <- add_metadata
    }
    return(x)
}

# This function can be used to add values to altExp
.add_to_altExps <- function(x, values, name = names(values), ...){
    # Check values
    if( !((is(values, "list") || is(values, "SimpleList")) &&
          length(values) > 0) ){
        stop("'values' must be non-empty list.", call. = FALSE)
    }
    # Check names
    if( !is.character(name) && length(name) > 1L ){
        stop("'name' must be a character value.", call. = FALSE)
    }
    # Names must match with list
    if( length(values) != length(name) ){
        stop("Lenght of 'name' must match with 'values'.", call. = FALSE)
    }
    #
    # If the object is SE, convert it to TreeSE
    if( !is(x, "SingleCellExperiment") ){
        x <- as(x, "TreeSummarizedExperiment")
        warning(
            "SummarizedExperiment does not have altExps slot. ",
            "Therefore, it is converted to TreeSummarizedExperiment.",
            call. = FALSE)
    }
    #
    # Add names to values
    names(values) <- name
    # Get altExps
    old_altexp <- altExps(x)
    # Check if names match with elements that are already present
    f <- names(old_altexp) %in% names(values)
    if( any(f) ){
        warning(
            "The following values are already present in `altExps` and will ",
            "be overwritten: '",
            paste(names(old_altexp)[f], collapse = "', '"),
            "'. Consider using the 'name' argument to specify alternative ",
            "names.", call. = FALSE)
    }
    # Keep only unique values
    values <- c( old_altexp[!f], values )
    # Add to altExps
    altExps(x) <- values
    return(x)
}

# This function can be used to add values to reducedDims
.add_values_to_reducedDims <- function(x, values, name, ...){
    # Check values
    if( !((is(values, "matrix") || is(values, "dist")) && length(values) > 0) ){
        stop("'values' must be a matrix.", call. = FALSE)
    }
    # Check names
    if( !.is_a_string(name) ){
        stop("'name' must be a character value.", call. = FALSE)
    }
    #
    # If the object is SE, convert it to TreeSE
    if( !is(x, "SingleCellExperiment") ){
        x <- as(x, "TreeSummarizedExperiment")
        warning(
            "SummarizedExperiment does not have reducedDims slot. ",
            "Therefore, it is converted to TreeSummarizedExperiment.",
            call. = FALSE)
    }
    if( !identical(rownames(as.matrix(values)), colnames(x)) ){
        stop("Rownames of the matrix should match with colnames(x).",
             " The result is not added to reducedDims.", call. = FALSE)
    }
    # Throw warning if values of reducedDim are overwritten
    if( name %in% names(reducedDims(x)) ){
        warning(
            "The following values are already present in `reducedDims` and",
            " will be overwritten: '", name,
            "'. Consider using the 'name' argument to specify alternative ",
            "names.", call. = FALSE)
    }
    reducedDim(x, name) <- values
    return(x)
}

################################################################################
# Other common functions

# keep dimnames of feature table (assay) consistent with the meta data
# of sample (colData) and feature (rowData)
.set_feature_tab_dimnames <- function(
        feature_tab, sample_meta, feature_meta, refseq = NULL,
        include.all = FALSE, ...) {
    #
    if( !.is_a_bool(include.all) ){
        stop("'include.all' must be TRUE or FALSE.", call. = FALSE)
    }
    # Sample and feature names must be present
    if( is.null(colnames(feature_tab)) || is.null(rownames(sample_meta)) ){
        stop("Sample ids must be present.", call. = FALSE)
    }
    if( is.null(rownames(feature_tab)) || is.null(rownames(feature_meta)) ){
        stop("Feature ids must be present.", call. = FALSE)
    }
    #
    # Metadata can include features that are not present in abundance table.
    # With 'include.all', it is possible to include them in the dataset.
    if( include.all ){
        samples <- union(colnames(feature_tab), rownames(sample_meta))
        features <- union(rownames(feature_tab), rownames(feature_meta))
        # Add additional columns and sort data
        feature_tab <- feature_tab[
            match(features, rownames(feature_tab)),
            match(samples, colnames(feature_tab)), drop = FALSE]
        sample_meta <- sample_meta[
            match(samples, rownames(sample_meta)), , drop = FALSE]
        feature_meta <- feature_meta[
            match(features, rownames(feature_meta)), , drop = FALSE]
        # Add correct names
        colnames(feature_tab) <- rownames(sample_meta) <- samples
        rownames(feature_tab) <- rownames(feature_meta) <- features
    }
    
    # The abundance table takes the precedence, and metadata is modified
    # accordingly. Give warning if there are samples or features that are not
    # included in the metadata.
    not_found <- sum(!colnames(feature_tab) %in% rownames(sample_meta))
    if( not_found != 0 ){
        warning("The dataset includes ", not_found, " samples that do not ",
                "have metadata. Please check for errors.", call. = FALSE)
    }
    not_found <- sum(!rownames(feature_tab) %in% rownames(feature_meta))
    if( not_found != 0 ){
        warning("The dataset includes ", not_found, " features that are not ",
                "included in the taxonomy table. Please check for errors.",
                call. = FALSE)
    }
    # If the metadata includes samples or features that will be removed give
    # warning for the user.
    not_found <- sum(!rownames(sample_meta) %in% colnames(feature_tab))
    if( not_found != 0 ){
        warning("The sample metadata includes ", not_found, " samples that ",
                "are not included in the abundance table and thus being ",
                "removed. Please check for errors.", call. = FALSE)
    }
    not_found <- sum(!rownames(feature_meta) %in% rownames(feature_tab))
    if( not_found != 0 ){
        warning("The taxonomy table includes ", not_found, " features that ",
                "are not present in the abundance table and thus being ",
                "removed. Please check for errors.", call. = FALSE)
    }
    
    # We order the metadata based on abundance table. Moreover, we subset
    # the metadata to match with abundance table if there are additional data.
    ind <- match(colnames(feature_tab), rownames(sample_meta))
    sample_meta <- sample_meta[ind, , drop = FALSE]
    rownames(sample_meta) <- colnames(feature_tab)
    ind <- match(rownames(feature_tab), rownames(feature_meta))
    feature_meta <- feature_meta[ind, , drop = FALSE]
    rownames(feature_meta) <- rownames(feature_tab)
    # Reference sequences are optional and all the features must have sequences.
    # This is because DNAStringSet object cannot have empty element for
    # features that were not included.
    if( !is.null(refseq) && all(rownames(feature_tab) %in% names(refseq)) ){
        refseq <- refseq[ rownames(feature_tab) ]
    } else if( !is.null(refseq) ){
        warning("Reference sequences are incompatible with the data.",
                call. = FALSE)
        refseq <- NULL
    } else{
        refseq <- NULL
    }
    # Return a list of these tables
    data_list <- list(
        assay = feature_tab,
        rowData = feature_meta,
        colData = sample_meta,
        referenceSeq = refseq
    )
    return(data_list)
}

#' Parse taxa in different taxonomic levels
#' @param taxa_tab `data.frame` object.
#'
#' @param sep character string containing a regular expression, separator
#'  between different taxonomic levels, defaults to one compatible with both
#'  GreenGenes and SILVA `; |;"`.
#'
#' @param col.name a single \code{character} value defining the column of
#' taxa_tab that includes taxonomical information.
#'
#' @param prefix.rm {\code{TRUE} or \code{FALSE}: Should
#'  taxonomic prefixes be removed? (default: \code{prefix.rm = FALSE})}
#'
#' @return  a `data.frame`.
#' @keywords internal
#' @importFrom IRanges CharacterList IntegerList
#' @importFrom S4Vectors DataFrame
#' @noRd
.parse_taxonomy <- function(
        taxa_tab, sep = "; |;", col.name = column_name, column_name = "Taxon",
        remove.prefix = prefix.rm, prefix.rm = removeTaxaPrefixes,
        removeTaxaPrefixes = FALSE, ...) {
    ############################### Input check ################################
    # Check sep
    if(!.is_non_empty_string(sep)){
        stop("'sep' must be a single character value.", call. = FALSE)
    }
    # Check col.name
    if( !(.is_non_empty_string(col.name) && col.name %in% colnames(taxa_tab)) ){
        stop("'col.name' must be a single character value defining column ",
             "that includes information about taxonomic levels.", call. = FALSE)
    }
    # Check remove.prefix
    if(!.is_a_bool(remove.prefix)){
        stop("'remove.prefix' must be TRUE or FALSE.", call. = FALSE)
    }
    ############################## Input check end #############################
    
    #  work with any combination of taxonomic ranks available
    all_ranks <- .taxonomy_rank_prefixes
    all_prefixes <- paste0(all_ranks, "__")
    names(all_prefixes) <- names(all_ranks)
    
    # split the taxa strings
    taxa_split <- CharacterList(strsplit(taxa_tab[, col.name],sep))
    # extract present prefixes
    taxa_prefixes <- lapply(
        taxa_split, gsub, pattern = "(^[a-zA-Z]+__).*", replacement = "\\1")
    # match them to the order given by present_prefixes
    taxa_prefixes_match <- lapply(taxa_prefixes, match, x = all_prefixes)
    taxa_prefixes_match <- IntegerList(taxa_prefixes_match)
    # get the taxa values without prefixes
    if(remove.prefix){
        pattern <- paste0("(", paste0(all_ranks, collapse = "|"), ")__")
        taxa_split <- lapply(
            taxa_split, gsub, pattern = pattern, replacement = "")
        taxa_split <- CharacterList(taxa_split)
    }
    # extract by order matches
    taxa_split <- taxa_split[taxa_prefixes_match]
    #
    if(length(unique(lengths(taxa_split))) != 1L){
        stop("Something went wrong while splitting taxonomic ",
             "levels. Please check that 'sep' is correct.", call. = FALSE)
    }
    taxa_tab <- DataFrame(as.matrix(taxa_split))
    colnames(taxa_tab) <- names(all_ranks)
    
    # Subset columns so that include only those columns that have some
    # information
    non_empty <- colSums(is.na(taxa_tab)) != nrow(taxa_tab)
    taxa_tab <- taxa_tab[ , non_empty, drop = FALSE]
    
    return(taxa_tab)
}

################################################################################
# internal wrappers for agglomerateByRank/agglomerateByVariable
.merge_features <- function(x, merge.by = rank, rank = NULL, ...) {
    # Check if merge.by parameter belongs to taxonomyRanks
    if( .is_a_string(merge.by) && merge.by %in% taxonomyRanks(x) ){
        # Merge using agglomerateByRank
        x <- agglomerateByRank(x, rank = merge.by, ...)
    } else if( !is.null(merge.by) ){
        # Merge using agglomerateByVariable
        x <- agglomerateByVariable(x, by = "rows", group = merge.by, ...)
    }
    return(x)
}

################################################################################
# This function sets taxonomy ranks based on rowData of TreeSE. With this,
# user can automatically set ranks based on imported data.
.set_ranks_based_on_rowdata <- function(
        tse, set.ranks = FALSE, verbose = TRUE,
        ignore.col = "taxonomy_unparsed", ...){
    #
    if( !.is_a_bool(set.ranks) ){
        stop("'set.ranks' must be TRUE or FALSE.", call. = FALSE)
    }
    #
    if( !.is_a_bool(verbose) ){
        stop("'verbose' must be TRUE or FALSE.", call. = FALSE)
    }
    #
    if( !(is.character(ignore.col) || is.null(ignore.col)) ){
        stop("'ignore.col' must be a character value or NULL.", call. = FALSE)
    }
    #
    # Get ranks from rowData
    rd <- rowData(tse)
    # Remove those columns that are ignored. By default, the column
    # containing unparsed taxonomy.
    rd <- rd[ , !colnames(rd) %in% ignore.col, drop = FALSE]
    # Ranks must be character columns
    is_char <- lapply(
        rd, function(x) is.character(x) || is.factor(x))
    is_char <- unlist(is_char)
    rd <- rd[ , is_char, drop = FALSE]
    # If user wants to set ranks and there are ranks after filtering out
    # those columns that are not characters.
    if( set.ranks && ncol(rd) > 0L ){
        # Finally, set ranks and give message
        ranks <- colnames(rd)
        temp <- setTaxonomyRanks(ranks)
        if( verbose ){
            message(
                "TAXONOMY_RANKS set to: '",
                paste0(ranks, collapse = "', '"), "'")
        }
    }
    # If user wanted to set ranks but there were no suitable columns in rowData,
    # give warning
    if( set.ranks && ncol(rd) == 0L ){
        warning(
            "Ranks cannot be set. rowData(x) does not include columns ",
            "specifying character values.", call. = FALSE)
    }
    return(NULL)
}

################################################################################
# This function converts vector of character values to capitalized.
.capitalize <- function(x){
    paste0(toupper(substr(x, 1, 1)), substr(x, 2, nchar(x)))
}

################################################################################
# This function adds dimension reduction to reducedDim
.add_object_to_reduceddim <- function(
        tse, res, name, subset.result = TRUE, ...){
    # Test subset
    if( !.is_a_bool(subset.result) ){
        stop("'subset.result' must be TRUE or FALSE.", call. = FALSE)
    }
    #
    # If samples do not match / there were samples without appropriate metadata
    # and they are now removed
    if( !all(colnames(tse) %in% rownames(res)) && subset.result ){
        # Get samples that are being removed
        samples_rm <- setdiff(colnames(tse), rownames(res))
        # Take a subset
        tse <- tse[ , rownames(res) ]
        # Give a message
        warning("The following samples are removed from the data as they ",
                "are not includes in dimension reduction results ",
                "(see 'subset.result' parameter): '",
                paste0(samples_rm, collapse = "', '"), "'", call. = FALSE)
    } else if( !all(colnames(tse) %in% rownames(res)) && !subset.result ){
        # If user do not want to subset the data
        # Save attributes from the object as they are removed when subsetting
        attr <- attributes(res)
        attr <- attr[ !names(attr) %in% c("dim", "dimnames")]
        # Add samples that were removed
        res <- res[match(colnames(tse), rownames(res)), , drop = FALSE]
        rownames(res) <- colnames(tse)
        # Add attributes
        attr <- c(attributes(res), attr)
        attributes(res) <- attr
    }
    # Add object to reducedDIm
    reducedDim(tse, name) <- res
    return(tse)
}

#' Joint Robust PCA on Multiple Compositional Tables
#'
#' Internal engine for Joint Robust Principal Component Analysis (RPCA) using
#' OptSpace on multiple compositional tables.
#'
#' This function assumes a list of already extracted tables and is typically
#' called via \code{jointRPCAuniversal()} or \code{getJointRPCA()}.
#'
#' @param tables A list of compositional data tables (matrices or data frames).
#' @param n.test.samples Integer specifying the number of samples to hold out for testing
#'   (only used if \code{sample.metadata} is \code{NULL}). Default is 10.
#' @param sample.metadata Optional data frame containing sample-level metadata.
#' @param train.test.column The name of the column in \code{sample.metadata}
#'   that defines training vs test samples.
#' @param n.components Integer specifying the number of principal components to compute.
#' @param transform Character string specifying preprocessing applied to each
#'   input table before ordination: \code{"rclr"} or \code{"none"}.
#' @param optspace.tol Numeric tolerance passed to \code{vegan::optspace()}.
#' @param center,scale Logical; whether to center/scale the reconstructed matrix
#'   prior to SVD/PCA steps.
#' @param min.sample.count Minimum total count required for a sample to be retained.
#' @param min.feature.count Minimum total count required for a feature to be retained.
#' @param min.feature.frequency Minimum percentage (0–100) of samples in which a
#'   feature must be non-zero to be retained.
#' @param max.iterations Maximum number of optimization iterations.
#'
#' @return A list with \code{ord_res}, \code{dist}, \code{cv_stats}, and
#'   \code{rclr_tables}.
#'
#' @keywords internal
#' @noRd

.joint_rpca <- function(tables,
                        n.test.samples = 10,
                        sample.metadata = NULL,
                        train.test.column = NULL,
                        n.components = 3,
                        transform = c("rclr", "none"),
                        min.sample.count = 0,
                        min.feature.count = 0,
                        min.feature.frequency = 0,
                        max.iterations = 5,
                        optspace.tol = 1e-5,
                        center = TRUE,
                        scale = FALSE) {
    
    transform <- match.arg(transform)
    
    if (is.null(names(tables))) {
        names(tables) <- paste0("view", seq_along(tables))
    }
    
    if (n.components < 2) {
        stop("n.components must be at least 2.", call. = FALSE)
    }
    if (max.iterations < 1) {
        stop("max.iterations must be at least 1.", call. = FALSE)
    }
    
    # Filtering (always done, independent of rclr)
    tables <- lapply(tables, function(tbl) {
        out <- .rpca_table_processing(
            tbl,
            min.sample.count      = min.sample.count,
            min.feature.count     = min.feature.count,
            min.feature.frequency = min.feature.frequency
        )
        if (nrow(out) == 0L || ncol(out) == 0L) {
            stop(
                "Filtering removed all data in at least one table (0 features or 0 samples). ",
                "Try relaxing filtering thresholds (min.sample.count / min.feature.count / ",
                "min.feature.frequency) or check your input.",
                call. = FALSE
            )
        }
        out
    })
    
    # Find shared samples across views
    sample.sets <- lapply(tables, colnames)
    shared.all.samples <- Reduce(intersect, sample.sets)
    if (length(shared.all.samples) == 0L) {
        stop(
            "No samples overlap between all tables. ",
            "Check that colnames are consistent across views, or",
            "if you are using pre-transformed tables set transform = 'none'.",
            call. = FALSE
        )
    }
    if (length(shared.all.samples) < (n.components + 1L)) {
        stop(
            "Too few shared samples across all tables after filtering (",
            length(shared.all.samples), "). Need at least n.components + 1 shared samples. ",
            "Try lowering n.components or relaxing filtering thresholds.",
            call. = FALSE
        )
    }
    unshared.samples <- setdiff(unique(unlist(sample.sets)), shared.all.samples)
    if (length(unshared.samples) > 0) {
        warning(sprintf("Removing %d sample(s) that do not overlap in tables.", length(unshared.samples)))
    }
    
    # Restrict each table to the shared sample set
    tables <- lapply(tables, function(tbl) {
        tbl[, shared.all.samples, drop = FALSE]
    })
    shared.all.samples <- Reduce(intersect, lapply(tables, colnames))
    
    # Transform tables: rCLR or masking
    rclr_tables <- lapply(tables, function(tbl) {
        mat <- as.matrix(tbl)
        rown <- rownames(mat)
        coln <- colnames(mat)
        
        if (transform == "rclr") {
            mat[!is.finite(mat)] <- 0
            mat[mat < 0]         <- 0
            out <- vegan::decostand(mat, method = "rclr", MARGIN = 2)
            dimnames(out) <- list(rown, coln)
            out
        } else {
            out <- .mask_value_only(mat)$data
            dimnames(out) <- list(rown, coln)
            out
        }
    })
    names(rclr_tables) <- names(tables)
    
    # Determine train/test split
    if (!is.null(sample.metadata) && !is.null(train.test.column)) {
        md <- as.data.frame(sample.metadata)
        md <- md[shared.all.samples, , drop = FALSE]
        train.samples <- rownames(md)[md[[train.test.column]] == "train"]
        test.samples  <- rownames(md)[md[[train.test.column]] == "test"]
    } else {
        ord.tmp <- .optspace_helper(
            rclr.table      = t(rclr_tables[[1]]),
            feature.ids     = rownames(rclr_tables[[1]]),
            sample.ids      = colnames(rclr_tables[[1]]),
            n.components    = n.components,
            max.iterations  = max.iterations,
            tol             = optspace.tol,
            center          = center,
            scale           = scale
        )$ord_res
        sorted.ids <- rownames(ord.tmp$samples[order(ord.tmp$samples[, 1]), ])
        idx <- round(seq(1, length(sorted.ids), length.out = n.test.samples))
        test.samples <- sorted.ids[idx]
        train.samples <- setdiff(shared.all.samples, test.samples)
    }
    
    # Run joint OptSpace
    result <- .joint_optspace_helper(
        tables         = rclr_tables,
        n.components   = n.components,
        max.iterations = max.iterations,
        test.samples   = test.samples,
        train.samples  = train.samples,
        sample.order   = shared.all.samples,
        tol            = optspace.tol,
        center         = center,
        scale          = scale
    )
    
    return(list(
        ord_res      = result$ord_res,
        dist         = result$dist,
        cv_stats     = result$cv_stats,
        rclr_tables  = rclr_tables
    ))
}

#' Joint RPCA Ordination Across Multiple Compositional Tables
#'
#' Internal function that performs Robust PCA via joint OptSpace decomposition across multiple compositional tables.
#' It splits each table into train/test sets, applies joint factorization, reconstructs sample and feature embeddings,
#' optionally projects test samples, computes a sample distance matrix, and returns cross-validation error statistics.
#'
#' @param tables A list of compositional matrices or data frames with features as rows and samples as columns.
#' @param n.components Number of principal components to compute.
#' @param max.iterations Maximum number of optimization iterations for OptSpace.
#' @param test.samples Character vector of sample IDs to be projected into the ordination space.
#' @param train.samples Character vector of sample IDs used to fit the ordination.
#'
#' @return A list with:
#' \describe{
#'   \item{ord_res}{An \code{OrdinationResults} object containing embeddings, loadings, and variance explained.}
#'   \item{dist}{A \code{DistanceMatrix} object for sample embeddings.}
#'   \item{cv_stats}{A data frame summarizing reconstruction error across iterations and tables.}
#' }
#'
#' @keywords internal
#' @noRd

.joint_optspace_helper <- function(tables,
                                   n.components,
                                   max.iterations,
                                   test.samples,
                                   train.samples,
                                   sample.order = NULL,
                                   tol = 1e-5,
                                   center = TRUE,
                                   scale = FALSE) {
    
    # Coerce to matrices and enforce colnames presence
    tables <- lapply(tables, function(tbl) {
        mat <- as.matrix(tbl)
        if (is.null(colnames(mat))) {
            stop("[.joint_optspace_helper] Input table is missing column names (sample IDs).")
        }
        mat
    })
    
    # Global set of samples present in all views
    all_samples <- Reduce(intersect, lapply(tables, colnames))
    
    # Align train/test to actually available samples
    test.samples  <- intersect(test.samples,  all_samples)
    train.samples <- intersect(train.samples, all_samples)
    
    if (!length(test.samples) || !length(train.samples)) {
        stop("[.joint_optspace_helper] Empty train/test split after aligning sample IDs.")
    }
    
    # Split and transpose training/test data per table
    tables.split <- lapply(tables, function(tbl) {
        list(
            t(tbl[, test.samples,  drop = FALSE]),
            t(tbl[, train.samples, drop = FALSE])
        )
    })
    
    # Format input for solver
    tables.for.solver <- lapply(tables.split, function(pair) {
        lapply(pair, as.matrix)
    })
    
    # Run joint OptSpace solver
    opt.result <- .joint_optspace_solve(
        train.test.pairs = tables.for.solver,
        n.components     = n.components,
        max.iter         = max.iterations,
        tol              = tol
    )
    
    U <- opt.result$U
    S <- opt.result$S
    V_list <- opt.result$V_list
    dists <- opt.result$dists
    
    # Assign row/column names to loadings
    pc.names <- paste0("PC", seq_len(n.components))
    
    # Combine feature loadings with table-derived row names
    vjoint <- do.call(rbind, Map(function(tbl, V) {
        rownames(V) <- rownames(tbl)
        colnames(V) <- pc.names
        V
    }, tables, V_list))
    
    U <- U[seq_along(train.samples), , drop = FALSE]
    rownames(U) <- train.samples
    colnames(U) <- pc.names
    
    # Recenter & re-factor via SVD
    X <- U %*% S %*% t(vjoint)
    
    if (center) {
        X <- sweep(X, 2, colMeans(X))
        X <- sweep(X, 1, rowMeans(X))
    }
    if (scale) {
        X <- scale(X, center = FALSE, scale = TRUE)
    }
    svd.res <- svd(X)
    u <- svd.res$u[, seq_len(n.components), drop = FALSE]
    v <- svd.res$v[, seq_len(n.components), drop = FALSE]
    s.eig <- svd.res$d[seq_len(n.components)]
    
    rownames(u) <- train.samples
    rownames(v) <- rownames(vjoint)
    pc.names <- paste0("PC", seq_len(n.components))
    colnames(u) <- colnames(v) <- pc.names
    
    # Build a named per-view features list
    features_list <- lapply(seq_along(tables), function(i) {
        rid <- rownames(tables[[i]])          
        v[rid, , drop = FALSE]                
    })
    names(features_list) <- names(tables)    
    
    prop.exp <- s.eig^2 / sum(s.eig^2)
    ord_res <- .ordination_results(
        method = "rpca",
        eigvals = setNames(s.eig, pc.names),
        samples = u,
        features = features_list,              
        proportion.explained = setNames(prop.exp, pc.names)
    )
    
    # Project test samples
    if (length(test.samples) > 0) {
        test.matrices <- lapply(tables, function(tbl) tbl[, test.samples, drop = FALSE])
        names(test.matrices) <- names(tables)  
        ord_res <- .transform(ord_res, test.matrices, apply.rclr = FALSE)
    }
    
    # Compute distance matrix and CV error summary
    dist.base <- as.matrix(dist(ord_res$samples))
    
    if (!is.null(sample.order)) {
        order_use <- intersect(sample.order, rownames(dist.base))
        dist.mat  <- dist.base[order_use, order_use, drop = FALSE]
    } else {
        dist.mat  <- dist.base
        order_use <- rownames(dist.base)
    }
    
    dist.res <- .distance_matrix(dist.mat, ids = order_use)
    
    cv.dist <- data.frame(t(dists))
    colnames(cv.dist) <- c("mean_CV", "std_CV")
    cv.dist$run <- sprintf("tables_%d.n.components_%d.max.iterations_%d.n.test_%d",
                           length(tables), n.components, max.iterations, length(test.samples))
    cv.dist$iteration <- seq_len(nrow(cv.dist))
    rownames(cv.dist) <- seq_len(nrow(cv.dist))
    
    return(list(ord_res = ord_res, dist = dist.res, cv_stats = cv.dist))
}

#' Apply Projection of New Compositional Tables to Existing Ordination
#'
#' Internal function that transforms and projects new sample tables into an existing Joint RPCA ordination space.
#' It handles feature alignment, optional rCLR preprocessing, padding of missing features,
#' and merges projected samples with existing ones.
#'
#' @param ordination A list containing previous ordination results: `samples`, `features`, and `eigvals`.
#' @param tables A named list of new compositional tables (matrices or data frames)
#'   with features as rows and samples as columns.
#' @param apply.rclr Logical; whether to apply rCLR transformation to new input tables before projection. Default is `TRUE`.
#'
#' @return An updated ordination list with the `samples` matrix extended to include new projected samples.
#' @keywords internal
#' @noRd

.transform <- function(ordination, tables,
                       apply.rclr = TRUE) {
    
    Udf    <- ordination$samples
    Vobj   <- ordination$features
    s.eig  <- ordination$eigvals
    
    # Ensure tables is a list of views
    if (!is.list(tables)) {
        stop("[.transform] 'tables' must be a list of view matrices (features x samples).")
    }
    
    if (is.list(Vobj) && !is.null(names(Vobj))) {
        if (is.null(names(tables))) {
            names(tables) <- names(Vobj)[seq_along(tables)]
        }
        if (is.null(names(tables))) {
            stop("[.transform] 'tables' must be a *named* list of view matrices (features x samples).")
        }
    }
    
    # rCLR if requested 
    prep_view <- function(tab) {
        mat <- as.matrix(tab)
        rown <- rownames(mat)
        coln <- colnames(mat)
        
        if (apply.rclr) {
            mat <- vegan::decostand(mat, method = "rclr", MARGIN = 2)
            
            dimnames(mat) <- list(rown, coln)
        }
        
        storage.mode(mat) <- "double"
        mat[!is.finite(mat)] <- 0
        mat
    }
    tables <- lapply(tables, prep_view)
    
    if (is.matrix(Vobj)) {
        all.features <- rownames(Vobj)
        tables <- lapply(tables, function(mat) {
            miss <- setdiff(all.features, rownames(mat))
            if (length(miss)) {
                pad <- matrix(0, nrow = length(miss), ncol = ncol(mat),
                              dimnames = list(miss, colnames(mat)))
                mat <- rbind(mat, pad)
            }
            mat[all.features, , drop = FALSE]
        })
        proj.mat <- do.call(cbind, tables)
        colnames(proj.mat) <- make.unique(colnames(proj.mat), sep = "_")
        ordination$samples <- .transform_helper(Udf, Vobj, s.eig, proj.mat)
        return(ordination)
    }
    
    # 3+-omic path: V is a named list per view
    if (!is.list(Vobj) || is.null(names(Vobj))) {
        stop("[.transform] ordination$features is neither a matrix nor a named list.")
    }
    
    # Intersect views by name, preserve training order
    views <- intersect(names(Vobj), names(tables))
    if (!length(views)) stop("[.transform] No overlapping view names between ordination and new tables.")
    
    test.matrices <- list()
    for (vw in views) {
        Vvw <- Vobj[[vw]]
        stopifnot(is.matrix(Vvw), !is.null(rownames(Vvw)))
        mat <- tables[[vw]]
        
        train_feats <- rownames(Vvw)
        miss <- setdiff(train_feats, rownames(mat))
        if (length(miss)) {
            pad <- matrix(0, nrow = length(miss), ncol = ncol(mat),
                          dimnames = list(miss, colnames(mat)))
            mat <- rbind(mat, pad)
        }
        mat <- mat[train_feats, , drop = FALSE]
        
        test.matrices[[vw]] <- mat
    }
    
    ordination$samples <- .transform_helper(Udf, Vobj, s.eig, test.matrices)
    return(ordination)
}

#' Project New Data into Existing Ordination Space
#'
#' Internal function to align rCLR-transformed samples to an existing RPCA ordination space.
#' Handles feature alignment, deduplication of sample names, double-centering normalization,
#' and projection into low-rank space using previously learned components.
#'
#' @param Udf Matrix of training sample embeddings (samples × components).
#' @param Vdf Matrix or list of feature loadings (features × components, or per-view list).
#' @param s.eig Singular values from the RPCA decomposition.
#' @param table.rclr.project New rCLR-transformed table(s) for projection (features × samples).
#'   If `Vdf` is a matrix, provide a single matrix; if `Vdf` is a list, provide a named list of matrices per view.
#' @param dedup.samples Logical; whether to merge samples with identical names (e.g. suffixes like `_1`, `_2`)
#'   by averaging their feature values. Default is `TRUE`. Set to `FALSE` to preserve all duplicate sample IDs.
#'
#' @return A combined matrix of training and projected samples (samples × components).
#' @keywords internal
#' @noRd

.transform_helper <- function(Udf, Vdf, s.eig, table.rclr.project,
                              dedup.samples = TRUE) {
    
    # Legacy path (single view)
    if (is.matrix(Vdf)) {
        stopifnot(is.matrix(table.rclr.project))
        # Align rows by name
        common <- intersect(rownames(Vdf), rownames(table.rclr.project))
        if (length(common) < ncol(Udf))
            stop(sprintf("[.transform_helper] Too few matching features: %d", length(common)))
        
        M <- t(as.matrix(table.rclr.project[common, , drop = FALSE]))   
        V <- as.matrix(Vdf[common, , drop = FALSE])                     
        
        # Dedup of sample IDs
        if (dedup.samples) {
            sid <- sub("_\\d+$", "", rownames(M))
            if (any(duplicated(sid))) {
                M <- rowsum(M, group = sid, reorder = FALSE) / as.vector(table(sid))
            } else {
                rownames(M) <- sid
            }
        }
        
        # Projection (match training scaling)
        Uproj <- M %*% V
        # Scale by singular values
        if (length(s.eig)) {
            Sinv <- diag(1 / s.eig, nrow = length(s.eig))
            Uproj <- Uproj %*% Sinv
        }
        
        colnames(Uproj) <- colnames(Udf)
        U.combined <- rbind(Udf[setdiff(rownames(Udf), rownames(Uproj)), , drop = FALSE], Uproj)
        return(U.combined)
    }
    
    # Multi-view path (named lists)
    stopifnot(is.list(Vdf), is.list(table.rclr.project))
    views <- intersect(names(Vdf), names(table.rclr.project))
    if (!length(views)) stop("[.transform_helper] No overlapping views.")
    
    # Project per view, then sum contributions in the shared latent space
    Usum <- NULL
    ncomp <- ncol(Udf)
    for (vw in views) {
        Vvw <- Vdf[[vw]]
        Tvw <- table.rclr.project[[vw]]
        stopifnot(is.matrix(Vvw), is.matrix(Tvw))
        
        common <- intersect(rownames(Vvw), rownames(Tvw))
        if (length(common) < ncomp) {
            stop(sprintf("[.transform_helper] View '%s': too few matching features (%d).", vw, length(common)))
        }
        
        M <- t(as.matrix(Tvw[common, , drop = FALSE]))     
        V <- as.matrix(Vvw[common, , drop = FALSE])        
        
        # Accumulate per-view U
        Uvw <- M %*% V                                     
        if (is.null(Usum)) {
            Usum <- Uvw
        } else {
            # Align rows (samples) by name before summing
            all_s <- union(rownames(Usum), rownames(Uvw))
            Utmp  <- matrix(0, nrow = length(all_s), ncol = ncol(Udf),
                            dimnames = list(all_s, colnames(Udf)))
            Utmp[rownames(Usum), ] <- Usum
            Utmp[rownames(Uvw), ]  <- Utmp[rownames(Uvw), ] + Uvw
            Usum <- Utmp
        }
    }
    
    # Sample dedup (after combining views)
    if (dedup.samples) {
        sid <- sub("_\\d+$", "", rownames(Usum))
        if (any(duplicated(sid))) {
            Usum <- rowsum(Usum, group = sid, reorder = FALSE) / as.vector(table(sid))
        } else {
            rownames(Usum) <- sid
        }
    }
    
    # Scale by S
    if (length(s.eig)) {
        Sinv <- diag(1 / s.eig, nrow = length(s.eig))
        Usum <- Usum %*% Sinv
    }
    colnames(Usum) <- colnames(Udf)
    
    # Merge with training U, avoiding duplicates
    keep_train <- setdiff(rownames(Udf), rownames(Usum))
    rbind(Udf[keep_train, , drop = FALSE], Usum)
}

#' RPCA Table Filtering and Preprocessing
#'
#' Internal function that performs filtering and cleanup on a compositional data table
#' prior to Robust PCA analysis. Removes low-count features/samples, enforces non-zero frequency thresholds,
#' checks for ID duplication, and returns a matrix suitable for transformation and ordination.
#'
#' @param table A matrix or data frame with features as rows and samples as columns.
#' @param min.sample.count Minimum total count required for a sample to be retained. Default is 0.
#' @param min.feature.count Minimum total count required for a feature to be retained. Default is 0.
#' @param min.feature.frequency Minimum percentage (0–100) of samples in which a feature must be non-zero. Default is 0.
#'
#' @return A filtered numeric matrix containing non-empty features and samples.
#' @keywords internal
#' @noRd

.rpca_table_processing <- function(table,
                                   min.sample.count = 0,
                                   min.feature.count = 0,
                                   min.feature.frequency = 0) {
    # Ensure the input is a matrix
    if (is.data.frame(table)) {
        table <- as.matrix(table)
    }
    
    n.features <- nrow(table)
    n.samples  <- ncol(table)
    
    # Filter features by total count
    if (!is.null(min.feature.count)) {
        feature.totals <- rowSums(table, na.rm = TRUE)
        keep.features <- feature.totals > min.feature.count
        table <- table[keep.features, , drop = FALSE]
    }
    
    # Filter features by frequency across samples
    if (!is.null(min.feature.frequency)) {
        freq.threshold <- min.feature.frequency / 100
        feature.freq <- rowMeans(table > 0, na.rm = TRUE)
        keep.features <- feature.freq > freq.threshold
        table <- table[keep.features, , drop = FALSE]
    }
    
    # Filter samples by total count
    if (!is.null(min.sample.count)) {
        sample.totals <- colSums(table, na.rm = TRUE)
        keep.samples <- sample.totals > min.sample.count
        table <- table[, keep.samples, drop = FALSE]
    }
    
    # Check for duplicate IDs
    if (any(duplicated(colnames(table)))) {
        stop("Data table contains duplicate sample (column) IDs.", call. = FALSE)
    }
    if (any(duplicated(rownames(table)))) {
        stop("Data table contains duplicate feature (row) IDs.", call. = FALSE)
    }
    
    # Remove empty rows and columns if sample filtering applied
    if (!is.null(min.sample.count)) {
        nonzero.features <- rowSums(table, na.rm = TRUE) > 0
        nonzero.samples  <- colSums(table, na.rm = TRUE) > 0
        table <- table[nonzero.features, nonzero.samples, drop = FALSE]
    }
    
    return(table)
}

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
#' @noRd

.mask_value_only <- function(mat) {
    # Ensure matrix is at least 2D
    if (is.vector(mat)) {
        mat <- matrix(mat, nrow = 1)
    }
    
    # Ensure matrix is not more than 2D
    if (length(dim(mat)) > 2) {
        stop("Input matrix can only have two dimensions or less")
    }
    
    # Generate logical mask: TRUE where values are missing
    mask <- !is.finite(mat)  
    
    # Create masked matrix
    masked.mat <- mat
    masked.mat[!is.finite(mat)] <- NA
    
    # Return as a masked matrix
    return(structure(list(
        data = masked.mat,
        mask = mask
    ), class = "MaskedMatrix"))
}

#' Internal constructor for ordination results
#' @keywords internal
#' @noRd
.ordination_results <- function(method, eigvals, samples, features,
                               proportion.explained, dist = NULL, metadata = list()) {
    return(structure(list(
        method = method,
        eigvals = eigvals,
        samples = samples,
        features = features,
        proportion.explained = proportion.explained,
        dist = dist,
        metadata = metadata
    ), class = "OrdinationResults"))
}

#' Internal constructor for a distance matrix object
#' @keywords internal
#' @noRd
.distance_matrix <- function(matrix, ids = NULL, method = "euclidean") {
    if (!is.matrix(matrix)) stop("Input must be a matrix.")
    if (!isSymmetric(matrix)) stop("Distance matrix must be symmetric.")
    if (!is.null(ids)) {
        if (length(ids) != nrow(matrix)) stop("Length of 'ids' must match matrix dimensions.")
        rownames(matrix) <- ids
        colnames(matrix) <- ids
    }
    return(structure(list(
        data = matrix,
        ids = rownames(matrix),
        method = method
    ), class = "DistanceMatrix"))
}

#' Extract per-experiment assay tables from a MultiAssayExperiment
#'
#' @param x A MultiAssayExperiment.
#' @param experiments Character vector of experiment names (or NULL for all).
#'
#' @return A list with `tables`, `experiments`, and `assay_names_used`.
#'
#' @keywords internal
#' @noRd
.extract_mae_tables <- function(x, experiments = NULL) {
    exps <- experiments(x)
    
    if (is.null(experiments)) {
        experiments <- names(exps)
    }
    if (length(experiments) == 0L) {
        stop("No experiments found in 'x'.", call. = FALSE)
    }
    
    assay_names_used <- setNames(character(length(experiments)), experiments)
    tables <- vector("list", length(experiments))
    names(tables) <- experiments
    
    for (i in seq_along(experiments)) {
        e <- experiments[[i]]
        exp_se <- exps[[e]]
        if (is.null(exp_se)) {
            stop(sprintf("Experiment '%s' not found in 'x'.", e), call. = FALSE)
        }
        
        anm <- assayNames(exp_se)
        default_assay <- if (length(anm)) anm[[1]] else 1L
        
        assay_names_used[[e]] <- if (is.character(default_assay)) default_assay else as.character(default_assay)
        tables[[e]] <- assay(exp_se, default_assay)
    }
    
    list(
        tables = tables,
        experiments = experiments,
        assay_names_used = assay_names_used
    )
}