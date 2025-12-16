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

#ORDINATION RESULTS FUNCTION
#
#constructor function
.OrdinationResults <- function(method, eigvals, samples, features,
                               proportion.explained, dist = NULL, metadata = list()) {
    structure(list(
        method = method,
        eigvals = eigvals,
        samples = samples,
        features = features,
        proportion.explained = proportion.explained,
        dist = dist,
        metadata = metadata
    ), class = "OrdinationResults")
}

#print method
.print.OrdinationResults <- function(x, ...) {
    cat("OrdinationResults (method:", x$method, ")\n")
    cat("Number of components:", length(x$eigvals), "\n")
    cat("Variance explained:\n")
    print(round(x$proportion.explained, 3))
    invisible(x)
}

#summary method
.summary.OrdinationResults <- function(object, ...) {
    print(object)
    cat("\nSample scores (first few rows):\n")
    print(head(object$samples))
    cat("\nFeature loadings (first few rows):\n")
    print(head(object$features))
    invisible(object)
}

#plot method
.plot.OrdinationResults <- function(x, comps = c(1, 2), ...) {
    if (length(comps) != 2) stop("Please select two components to plot.")
    plot(x$samples[, comps], col = "blue", pch = 19,
         xlab = paste0("PC", comps[1]),
         ylab = paste0("PC", comps[2]),
         main = paste("Ordination (", x$method, ")", sep = ""))
    points(x$features[, comps], col = "red", pch = 4)
    legend("topright", legend = c("Samples", "Features"),
           col = c("blue", "red"), pch = c(19, 4))
}

#DISTANCE MATRIX FUNCTION
#
.DistanceMatrix <- function(matrix, ids = NULL, method = "euclidean") {
    if (!is.matrix(matrix)) stop("Input must be a matrix.")
    if (!isSymmetric(matrix)) stop("Distance matrix must be symmetric.")
    if (!is.null(ids)) {
        if (length(ids) != nrow(matrix)) stop("Length of 'ids' must match matrix dimensions.")
        rownames(matrix) <- ids
        colnames(matrix) <- ids
    }
    structure(list(
        data = matrix,
        ids = rownames(matrix),
        method = method
    ), class = "DistanceMatrix")
}

.print.DistanceMatrix <- function(x, ...) {
    cat("DistanceMatrix (", x$method, ")\n", sep = "")
    cat("Number of objects:", length(x$ids), "\n")
    print(head(x$data, 6))  # Show only top part
    invisible(x)
}

.summary.DistanceMatrix <- function(object, ...) {
    cat("Summary of DistanceMatrix\n")
    cat("Method:", object$method, "\n")
    cat("Size:", nrow(object$data), "x", ncol(object$data), "\n")
    cat("IDs:\n")
    print(head(object$ids, 6))
    cat("\nDistance Summary Stats:\n")
    print(summary(as.vector(object$data[upper.tri(object$data)])))
    invisible(object)
}

#wrapper to store dataset-specific sample scores
.dataset_specific_scores <- function(rclr.tables, n.components = 2, max.iterations = 5) {
    scores <- lapply(seq_along(rclr.tables), function(i) {
        tbl <- rclr.tables[[i]]
        res <- .optspace_helper(
            rclr.table     = t(tbl),
            feature.ids    = rownames(tbl),
            subject.ids    = colnames(tbl),
            n.components   = n.components,
            max.iterations = max.iterations
        )
        res$ord.res$samples
    })
    
    names(scores) <- paste0("Dataset_", seq_along(scores))
    
    return(scores)
}

#wrapper to store dataset-specific feature loadings
.dataset_specific_loadings <- function(rclr.tables, n.components = 2, max.iterations = 5) {
    loadings <- lapply(seq_along(rclr.tables), function(i) {
        tbl <- rclr.tables[[i]]
        res <- .optspace_helper(
            rclr.table     = t(tbl),
            feature.ids    = rownames(tbl),
            subject.ids    = colnames(tbl),
            n.components   = n.components,
            max.iterations = max.iterations
        )
        res$ord.res$features
    })
    
    names(loadings) <- paste0("Dataset_", seq_along(loadings))
    
    return(loadings)
}

################################################################################
# Joint-RPCA benchmarking helpers

.keep_finite_cols <- function(X) {
    ok <- apply(X, 2, function(v) all(is.finite(v)))
    if (!any(ok)) stop("All columns removed by finite filter.")
    X[, ok, drop = FALSE]
}

.drop_constant_cols <- function(X) {
    sds <- apply(X, 2, function(v) sd(v, na.rm = TRUE))
    keep <- is.finite(sds) & (sds > 0)
    if (!any(keep)) stop("No non-constant columns remain after filtering.")
    X[, keep, drop = FALSE]
}

.prep_train_test <- function(X_train, X_test) {
    m <- colMeans(X_train, na.rm = TRUE)
    s <- apply(X_train, 2, sd, na.rm = TRUE)
    s[s == 0 | !is.finite(s)] <- 1
    list(
        Xtr = sweep(sweep(X_train, 2, m, "-"), 2, s, "/"),
        Xte = sweep(sweep(X_test,  2, m, "-"), 2, s, "/")
    )
}

.evaluate_model_cv <- function(features, labels, folds = 5, ntree = 500, seed = 42) {
    set.seed(seed)
    tab <- table(labels)
    if (length(labels) < 2L || length(tab) < 2L) stop("Need >=2 samples and >=2 classes.")
    folds <- max(2L, min(as.integer(folds), as.integer(min(tab)), length(labels) - 1L))
    folds_idx <- caret::createFolds(labels, k = folds, list = TRUE, returnTrain = FALSE)
    
    accs <- numeric(length(folds_idx)); aucs <- numeric(length(folds_idx))
    for (i in seq_along(folds_idx)) {
        test_idx  <- folds_idx[[i]]
        train_idx <- setdiff(seq_along(labels), test_idx)
        Xtr <- features[train_idx, , drop = FALSE]; Xte <- features[test_idx, , drop = FALSE]
        ytr <- labels[train_idx]; yte <- labels[test_idx]
        if (length(unique(ytr)) < 2L) { accs[i] <- NA_real_; aucs[i] <- NA_real_; next }
        
        pp <- .prep_train_test(Xtr, Xte)
        rf <- randomForest(x = pp$Xtr, y = ytr, ntree = ntree)
        
        yhat <- predict(rf, pp$Xte, type = "response")
        accs[i] <- mean(yhat == yte)
        
        probs <- predict(rf, pp$Xte, type = "prob")
        all_lvls <- levels(labels)
        miss <- setdiff(all_lvls, colnames(probs))
        if (length(miss)) for (mm in miss) probs <- cbind(probs, setNames(rep(0, nrow(probs)), mm))
        probs <- probs[, all_lvls, drop = FALSE]
        
        if (length(unique(yte)) < 2L) {
            aucs[i] <- NA_real_
        } else {
            aucs[i] <- tryCatch(as.numeric(pROC::multiclass.roc(yte, probs)$auc), error = function(e) NA_real_)
        }
    }
    list(accuracy = mean(accs, na.rm = TRUE), auc = mean(aucs, na.rm = TRUE))
    
}

.get_fold_metrics <- function(X) {
    set.seed(42)
    idx <- caret::createFolds(labels, k = safe_k, list = TRUE, returnTrain = FALSE)
    acc <- auc <- numeric(length(idx))
    for (i in seq_along(idx)) {
        te <- idx[[i]]; tr <- setdiff(seq_along(labels), te)
        if (length(unique(labels[tr])) < 2L) {acc[i] <- NA; auc[i] <- NA; next}
        pp <- .prep_train_test(X[tr, , drop = FALSE], X[te, , drop = FALSE])
        rf <- randomForest(pp$Xtr, labels[tr], ntree = 500)
        yhat <- predict(rf, pp$Xte)
        acc[i] <- mean(yhat == labels[te])
        probs <- predict(rf, pp$Xte, type = "prob")
        miss <- setdiff(levels(labels), colnames(probs))
        if (length(miss)) for (mm in miss) probs <- cbind(probs, setNames(rep(0, nrow(probs)), mm))
        probs <- probs[ , levels(labels), drop = FALSE]
        auc[i] <- tryCatch(as.numeric(pROC::multiclass.roc(labels[te], probs)$auc), error = function(e) NA)
    }
    tibble::tibble(Fold = seq_along(idx), Accuracy = acc, MacroAUC = auc)
}

.rep_dim <- function(X) ncol(X)

.timeit <- function(expr) { t0 <- proc.time(); force(expr); as.numeric((proc.time()-t0)["elapsed"]) }

.ci95 <- function(x){ x <- x[is.finite(x)]; m <- mean(x); s <- sd(x); n <- length(x); if(n <= 1||!is.finite(s)||s == 0) c(m,m,m) else c(m, m-1.96*s/sqrt(n), m+1.96*s/sqrt(n)) }

#helper to evaluate a score matrix fairly
.eval_method <- function(U_scores, meta, prefix = "AX") {
    U <- as.data.frame(U_scores)
    k <- ncol(U); colnames(U) <- paste0(prefix, seq_len(k))
    U$sample_id <- rownames(U_scores)
    
    meta_tmp <- meta
    if ("sample_id" %in% colnames(meta_tmp)) {
        meta_tmp$sample_id <- as.character(meta_tmp$sample_id)
    } else {
        meta_tmp <- tibble::rownames_to_column(meta_tmp, "sample_id")
    }
    df <- dplyr::left_join(U, meta_tmp, by = "sample_id")
    
    #Wilcoxon on first two axes
    w1p <- if (k >= 1)
        suppressWarnings(wilcox.test(df[[paste0(prefix, 1)]] ~ df$Group, exact = FALSE)$p.value)
    else NA_real_
    w2p <- if (k >= 2)
        suppressWarnings(wilcox.test(df[[paste0(prefix, 2)]] ~ df$Group, exact = FALSE)$p.value)
    else NA_real_
    
    #PERMANOVA on first up to 3 axes
    axes <- paste0(prefix, seq_len(min(3, k)))
    perm_R2 <- perm_F <- perm_p <- NA_real_
    if (length(axes) >= 2) {
        perm <- vegan::adonis2(df[, axes] ~ Group, data = df, method = "euclidean")
        perm_R2 <- perm$R2[1]; perm_F <- perm$F[1]; perm_p <- perm$`Pr(>F)`[1]
    }
    
    #weighted RF AUROC
    aucv <- NA_real_
    if (length(axes) >= 2) {
        rf_df <- na.omit(df[, c("Group", axes)])
        rf_df$Group <- factor(rf_df$Group, levels = c("non-IBD", "IBD"))
        if (nlevels(rf_df$Group) == 2 && all(table(rf_df$Group) >= 5)) {
            cls_tab <- table(rf_df$Group)
            wts <- as.numeric(1 / cls_tab); names(wts) <- names(cls_tab)
            set.seed(42)
            rf_prob <- ranger(
                Group ~ ., data = rf_df,
                num.trees     = 1000,
                probability   = TRUE,
                class.weights = wts,
                oob.error     = TRUE
            )
            p_ibd <- rf_prob$predictions[, "IBD"]
            roc_obj <- pROC::roc(rf_df$Group, p_ibd, levels = c("non-IBD", "IBD"))
            aucv <- as.numeric(pROC::auc(roc_obj))
        }
    }
    
    tibble(
        method       = prefix,
        wilcox_PC1_p = w1p,
        wilcox_PC2_p = w2p,
        permanova_R2 = perm_R2,
        permanova_F  = perm_F,
        permanova_p  = perm_p,
        AUROC        = aucv
    )
}

.lr <- function(mat, top, bot, pcnt = 0.5) {
    log(
        (colSums(mat[top, , drop = FALSE]) + pcnt) /
            (colSums(mat[bot, , drop = FALSE]) + pcnt)
    )
}

.make_groups_autodetect <- function(meta_df, sample_ids, min_frac = 0.01, min_abs = 10L) {
    out <- data.frame(
        sample_id = sample_ids,
        Group = factor(NA, levels = c("IBD", "non-IBD"))
    )
    if (is.null(meta_df) || !nrow(meta_df)) return(out)
    
    md <- as.data.frame(meta_df, stringsAsFactors = FALSE)
    names(md) <- tolower(trimws(names(md)))
    
    sid <- tolower(trimws(as.character(sample_ids)))
    thresh <- max(min_abs, floor(length(sid) * min_frac))
    overlaps <- vapply(md, function(col) {
        x <- tolower(trimws(as.character(col)))
        sum(!is.na(x) & x %in% sid)
    }, FUN.VALUE = integer(1))
    max_ov <- suppressWarnings(max(overlaps, na.rm = TRUE))
    if (!is.finite(max_ov) || max_ov < thresh) return(out)
    best <- names(overlaps)[which.max(overlaps)]
    
    if (!"diagnosis" %in% names(md)) return(out)
    
    dx  <- tolower(trimws(as.character(md$diagnosis)))
    grp <- ifelse(grepl("\\b(uc|cd|ibd)\\b", dx), "IBD",
                  ifelse(grepl("^\\s*non", dx), "non-IBD", NA_character_))
    
    md$sample_id <- tolower(trimws(as.character(md[[best]])))
    md$Group <- factor(grp, levels = c("IBD", "non-IBD"))
    join_tbl <- unique(md[, c("sample_id", "Group")])
    
    joined <- dplyr::left_join(
        data.frame(sample_id = sid, stringsAsFactors = FALSE),
        join_tbl, by = "sample_id"
    )
    joined$sample_id <- sample_ids
    joined
}

.eval_scores <- function(scores_df) {
    out <- list()
    if (!("Group" %in% names(scores_df))) return(out)
    
    #detect available component columns
    comp_cols <- grep("^V\\d+$", names(scores_df), value = TRUE)
    if (!length(comp_cols)) return(out)
    use_cols <- comp_cols[seq_len(min(3L, length(comp_cols)))]
    
    #Wilcoxon tests
    if ("V1" %in% names(scores_df)) {
        res_w1 <- try(wilcox.test(scores_df$V1 ~ scores_df$Group, exact = FALSE), silent = TRUE)
        out$wilcox_PC1_p <- if (!inherits(res_w1, "try-error")) res_w1$p.value else NA_real_
    }
    if ("V2" %in% names(scores_df)) {
        res_w2 <- try(wilcox.test(scores_df$V2 ~ scores_df$Group, exact = FALSE), silent = TRUE)
        out$wilcox_PC2_p <- if (!inherits(res_w2, "try-error")) res_w2$p.value else NA_real_
    }
    
    #PERMANOVA (only if >=2 components exist)
    if (length(use_cols) >= 2) {
        perm_df <- na.omit(scores_df[, c("Group", use_cols), drop = FALSE])
        if (nrow(perm_df) > 5 &&
            is.factor(perm_df$Group) &&
            nlevels(perm_df$Group) >= 2 &&
            all(table(perm_df$Group) >= 3)) {
            
            comp_mat <- as.matrix(perm_df[, use_cols, drop = FALSE])
            colnames(comp_mat) <- use_cols
            
            perm <- try(
                vegan::adonis2(comp_mat ~ Group, data = perm_df, method = "euclidean"),
                silent = TRUE
            )
            
            if (!inherits(perm, "try-error")) {
                out$permanova_R2 <- perm$R2[1]
                out$permanova_F  <- perm$F[1]
                out$permanova_p  <- perm$`Pr(>F)`[1]
            } else {
                out$permanova_R2 <- NA_real_
                out$permanova_F  <- NA_real_
                out$permanova_p  <- NA_real_
            }
        }
    }
    
    #AUROC with ranger if at least 1 component exists
    rf_df <- na.omit(scores_df[, c("Group", use_cols), drop = FALSE])
    if (nrow(rf_df) && is.factor(rf_df$Group) && nlevels(rf_df$Group) >= 2) {
        cls_tab <- table(rf_df$Group)
        wts <- as.numeric(1 / cls_tab)
        names(wts) <- names(cls_tab)
        set.seed(42)
        rf_prob <- ranger::ranger(
            Group ~ ., data = rf_df, num.trees = 500,
            probability = TRUE, class.weights = wts, oob.error = TRUE
        )
        if ("IBD" %in% colnames(rf_prob$predictions)) {
            p_ibd <- rf_prob$predictions[, "IBD"]
            roc_obj <- pROC::roc(rf_df$Group, p_ibd, levels = c("non-IBD", "IBD"))
            out$AUROC <- as.numeric(pROC::auc(roc_obj))
        }
    }
    
    out
}

#convenience: run Joint-RPCA and return scores + metrics
.fit_and_score <- function(mae, k, grp_df) {
    set.seed(42)
    mae <- runJointRPCA(
        x = mae,
        n.components = k,
        max.iterations = 5,
        rclr.transform.tables = TRUE,
        min.sample.count = 1,
        min.feature.count = 0,
        min.feature.frequency = 0
    )
    
    fit <- metadata(mae)$JointRPCA[["JointRPCA"]]
    
    U <- as.data.frame(fit$ord.res$samples)
    colnames(U) <- paste0("V", seq_len(ncol(U)))
    U$sample_id <- rownames(U)
    U2 <- dplyr::left_join(U, grp_df, by = "sample_id")
    
    list(
        scores  = U2,
        metrics = .eval_scores(U2),
        fit     = fit
    )
}

.find_subject_col <- function(meta) {
    nm <- tolower(trimws(names(meta)))
    hits <- c(
        "participant.id", "participant_id", "participantid",
        "subject", "subject_id", "host_subject_id", "host.subject.id",
        "participant", "host_subject"
    )
    ix <- intersect(nm, hits)
    if (length(ix)) ix[1] else NULL
}

.build_sample_subject_map <- function(meta_df, sample_ids, min_frac = 0.01, min_abs = 10L) {
    if (is.null(meta_df) || !nrow(meta_df)) return(NULL)
    md <- as.data.frame(meta_df, stringsAsFactors = FALSE)
    names(md) <- tolower(trimws(names(md)))
    
    sid <- tolower(trimws(as.character(sample_ids)))
    thresh <- max(min_abs, floor(length(sid) * min_frac))
    overlaps <- vapply(md, function(col) {
        x <- tolower(trimws(as.character(col)))
        sum(!is.na(x) & x %in% sid)
    }, FUN.VALUE = integer(1))
    max_ov <- suppressWarnings(max(overlaps, na.rm = TRUE))
    if (!is.finite(max_ov) || max_ov < thresh) return(NULL)
    best_id_col <- names(overlaps)[which.max(overlaps)]
    
    subj_col <- .find_subject_col(md)
    if (is.null(subj_col)) return(NULL)
    
    md$sample_id  <- tolower(trimws(as.character(md[[best_id_col]])))
    md$subject_id <- as.character(md[[subj_col]])
    out <- unique(md[, c("sample_id", "subject_id")])
    out <- out[!is.na(out$sample_id) & nzchar(out$sample_id) &
                   !is.na(out$subject_id) & nzchar(out$subject_id),
               , drop = FALSE]
    if (!nrow(out)) return(NULL)
    out
}

.clr_transform <- function(mat, pseudo = 1e-6) {
    x <- log(mat + pseudo)
    x <- sweep(x, 2, colMeans(x), FUN = "-")
    x[!is.finite(x)] <- 0
    x
}

.zscore_rows <- function(mat) {
    m <- rowMeans(mat)
    s <- matrixStats::rowSds(mat)
    s[s == 0 | !is.finite(s)] <- 1
    sweep(sweep(mat, 1, m, "-"), 1, s, "/")
}

.hellinger <- function(mat) {
    cs <- colSums(mat)
    cs[cs <= 0 | !is.finite(cs)] <- 1
    p  <- sweep(mat, 2, cs, "/")
    x  <- sqrt(p)
    x[!is.finite(x)] <- 0
    x
}

.make_scores_df <- function(S, sample_ids) {
    S <- as.matrix(S)
    colnames(S) <- paste0("V", seq_len(ncol(S)))
    out <- as.data.frame(S)
    out$sample_id <- sample_ids
    dplyr::left_join(out, grp_df, by = "sample_id")
}

.eval_wrapper <- function(scores_df, label) {
    list(model = label, metrics = .eval_scores(scores_df), scores = scores_df)
}

#collect metrics
.grab <- function(x) {
    m <- x$metrics
    get_num <- function(z) if (is.null(z)) NA_real_ else as.numeric(z)
    c(
        wilcox_PC1_p = get_num(m$wilcox_PC1_p),
        wilcox_PC2_p = get_num(m$wilcox_PC2_p),
        permanova_R2 = get_num(m$permanova_R2),
        permanova_p  = get_num(m$permanova_p),
        AUROC        = get_num(m$AUROC)
    )
}

.make_mae_for <- function(cols) {
    cd2 <- S4Vectors::DataFrame(row.names = samps[cols])
    se_mgx2 <- SummarizedExperiment::SummarizedExperiment(
        list(counts = X_mgx[, cols, drop = FALSE]),
        colData = cd2
    )
    se_mtx2 <- SummarizedExperiment::SummarizedExperiment(
        list(counts = X_mtx[, cols, drop = FALSE]),
        colData = cd2
    )
    mae2 <- MultiAssayExperiment::MultiAssayExperiment(list(MGX = se_mgx2, MTX = se_mtx2))
    MultiAssayExperiment::intersectColumns(mae2)
}

.plt_ord <- function(scores, title) {
    ggplot2::ggplot(scores, ggplot2::aes(V1, V2, color = Group)) +
        ggplot2::geom_point(alpha = 0.8, size = 1.1) +
        ggplot2::labs(title = title, x = "PC1", y = "PC2", color = NULL) +
        ggplot2::theme_minimal()
}

.plot_comp <- function(obj) {
    if (is.null(obj)) return(invisible(NULL))
    ggplot2::ggplot(obj$scores, ggplot2::aes(V1, V2, color = Group)) +
        ggplot2::geom_point(alpha = 0.8, size = 1.0) +
        ggplot2::labs(title = obj$model, x = "Comp1", y = "Comp2", color = NULL) +
        ggplot2::theme_minimal()
}

.get_view <- function(obj, keys) {
    if (is.null(obj)) return(NULL)
    if (is.list(obj)) {
        for (k in keys) {
            if (!is.null(obj[[k]])) return(as.data.frame(obj[[k]]))
        }
        return(NULL)
    }
    if (is.matrix(obj) || is.data.frame(obj)) return(as.data.frame(obj))
    NULL
}

.show_top <- function(V, label, top_k = 15) { 
    if (is.null(V) || !ncol(V)) return(invisible(NULL)) 
    ord <- order(V[, 1], decreasing = TRUE)
    top_ix <- seq_len(min(top_k, length(ord))) 
    cat(sprintf("\nTop %s features on PC1:\n", label)) 
    print(data.frame( 
        feature = rownames(V)[ord][top_ix], 
        loading = V[ord, 1][top_ix] 
    ))
}