#' Add information on microbial modules
#'
#' \code{getModules} and \code{addModules} generate a modules table from a list
#' of microbial signatures.
#'
#' @param x a
#' \code{\link[TreeSummarizedExperiment:TreeSummarizedExperiment-constructor]{TreeSummarizedExperiment}}.
#'
#' @param sigs \code{Character list}. List of microbial signatures in the
#' metaphlan taxonomy format.
#'
#' @param exact.tax.level \code{Logical scalar}. Should only the last
#' taxonomic rank be used to determine whether a feature belongs to a module?
#' if \code{FALSE}, all ranks are considered. (Default: \code{FALSE}).
#' 
#' @details
#' \code{getModules} and \code{addModules} allow to integrate information on
#' microbial modules into \code{rowData(x)} from databases like BugSigDB or
#' custom tables of microbes like \code{butyrate}. \code{sigs} must be a single
#' character vector or a list of character vectors, where each vector
#' corresponds to a module. Each microbial signature should follow the metaphlan
#' taxonomy format, e.g.,
#' \code{"k__Bacteria|p__Actinobacteria|c__Actinomycetia|o__Corynebacteriales"}.
#' 
#' @return
#' \code{getModules} returns a modules table where rows and columns represent
#' features and modules, respectively. \code{addModules} returns an updated
#' \code{\link[TreeSummarizedExperiment:TreeSummarizedExperiment-constructor]{TreeSummarizedExperiment}}
#' object with the modules table in the \code{rowData} slot.
#' 
#' @name getModules
#'
#' @examples
#' \dontrun{
#' # Load butyrate module
#' data("butyrate", package = "ariadne")
#' 
#' # Load dataset
#' data("Tengeler2020", package = "mia")
#' tse <- Tengeler2020
#' 
#' # Convert taxonomy to fulltax labels (add to getTaxonomyLabels?)
#' sigs <- mia:::.tax_table2label(butyrate)
#' 
#' # Generate modules table
#' modules <- getModules(tse, sigs)
#' 
#' # Generate & store modules table
#' tse <- addModules(tse, sigs)
#' }
NULL

#' @rdname getModules
#' @export
#' @importFrom SummarizedExperiment rowData
setMethod("addModules", signature = c(x = "TreeSummarizedExperiment"),
    function(x, sigs, exact.tax.level = FALSE){
        # Make modules table
        modules <- getModules(
            x,
            sigs,
            exact.tax.level = exact.tax.level
        )
        # Bind modules table with rowData
        rowData(x) <- cbind(rowData(x), modules)
        return(x)
    }
)

#' @rdname getModules
#' @export
setMethod("getModules", signature = c(x = "TreeSummarizedExperiment"),
    function(x, sigs, exact.tax.level = FALSE){
        # Check sigs
        if( !is.vector(sigs) ){
            stop("sigs must be a character vector or list of character ",
                "vectors, where each vector corresponds to a module",
                call. = FALSE)
        }
        # Check exact.tax.level
        if( !.is_a_bool(exact.tax.level) ){
            stop("verbose must be TRUE or FALSE.", call. = FALSE)
        }
        # Make modules table
        modules <- .make_modules_table(
            x,
            sigs,
            exact.tax.level = exact.tax.level
        )
        return(modules)
    }
)

# Define function to construct modules table based on bugsigdb signatures
#' @importFrom stringr fixed str_detect str_escape str_remove
.make_modules_table <- function(x, sigs, exact.tax.level = FALSE){
    # Retrieve taxonomic labels for features
    tax.labs <- .tax_table2label(rowData(x)[taxonomyRanks(x)])
    # Reduce to deepest rank if exact.tax.level is on
    if( exact.tax.level ){
        tax.labs <- str_remove(tax.labs, ".*\\|")
    }
    # Convert sigs to list in case of only one module
    if( !is(sigs, "list") ){
        sigs <- list(module = sigs)
    }
    # Initialise all-FALSE modules table
    modules <- matrix(
        FALSE,
        nrow = length(tax.labs),
        ncol = length(sigs),
        dimnames = list(rownames(x), names(sigs))
    )
    # For every signature in modules list
    for( i in seq_along(sigs) ){
        # Extract deepest taxonomic rank
        sig <- sigs[[i]] |>
            str_remove(".*\\|") |>
            str_escape() |>
            paste0(collapse = "|")
        # Find which features belong to the current signature
        modules[ , i] <- str_detect(tax.labs, sig)
    }
    return(modules)
}

### HELPER FUNCTIONS ###

# Add taxrank prefixes to taxonomy table
.add_prefix_to_taxtable <- function(tax){
    # Retrieve taxrank prefixes
    tax.prefix <- paste0(getTaxonomyRankPrefixes()[tolower(names(tax))], "__")
    # Add prefix unless rank is NA
    tax <- mapply(function(rank, prefix)
        ifelse(is.na(rank), NA, paste0(prefix, rank)), tax, tax.prefix)
    return(tax)
}

# Reduce taxcols of rowData to taxstring in metaphlan format
#' @importFrom SummarizedExperiment rowData
.tax_table2label <- function(x){
    # Add taxrank prefixes to taxcols of rowData
    tax <- .add_prefix_to_taxtable(x)
    # Collapse taxcols to taxstring in metaphlan format
    tax <- apply(tax, 1L, function(row) paste(row, collapse = "|"))
    # Remove empty taxranks
    tax <- gsub("(?:\\|[a-z]__)+$", "", tax)
    return(tax)
}
