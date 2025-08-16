#' Add information on microbial modules
#'
#' \code{getModules} and \code{addModules} generate a modules table from a list
#' of microbial signatures.
#'
#' @param x a
#' \code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}.
#'
#' @param sigs \code{Character list}. List of microbial signatures output of
#' bugsigdbr getSignatures in the metaphlan format.
#'
#' @param exact.tax.level \code{Logical scalar}. Should only the last
#' taxonomic rank be used to determine whether a feature belongs to a module?
#' if \code{FALSE}, all ranks are considered. (Default: \code{FALSE}).
#' 
#' @param ... additional parameters.
#' 
#' @details
#' text
#' 
#' @return
#' \code{getModules} returns a modules table where rows and columns represent
#' features and modules, respectively. \code{addModules} returns an updated
#' object of the same class as \code{x} with the modules table in the
#' \code{rowData} slot.
#' 
#' @name getModules
#'
#' @examples
#' data("Tengeler2020", package = "mia")
#' tse <- Tengeler2020
#' 
NULL

#' @rdname getModules
#' @export
#' @importFrom SummarizedExperiment rowData
setMethod("addModules", signature = c(x = "SummarizedExperiment"),
    function(x, sigs, exact.tax.level = FALSE, ...){
      
        modules <- getModules(
            x,
            sigs,
            exact.tax.level = exact.tax.level,
            ...
        )
        
        rowData(x) <- cbind(rowData(x), modules)
        return(x)
    }
)

#' @rdname getModules
#' @export
setMethod("getModules", signature = c(x = "SummarizedExperiment"),
    function(x, sigs, exact.tax.level = FALSE, ...){
        
        modules <- .make_modules_table(
            x,
            sigs,
            exact.tax.level = exact.tax.level,
            ...
        )
        
        return(modules)
    }
)

#' @importFrom stringr str_detect
# Define function to construct modules table based on bugsigdb signatures
.make_modules_table <- function(x, sigs, exact.tax.level = FALSE, ...){
    # Retrieve prefixes for taxonomic ranks
    tax.ranks <- getTaxonomyRankPrefixes()
    # Retrieve taxonomic labels for features
    tax.labs <- .rowdata2taxonomy(x)
    # Reduce to deepest rank if exact.tax.level is on
    if( exact.tax.level ){
        tax.labs <- gsub(".*\\|", "", tax.labs)
    }
    # Initialise empty list for signatures
    sig.list <- list()
    # For every signature in output from bugsigdbr::getSignatures
    for( i in seq_along(sigs) ){
        # Extract deepest taxonomic rank
        sig <- gsub(".*\\|", "", sigs[[i]])
        # Find which features belong to the current signature
        sig.list[[i]] <- unlist(lapply(tax.labs, function(x) any(str_detect(x, sig))))
    }
    # Build modules table from signature list
    modules <- do.call(cbind, sig.list)
    # Add names to rows and columns
    rownames(modules) <- rownames(x)
    colnames(modules) <- names(sigs)
    return(modules)
}

### HELPER FUNCTIONS ###

#' @importFrom SummarizedExperiment rowData
.rowdata2taxonomy <- function(x){
  
    taxa.prefix <- paste0(getTaxonomyRankPrefixes()[tolower(taxonomyRanks(x))], "__")
    full.tax <- rowData(x)[ , taxonomyRanks(x)]
    full.tax <- mapply(function(rank, prefix) paste0(prefix, rank), full.tax, taxa.prefix)
    full.tax <- apply(full.tax, 1L, function(row) paste(row, collapse = "|"))
    full.tax <- gsub("\\|\\w__NA", "", full.tax)
    
    return(full.tax)
}

.add_prefix_to_rowdata <- function(x){
  
    taxa.prefix <- paste0(getTaxonomyRankPrefixes()[tolower(taxonomyRanks(x))], "__")
    
    full.tax <- rowData(x)[ , taxonomyRanks(x)]
    
    full.tax <- mapply(function(rank, prefix)
      ifelse(is.na(rank), NA, paste0(prefix, rank)), full.tax, taxa.prefix)
    
    return(full.tax)
}