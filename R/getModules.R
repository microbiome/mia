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
    function(x, sigs, ...){
        modules <- getModules(x, sigs, ...)
        rowData(x) <- cbind(rowData(x), modules)
        return(x)
    }
)

#' @rdname getModules
#' @export
setMethod("getModules", signature = c(x = "SummarizedExperiment"),
    function(x, sigs, ...){
        modules <- .make_modules_table(x, sigs, ...)
        return(modules)
    }
)

# Define function to construct modules table based on bugsigdb signatures
.make_modules_table <- function(x, sigs, ...){
    # Retrieve prefixes for taxonomic ranks
    tax.ranks <- getTaxonomyRankPrefixes()
    # Retrieve taxonomic labels for features
    tax.labs <- x |>
        getTaxonomyLabels(make.unique = FALSE, with.rank = TRUE) |>
        tolower()
    # Initialise empty list for signatures
    sig.list <- list()
    # For every signature in output from bugsigdbr::getSignatures
    for( i in seq_along(sigs) ){
        # Extract deepest taxonomic rank
        sig <- gsub(".*\\|", "", sigs[[i]])
        # Get letter for taxonomic rank
        sig.labs <- substr(sig, 1, 1)
        # Match letter to taxonomic rank
        sig.ranks <- names(tax.ranks)[match(sig.labs, tax.ranks)]
        # Replace letter with taxonomic rank in signature name
        sig.labs <- tolower(paste(sig.ranks, gsub(".*\\w__", "", sig), sep = ":"))
        # Find which features belong to the current signature
        sig.list[[i]] <- tax.labs %in% sig.labs
    }
    # Build modules table from signature list
    modules <- do.call(cbind, sig.list)
    # Add names to rows and columns
    rownames(modules) <- rownames(x)
    colnames(modules) <- names(sigs)
    return(modules)
}
