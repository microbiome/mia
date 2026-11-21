#' Split taxonomy strings into taxonomic ranks
#'
#' Some pipelines report the taxonomy of each feature as a single string, such
#' as \code{"k__Bacteria; p__Firmicutes; c__Clostridia"}.
#' \code{parseTaxonomy} splits these strings into one column per taxonomic
#' rank. The same parser is used by the importers of mia, e.g.
#' \code{\link[=importQIIME2]{importQIIME2()}} and
#' \code{\link[=importMetaPhlAn]{importMetaPhlAn()}}.
#'
#' @param x \code{Character vector} of taxonomy strings, or a
#' \code{data.frame} or \code{DataFrame} (e.g. \code{rowData(tse)}) with
#' the taxonomy strings in column \code{col.name}.
#'
#' @param col.name \code{Character scalar}. The column of \code{x} with the
#' taxonomy strings. Ignored if \code{x} is a character vector.
#' (Default: \code{"Taxon"})
#'
#' @param sep \code{Character scalar}. A regular expression that separates the
#' ranks. The default fits Greengenes and QIIME 2 (\code{"; "} or \code{";"}).
#' Use, e.g., \code{"\\\\|"} for MetaPhlAn. (Default: \code{"; |;"})
#'
#' @param remove.prefix \code{Logical scalar}. Should the rank prefixes (e.g.
#' \code{"p__"}) be removed from the values? (Default: \code{FALSE})
#'
#' @details
#' The ranks are recognised by their prefixes: \code{d__} (domain),
#' \code{sk__} (superkingdom), \code{k__} (kingdom), \code{p__} (phylum),
#' \code{c__} (class), \code{o__} (order), \code{f__} (family),
#' \code{g__} (genus), \code{s__} (species) and \code{t__} (strain).
#' A rank that is missing from a string gets \code{NA}, and ranks that are
#' missing from all strings are left out.
#'
#' @return A \code{DataFrame} with one column per rank, named in lower case
#' (\code{"kingdom"}, \code{"phylum"}, ...) as in the importers. The row
#' names are the names of \code{x}, or the row names if \code{x} is a
#' table.
#'
#' @seealso
#' \code{\link[=taxonomy-methods]{setTaxonomyRanks}},
#' \code{\link[=taxonomy-methods]{getTaxonomyLabels}}
#'
#' @name parseTaxonomy
#' @export
#'
#' @examples
#' taxa <- c(
#'     ASV1 = "k__Bacteria; p__Firmicutes; c__Bacilli; g__Lactobacillus",
#'     ASV2 = "k__Bacteria; p__Bacteroidetes; c__Bacteroidia")
#' parseTaxonomy(taxa)
#' parseTaxonomy(taxa, remove.prefix = TRUE)
#'
#' # MetaPhlAn-style strings
#' parseTaxonomy(
#'     "k__Bacteria|p__Firmicutes|c__Clostridia", sep = "\\|",
#'     remove.prefix = TRUE)
#'
#' # Taxonomy stored in one column of rowData
#' library(SummarizedExperiment)
#' tse <- SummarizedExperiment(
#'     assays = list(counts = matrix(1:4, 2L, dimnames = list(names(taxa),
#'         c("S1", "S2")))),
#'     rowData = DataFrame(Taxon = taxa))
#' rowData(tse) <- cbind(
#'     parseTaxonomy(rowData(tse), remove.prefix = TRUE), rowData(tse))
#' rowData(tse)
#'
NULL

parseTaxonomy <- function(
        x, col.name = "Taxon", sep = "; |;", remove.prefix = FALSE){
    if( is.character(x) ){
        x <- DataFrame(Taxon = unname(x), row.names = names(x))
        col.name <- "Taxon"
    }
    if( !(is.data.frame(x) || is(x, "DataFrame")) ){
        stop("'x' must be a character vector, data.frame or DataFrame.",
            call. = FALSE)
    }
    if( !(.is_non_empty_string(col.name) && col.name %in% colnames(x)) ){
        stop("'col.name' must be a column of 'x'.", call. = FALSE)
    }
    taxa <- as.character(x[[col.name]])
    if( anyNA(taxa) ){
        stop("Taxonomy strings must not contain NA values.", call. = FALSE)
    }
    res <- .parse_taxonomy(
        data.frame(Taxon = taxa), col.name = "Taxon", sep = sep,
        remove.prefix = remove.prefix)
    rownames(res) <- rownames(x)
    return(res)
}
