#' Import Metaphlan results to \code{SummarizedExperiment}
#
#' @param file a single \code{character} value defining the file
#'   path of the file. The file must be in merged Metaphlan format.
#'
#' @param ... additional arguments:
#' \itemize{
#'   \item{\code{removeTaxaPrefixes}:} {\code{TRUE} or \code{FALSE}: Should
#'     taxonomic prefixes be removed? (default:
#'     \code{removeTaxaPrefixes = FALSE})}
#' }
#'
#' @details
#' Import Metaphlan results. Input must be in merged Metaphlan format.
#'
#' @return  A
#' \code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#' object
#'
#' @name loadFromMetaphlan
#' @seealso
#' \code{\link[=makeTreeSummarizedExperimentFromPhyloseq]{makeTreeSummarizedExperimentFromPhyloseq}}
#' \code{\link[=makeSummarizedExperimentFromBiom]{makeSummarizedExperimentFromBiom}}
#' \code{\link[=makeTreeSummarizedExperimentFromDADA2]{makeTreeSummarizedExperimentFromDADA2}}
#' \code{\link[=loadFromQIIME2]{loadFromQIIME2}}
#' \code{\link[=loadFromMothur]{loadFromMothur}}
#'
#' @export
#' @author Leo Lahti and Tuomas Borman. Contact: \url{microbiome.github.io}
#' 
#' @references
#'
#' @examples
#' \dontrun{
#' # File path
#' file_path <- "/data/metaphlan_result.txt"
#' # Import data
#' se <- loadFromMetaphlan(file)
#' se
#' }
#' 

loadFromMetaphlan <- function(file, ...){
    ################################ Input check ################################
    if(!.is_non_empty_string(file)){
        stop("'file' must be a single character value.",
             call. = FALSE)
    }
    if (!file.exists(file)) {
        stop(file, " does not exist", call. = FALSE)
    }
    ############################## Input check end #############################
    # Read the table. Catch error and give more informative message
    table <- tryCatch(
        {
            read.table(file, header = TRUE, comment.char = "#")
        },
            error = function(condition){
                stop("Error while reading ", file,
                     "\nPlease check that the file is in merged Metaphlan file format.",
                     call. = FALSE)
        }
    )
    # Subset so that only those rows are included that include all taxonomic levels
    table <- .get_rows_that_include_lowest_level(table)
    # Get those columns that belong to rowData
    rowdata <- table[, 1:2, drop = FALSE]
    # Get those columns that belong to assay
    assay <- table[, 3:ncol(table), drop = FALSE]
    # Parse taxonomic levels
    taxonomy <- .parse_taxonomy(rowdata[ , 1, drop = FALSE], sep = "\\|", column_name = "clade_name", ...)
    # Add parsed taxonomy level information to rowdata
    rowdata <- cbind(taxonomy, rowdata)
    
    # Create SE
    x <- SummarizedExperiment::SummarizedExperiment(assays = list(counts = assay), 
                                                    rowData = rowdata)
    # Add rownames from the lowest taxonomic level
    rownames(x) <- .get_taxonomic_label(x)
    return(x)
}

################################ HELP FUNCTIONS ################################
# Get the lowest level of the string that contains multiple taxonomic levels with prefixes
# Output is single character that specifies the rank, e.g, "s" == "Species"
.get_lowest_taxonomic_level <- function(string){
    # Get indices that specify location of rank prefixes 
    levels <- gregexpr("([kpcofgs]+)__", string)[[1]]
    # Get the location of lowest rank
    lowest_level_ind <- levels[length(levels)]
    # Get the lowest rank that was found
    lowest_level <- substr(string, start = lowest_level_ind, stop = lowest_level_ind)
    return(lowest_level)
}

# Subset rows so that only those rows are included that have all the taxonomic levels
# that are present in the data
.get_rows_that_include_lowest_level <- function(table){
    # Get the lowest level of each row
    levels <- sapply(table[["clade_name"]], FUN = .get_lowest_taxonomic_level)
    # Order the data and get the lowest level of that the data includes
    order <- c("s", "g", "f", "o", "c", "p", "k")
    lowest_level_found <- levels[order(match(levels, order))][1]
    # Get those rows that include information at lowest level
    table <- table[grepl(paste0(lowest_level_found, "__"), table[["clade_name"]]), ]
    return(table)
}
