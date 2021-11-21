#' Import Metaphlan results to \code{TreeSummarizedExperiment}
#
#' @param metaphlan a single \code{character} value defining the file
#'   path of the Metaphlan file. The file must be in merged Metaphlan format.
#'
#' @param sample_meta a single \code{character} value defining the file
#'   path of the sample metadata file. The file must be in \code{tsv} format
#'   (default: \code{sample_meta = NULL}).
#'   
#' @param phy_tree a single \code{character} value defining the file
#'   path of the phylogenetic tree.
#'   (default: \code{phy_tree = NULL}).
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
#' \code{\link[TreeSummarizedExperiment:TreeSummarizedExperiment-class]{TreeSummarizedExperiment}}
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
#' tse <- loadFromMetaphlan(file_path)
#' tse
#' }
#' 

loadFromMetaphlan <- function(metaphlan, sample_meta = NULL, phy_tree = NULL, ...){
    ################################ Input check ################################
    if(!.is_non_empty_string(metaphlan)){
        stop("'metaphlan' must be a single character value.",
             call. = FALSE)
    }
    if (!file.exists(metaphlan)) {
        stop(metaphlan, " does not exist", call. = FALSE)
    }
    if(!is.null(sample_meta) && !.is_non_empty_string(sample_meta)){
        stop("'sample_meta' must be a single character value or NULL.",
             call. = FALSE)
    }
    if(!is.null(phy_tree) && !.is_non_empty_string(phy_tree)){
        stop("'phy_tree' must be a single character value or NULL.",
             call. = FALSE)
    }
    ############################## Input check end #############################
    # Parse assay and rowdata from the file
    assay_and_rowdata <- .parse_assay_and_rowdata_from_metaphlan(metaphlan, ...)
    assay <- assay_and_rowdata$assay
    rowdata <- assay_and_rowdata$rowdata
    
    # Load sample meta data if it is provided, otherwise initialize empty table
    if (!is.null(sample_meta)) {
        coldata <- read.table(file = sample_meta, header = TRUE, sep = "\t")
    } else {
        coldata <- S4Vectors:::make_zero_col_DataFrame(ncol(assay))
        rownames(coldata) <- colnames(assay)
    }
    
    # Load tree if it is provided
    if (!is.null(phy_tree)) {
        tree <- ape::read.tree(phy_tree)
    } else {
        tree <- NULL
    }
    
    # Check that features and samples match
    assay <- .set_feature_tab_dimnames(assay, coldata, rowdata)
    # Create TSE
    x <- TreeSummarizedExperiment(assays = list(counts = assay), 
                                    rowData = rowdata, 
                                    colData = coldata, 
                                    rowTree = tree)
    return(x)
}

################################ HELP FUNCTIONS ################################
# Parse assay and rowdata from metaphlan file 
.parse_assay_and_rowdata_from_metaphlan <- function(file, ...){
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
    # Lowest level includes rownames
    rownames(rowdata) <- taxonomy[ , ncol(taxonomy)]
    rownames(assay) <- taxonomy[ , ncol(taxonomy)]
    
    return(list(assay = assay, rowdata = rowdata))
}

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
