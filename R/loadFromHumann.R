#' Import HUMAnN results to \code{TreeSummarizedExperiment}
#
#' @param file a single \code{character} value defining the file
#' path of the HUMAnN file. The file must be in merged HUMAnN format.
#'
#' @param colData a DataFrame-like object that includes sample names in
#' rownames, or a single \code{character} value defining the file
#' path of the sample metadata file. The file must be in \code{tsv} format
#' (default: \code{colData = NULL}).
#' 
#' @param ... additional arguments:
#' \itemize{
#'   \item{\code{assay.type}:} {A single character value for naming 
#'   \code{\link[TreeSummarizedExperiment:TreeSummarizedExperiment-class]{assay}} 
#'   (default: \code{assay.type = "counts"})}
#'   \item{\code{removeTaxaPrefixes}:} {\code{TRUE} or \code{FALSE}: Should
#'     taxonomic prefixes be removed? (default:
#'     \code{removeTaxaPrefixes = FALSE})}
#' }
#'
#' @details
#' Import HUMAnN (version 3.0) results of pathways or gene families. Input must
#' be in merged HUMAnN format. (See the HUMAnN documentation and
#' \code{humann_join_tables} method.)
#' 
#' The function parses gene/pathway information along with taxonomy information
#' from the input file. This information is stored to \code{rowData}. Abundances
#' are stored to \code{assays}.
#' 
#' Usually the workflow includes also taxonomy data from Metaphlan. See
#' \link[=loadFromMetaphlan]{loadFromMetaphlan} to load the data to \code{TreeSE}.
#'
#' @return  A
#' \code{\link[TreeSummarizedExperiment:TreeSummarizedExperiment-class]{TreeSummarizedExperiment}}
#' object
#'
#' @name loadFromMetaphlan
#' @seealso
#' \code{\link[=loadFromMetaphlan]{loadFromMetaphlan}}
#' \code{\link[=makeTreeSEFromPhyloseq]{makeTreeSEFromPhyloseq}}
#' \code{\link[=makeTreeSEFromBiom]{makeTreeSEFromBiom}}
#' \code{\link[=makeTreeSEFromDADA2]{makeTreeSEFromDADA2}}
#' \code{\link[=loadFromQIIME2]{loadFromQIIME2}}
#' \code{\link[=loadFromMothur]{loadFromMothur}}
#'
#' @export
#' @author Leo Lahti and Tuomas Borman. Contact: \url{microbiome.github.io}
#' 
#' @references
#' Beghini F, McIver LJ, Blanco-Míguez A, Dubois L, Asnicar F, Maharjan S,
#' Mailyan A, Manghi P, Scholz M, Thomas AM, Valles-Colomer M, Weingart G,
#' Zhang Y, Zolfo M, Huttenhower C, Franzosa EA, & Segata N (2021)
#' Integrating taxonomic, functional, and strain-level profiling of diverse
#' microbial communities with bioBakery 3 eLife 10:e65088.
#' 
#' @examples
#' \dontrun{
#' 
#' tse <- loadFromHumann("pathabundance.tsv")
#' 
#' }
#' 
NULL

loadFromHumann <- function(file, colData = NULL, ...){
    ################################ Input check ###############################
    if(!mia:::.is_non_empty_string(file)){
        stop("'file' must be a single character value.",
             call. = FALSE)
    }
    if (!file.exists(file)) {
        stop(file, " does not exist", call. = FALSE)
    }
    if(!is.null(colData) &&
       !(.is_non_empty_string(colData) || is.data.frame(colData) ||
         is.matrix(colData) || is(colData, "DataFrame")) ){
        stop("'colData' must be a single character value, DataFrame or NULL.",
             call. = FALSE)
    }
    ############################## Input check end #############################
    # Read metaphlan data
    data <- .read_humann(file)
    # Create TreeSE from the data
    tse <- .create_tse_from_humann(data, ...)
    # Add colData if provided
    if( !is.null(colData) ){
        tse <- .add_coldata(tse, colData)
    }
    return(tse)
}

################################ HELP FUNCTIONS ################################

# Read Humann file, catch error if it occurs
.read_humann <- function(file){
    # Read the table. Catch error and give more informative message
    table <- tryCatch(
        {
            read.delim(file, check.names = FALSE)
        },
        error = function(condition){
            stop("Error while reading ", file,
                 "\nPlease check that the file is in merged Metaphlan file format.",
                 call. = FALSE)
        }
    )
    # In the first column name, there is "# " prefix. Remove it
    colnames(table)[1] <- gsub("# ", "", colnames(table)[1])
    # Replace spaces with underscore
    colnames(table)[1] <- gsub(" ", "_", colnames(table)[1])
    # Add rownames
    rownames(table) <- table[, 1] 
    # Check that file is in right format
    if( .check_humann(table) ){
        stop("Error while reading ", file,
             "\nPlease check that the file is in merged HUMAnN file format.",
             call. = FALSE)
    }
    return(table)
}

# Check that metaphlan file contains correct information
.check_humann <- function(data){
    # Get rowdata columns
    rowdata_col <- c("Pathway", "Gene_Family")
    rowdata_id <- unlist(lapply(rowdata_col, grep, colnames(data)))
    rowdata_columns <- data[ , rowdata_id, drop = FALSE]
    # Get columns that go to assay
    assay_columns <- data[ , -rowdata_id, drop = FALSE]
    # Initialize result 
    result <- TRUE
    
    # Check rowdata column names that they contain right information, and check that 
    # rest of the columns represents abundances in samples.
    # If these requirements are met, give FALSE. Otherwise, give TRUE.
    if( any(colnames(rowdata_columns) %in% c("Pathway", "Gene_Family")) && 
        is.numeric(unlist(assay_columns)) ){
        result <- FALSE
    }
    return(result)
}

# This function parses humann file and creates tse from it.
.create_tse_from_humann <- function(data, colData, assay.type = "counts", ...){
    # Get rowdata columns
    rowdata_col <- c("Pathway", "Gene_Family")
    rowdata_id <- unlist(lapply(rowdata_col, grep, colnames(data)))
    rowdata <- data[ , rowdata_id, drop = FALSE]
    # Get columns that go to assay
    assay <- data[ , -rowdata_id, drop = FALSE]
    
    # Parse rowdata. The data includes gene/pathway info along with taxonomy
    # info. Separate the information to own columns.
    rowdata_id <- unlist(lapply(rowdata_col, grep, colnames(rowdata)))
    # Get the column and split gene/pathway and taxonomy info
    rowdata_temp <- strsplit(rowdata[[1]], "\\|",)
    # Are some rows missing taxonomy info? Get their indices.
    missing_taxa <- lengths(rowdata_temp) == 1
    # Create a df from the list of splitted data.
    rowdata_temp <- data.frame(t(data.frame(rowdata_temp)))
    # Replace missing taxonomy info with NA
    rowdata_temp[missing_taxa, 2] <- NA
    # Replace column names
    colnames(rowdata_temp) <- c(colnames(rowdata), "Taxonomy")
    
    # Now we have rowdata that includes gene/pathway info in one column and
    # taxonomy info in other. Let's parse the taxonomy info so that species
    # genus etc levels are in unique columns.
    taxonomy <- mia:::.parse_taxonomy(
        rowdata_temp, column_name = "Taxonomy", sep = "\\.", ...)
    
    # Convert all data to DataFrame
    rowdata <- DataFrame(rowdata)
    rowdata_temp <- DataFrame(rowdata_temp)
    
    # Create rowdata from the information that is parsed
    colnames(rowdata) <- paste0(colnames(rowdata), "_long")
    rowdata <- cbind(rowdata, rowdata_temp[, 1, drop = FALSE])
    rowdata <- cbind(rowdata, taxonomy)
    
    # Create assays
    assay <- as.matrix(assay)
    assays <- SimpleList(counts = assay)
    names(assays) <- assay.type

    # Create TreeSE
    tse <- TreeSummarizedExperiment(assays = assays, rowData = rowdata)
    return(tse)
}
