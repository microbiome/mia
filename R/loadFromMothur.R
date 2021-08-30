#' Import Mothur results as a \code{SummarizedExperiment}
#' 
#' This method creates a \code{SummarizedExperiment} object from \code{Mothur}
#' files provided as input. 
#'
#' @param sharedFile a single \code{character} value defining the file
#'   path of the feature table to be imported. The File has to be in 
#'   \code{shared file} format as defined in Mothur documentation.
#'
#' @param taxonomyFile a single \code{character} value defining the file path of
#'   the taxonomy table to be imported. The File has to be in \code{taxonomy
#'   file} or \code{constaxonomy file} format  as defined in Mothur
#'   documentation. (default: \code{taxonomyFile = NULL}).
#'
#' @param designFile a single \code{character} value defining the file path of
#'   the sample metadata to be imported. The File has to be in \code{desing
#'   file} format as defined in Mothur documentation. (default: \code{designFile
#'   = NULL}).
#'
#' @details
#' Results exported from Mothur can be imported as a
#' \code{SummarizedExperiment} using \code{loadFromMothur}. Except for the
#' \code{sharedFile}, the other data types, \code{taxonomyFile}, and
#' \code{designFile}, are optional, but are highly encouraged to be provided.
#'
#' @return  A
#' \code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#' object
#'
#' @name loadFromMothur
#' @seealso
#' \code{\link[=makeTreeSummarizedExperimentFromPhyloseq]{makeTreeSummarizedExperimentFromPhyloseq}}
#' \code{\link[=makeSummarizedExperimentFromBiom]{makeSummarizedExperimentFromBiom}}
#' \code{\link[=makeTreeSummarizedExperimentFromDADA2]{makeTreeSummarizedExperimentFromDADA2}}
#' \code{\link[=loadFromQIIME2]{loadFromQIIME2}}
#'
#' @author Leo Lahti and Tuomas Borman. Contact: \url{microbiome.github.io}
#'
#' @references
#' \url{https://mothur.org/}
#' \url{https://mothur.org/wiki/shared_file/}
#' \url{https://mothur.org/wiki/taxonomy_file/}
#' \url{https://mothur.org/wiki/constaxonomy_file/}
#' \url{https://mothur.org/wiki/design_file/}
#' 
#' @examples
#' # Abundance table
#' counts <- system.file("extdata", "mothur_example.shared", package = "mia")
#' # Taxa table (in "cons.taxonomy" or "taxonomy" format)
#' taxa <- system.file("extdata", "mothur_example.cons.taxonomy", package = "mia")
#' #taxa <- system.file("extdata", "mothur_example.taxonomy", package = "mia")
#' # Sample meta data
#' meta <- system.file("extdata", "mothur_example.design", package = "mia")
#' 
#' # Creates se object from files
#' se <- loadFromMothur(counts, taxa, meta)
#' se
NULL

#' @rdname loadFromMothur
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @export
loadFromMothur <- function(sharedFile,
                           taxonomyFile = NULL,
                           designFile = NULL) {

    # input check
    if(!.is_non_empty_string(sharedFile)){
        stop("'sharedFile' must be a single character value.", 
            call. = FALSE)
    }
    if(!is.null(taxonomyFile) && !.is_non_empty_string(taxonomyFile)){
        stop("'taxonomyFile' must be a single character value or NULL.", 
            call. = FALSE)
    }
    if(!is.null(designFile) && !.is_non_empty_string(designFile)){
        stop("'designFile' must be a single character value or NULL.",
            call. = FALSE)
    }
    
    # Reads the sharedFile 
    feature_tab_and_data_to_colData <- .read_mothur_feature(sharedFile)
    # Extracts feature_tab
    feature_tab <- feature_tab_and_data_to_colData$assay
    # Extracts data that goes to colData
    data_to_colData <- feature_tab_and_data_to_colData$colData
    
    # If rowData information exists, gets that. Otherwise, tax_tab is just 
    # data.frame without information
    if (!is.null(taxonomyFile)) {
        taxa_tab <- .read_mothur_taxonomy(taxonomyFile, feature_tab)
    } else {
        taxa_tab <- S4Vectors:::make_zero_col_DataFrame(nrow(feature_tab))
        rownames(taxa_tab) <- rownames(feature_tab)
    }
    
    # If colData informationor data_to_colData exists, gets that. Otherwise, 
    # sample_tab is just data frame without information
    if (!is.null(designFile) && !is.null(data_to_colData)) {
        sample_meta <- .read_mothur_sample_meta(designFile, data_to_colData)
    } else {
        sample_meta <- S4Vectors:::make_zero_col_DataFrame(ncol(feature_tab))
        rownames(sample_meta) <- colnames(feature_tab)
    }
    
    feature_tab <- .set_feature_tab_dimnames(feature_tab, sample_meta, taxa_tab)
    SummarizedExperiment(assays = S4Vectors::SimpleList(counts = feature_tab),
                         rowData = taxa_tab,
                         colData = sample_meta)
}

# These extra information must be added to colData. Return list of assay and 
# extra info
.read_mothur_feature <- function(sharedFile){
  
    if (!.is_mothur_shared_file(sharedFile)) {
        stop("The input '", sharedFile, "' must be in `shared` format.",
             call. = FALSE)
    }
  
    # Stores name of columns will be included in colData not in assays
    MOTHUR_NON_ASSAY_COLS <- c("label","numOtus","Group")
    data <- read.table(sharedFile, check.names=FALSE, header=TRUE,
                       sep="\t", stringsAsFactors=FALSE)
    # Checks that colnames contain information and it is not NULL
    if ( !(length(colnames(data)) > 0) || is.null(colnames(data)) ){
        stop("'shared' does not include names of taxa.",
           call. = FALSE)
    }
    # Takes all columns but not those that goes to colData, 
    # and transforms the data frame to matrix
    assay <- as.matrix(data[,!(colnames(data) %in% MOTHUR_NON_ASSAY_COLS )])
    # Transposes the matrix --> taxa to rows
    assay <- t(assay)
    # Gets those data that goes to colData, and creates a data frame it
    colData <- DataFrame(data[,MOTHUR_NON_ASSAY_COLS])
    return(list(assay = assay,
                colData = colData))
}

.read_mothur_taxonomy <- function(taxonomyFile, feature_tab){
    
    # If the file is in "cons.taxonomy" format
    if (.is_mothur_constaxonomy_file(taxonomyFile, feature_tab)) {
        data <- read.table(taxonomyFile, check.names=FALSE,
                           header=TRUE, sep="\t", stringsAsFactors=FALSE)
    } 
    # If the file is in "taxonomy" format, adds column names
    else if (.is_mothur_taxonomy_file(taxonomyFile, feature_tab)){
        data <- read.table(taxonomyFile, check.names=FALSE,
                           header=FALSE, sep="\t", 
                           stringsAsFactors=FALSE, 
                           col.names = c("OTU", "Taxonomy"))
    }
    # Else the file is not either gives an error
    else{
        stop("The input '", taxonomyFile, "' must be provided in the ",
             "`taxonomy` or `cons.taxonomy` format. In addition, it must ",
             "match the data of the 'sharedFile'",
             call. = FALSE)
    }
  
    # Column that includes taxonomical information
    MOTHUR_TAX_COL <- "Taxonomy"

    # Checks that colnames contain information, it is not NULL, and taxonomical 
    # information is present 
    if ( !(length(colnames(data)) > 0) || 
         is.null(colnames(data)) || 
         is.null(data[[MOTHUR_TAX_COL]]) ){
        stop("'taxonomy' does not include taxonomical information.",
             call. = FALSE)
    }
    
    # Removes additional characters between taxa
    data [,MOTHUR_TAX_COL] <- gsub("[\"]", "", data [,MOTHUR_TAX_COL])
    data [,MOTHUR_TAX_COL] <- gsub("[(1-100)]", "", data [,MOTHUR_TAX_COL])
    # Splits taxa level into separate columns
    into <- c("Kingdom", "Phylum", "Order", "Class", "Family", "Genus")
    tax <- tidyr::separate(data,
                           MOTHUR_TAX_COL,
                           into=into,
                           sep=";",
                           extra="merge")
    # Removes ";" from the end of genus level names
    tax$Genus <- gsub(";", "", tax$Genus)
    rowData <- tax
    # Adds rownames
    rownames(rowData) <- rowData$OTU
    return(rowData)
}

.read_mothur_sample_meta <- function(designFile, data_to_colData){
    # Checks if file is in "design" format. data_to_colData$Group includes 
    # sample names that were extracted from assay, i.e. sharedFile
    if (!.is_mothur_design_file(designFile, data_to_colData$Group)) {
      stop("The input '", designFile, "' must be in `design` format, 
           and it must inlude same sample names as 'sharedFile'.",
           call. = FALSE)
    }
  
    # Reads the file
    colData <- read.table(designFile, check.names=FALSE,
                          header=TRUE, sep="\t",
                          stringsAsFactors=FALSE)
    
    # Combines the extracted colData and data from the assay
    colData <- cbind(colData, data_to_colData)
    
    # Adds sample names to rownames of the data 
    rownames(colData) <- colData$group
    return(colData)
}

#' extract file extension
#' @noRd
.is_mothur_shared_file <- function(file){
    result <- FALSE
    # Columns that every shared file should have
    columns_that_must_be_found <- c("label","Group", "numOtus")
    
    # Reads the data
    data <- read.table(file, check.names=FALSE, header=TRUE,
                       sep="\t", stringsAsFactors=FALSE)
    
    # If data contains column names, then it is shared file
    if( identical(colnames(data)[1:3], columns_that_must_be_found) ){
        result <- TRUE
    }
    return(result)
}

#' extract file extension
#' @noRd
.is_mothur_taxonomy_file <- function(file, feature_tab){
    result <- FALSE
    
    # Reads the data
    data <- read.table(file, check.names=FALSE, header=FALSE,
                       sep="\t", stringsAsFactors=FALSE)
    
    # If data contains 2 columns and the first column includes same taxa as 
    # feature_tab in rownames, then it is taxonomy file
    if( ncol(data) == 2 && identical(data[,1], rownames(feature_tab)) ){
        result <- TRUE
    }
    return(result)
}

#' extract file extension
#' @noRd
.is_mothur_constaxonomy_file <- function(file, feature_tab){
    result <- FALSE
    # Columns that every constaxonomy file should have
    columns_that_must_be_found <- c("OTU", "Size", "Taxonomy")
    
    # Reads the data
    data <- read.table(file, check.names=FALSE, header=TRUE,
                       sep="\t", stringsAsFactors=FALSE)
    
    # If data contains column names, and "OTU" column that includes same taxa as 
    # feature_tab, 
    # then it is constaxonomy file
    if( identical(colnames(data), columns_that_must_be_found) && 
        identical(data$OTU, rownames(feature_tab)) ){
        result <- TRUE
    }
    return(result)
}

.is_mothur_design_file <- function(file, sample_names){
    result <- FALSE
    # Reads the data
    data <- read.table(file, check.names=FALSE, header=TRUE,
                       sep="\t", stringsAsFactors=FALSE)
    
    # If data contains "group" column that include sample names, then it is 
    # design file
    if( identical(data$group, sample_names) ){
        result <- TRUE
    }
    return(result)
}
