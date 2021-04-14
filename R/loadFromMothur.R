#' Import Mothur results to \code{SummarizedExperiment}
#' 
#' This method creates a \code{SummarizedExperiment} object from \code{Mothur}
#' files that are provided as an input. 
#'
#' @param sharedFile a single \code{character} value defining the file
#'   path of the feature table to be imported. 
#'   The File has to be in \code{shared} file format.
#'
#' @param taxonomyFile a single \code{character} value defining the file
#'   path of the taxonomy table to be imported. 
#'   The File has to be in \code{cons.taxonomy} file format. (default: \code{taxonomyFile = NULL}).
#'
#' @param designFile a single \code{character} value defining the file path
#'   of the sample metadata to be imported. 
#'   The File has to be in \code{desing} file format. (default: \code{designFile = NULL}).
#'
#' @param ... additional arguments:
#'
#' @details
#' Results exported from Mothur can be imported as a
#' \code{SummarizedExperiment} using \code{loadFromMothur}. Except for the
#' \code{sharedFile}, the other data types, \code{taxonomyFile}, and
#' \code{designFile}, are optional, but are highly -------------------------Check this
#' encouraged to be provided.
#'
#' @return  A
#' \code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#' object
#'
#' @name loadFromMothur
#' @seealso
#' \code{\link[=makeTreeSummarizedExperimentFromphyloseq]{makeTreeSummarizedExperimentFromphyloseq}}
#' \code{\link[=makeTreeSummarizedExperimentFromBiom]{makeTreeSummarizedExperimentFromBiom}}
#' \code{\link[=makeTreeSummarizedExperimentFromDADA2]{makeTreeSummarizedExperimentFromDADA2}}
#' \code{\link[=loadFromQIIME2}]{loadFromQIIME2}}
#'
#' @export
#' @author Leo Lahti and Tuomas Borman. Contact: \url{microbiome.github.io}
#'
#' @references
#' \url{https://mothur.org/}
#' 
#' @examples
#' # Abundance table
#' counts <- "~/Desktop/Baxter_FITs_Microbiome_2016_fit.final.tx.1.subsample.shared"
#' # Taxa table
#' taxa <- "~/Desktop/Baxter_FITs_Microbiome_2016_fit.final.tx.1.cons.taxonomy"
#' # Sample meta data
#' meta <- "~/Desktop/Baxter_FITs_Microbiome_2016_mapping.csv"
#' 
#' # Creates tse object from files
#' \dontrun{tse <- loadFromMothur(counts, taxa, meta)}
#' 

loadFromMothur <- function(sharedFile,
                           taxonomyFile = NULL,
                           designFile = NULL,
                           ...) {
  
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
    feature_tab_and_data_to_colData <- .read_mothur_feature(sharedFile, ...)
    # Extracts feature_tab
    feature_tab <- feature_tab_and_data_to_colData$assay
    # Extracts data that goes to colData
    data_to_colData <- feature_tab_and_data_to_colData$colData
    
    # If rowData information exists, gets that. Otherwise, tax_tab is just data frame without information
    if (!is.null(taxonomyFile)) {
        taxa_tab <- .read_mothur_taxonomy(taxonomyFile, feature_tab, ...)
    } else {
        taxa_tab <- S4Vectors:::make_zero_col_DataFrame(nrow(feature_tab))
    }
    
    # If colData informationor data_to_colData exists, gets that. Otherwise, sample_tab is just data frame without information
    if (!is.null(designFile) %% !is.null(data_to_colData)) {
        sample_meta <- .read_mothur_sample_meta(designFile, data_to_colData)
    } else {
        sample_meta <- S4Vectors:::make_zero_col_DataFrame(ncol(feature_tab))
    }
    
    return(SummarizedExperiment(
        assays = S4Vectors::SimpleList(counts = feature_tab),
        rowData = taxa_tab,
        colData = sample_meta
    ))
    
}
# These extra information must be added to colData. Return list of assay and extra info
.read_mothur_feature <- function(sharedFile, ...){
  
    if (.get_mothur_file_type(sharedFile) != "shared") {
      stop("The input '", sharedFile, "' must be in `shared` format.",
           call. = FALSE)
    }
    # 
    # # Reads the file
    # assay <- read.table(sharedFile, check.names=FALSE, header=TRUE,
    #                     sep="\t", stringsAsFactors=FALSE)
    # 
    # # File contains additional columns "label", "numOtus", and "Group". They are removed
    # assay$label <- NULL
    # assay$numOtus <- NULL
    # # Information of "Group" column is saved, it includes sample information
    # sample_names <- assay$Group
    # assay$Group <- NULL
    # 
    # # Initializes rownames
    # rownames(assay) <- NULL
    # # Saves sample information to rownames of table
    # rownames(assay) <- sample_names
    # 
    # # Transposes and converts assay from dataframe to matrix
    # assay <- t(as.matrix(assay))
    # return(assay)
  
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

.read_mothur_taxonomy <- function(taxonomyFile, feature_tab, ...){
  
    if (.get_mothur_taxonomy_file_type(taxonomyFile) != "cons.taxonomy") {
      stop("The input '", taxonomyFile, "' must be in `cons.taxonomy` format.",
           call. = FALSE)
    }
    
    # # Reads the file
    # rowData <- read.table(taxonomyFile, check.names=FALSE,
    #                       header=TRUE, sep="\t", stringsAsFactors=FALSE)
    # # Deletes additional information between taxa levels.
    # rowData$Taxonomy <- gsub("[\"]", "", rowData$Taxonomy)
    # rowData$Taxonomy <- gsub("[(1-100)]", "", rowData$Taxonomy)
    # # Separate taxa to own columns
    # rowData <- tidyr::separate(rowData, 
    #                            "Taxonomy", 
    #                            into=c("Kingdom", "Phylum", "Order", "Class", 
    #                                   "Family", "Genus"), sep=";", 
    #                            extra="merge")
    # 
    # # Deletes ";" from the end of Genus names
    # rowData$Genus <- gsub(";", "", rowData$Genus)
    # # Deletes addition "Size" column
    # rowData$Size <- NULL
    # 
    # # Adds taxa to rownames
    # rownames(rowData) <- rowData$OTU
    # # Deletes addition "OTU" column
    # rowData$OTU <- NULL
    # 
    # # # Creates a matrix from the table
    # # rowData <- as.matrix(rowData)
    # # rowData <- (rowData)
    # 
    # # Gets only those taxa that are included in feature_tab
    # #rowData <- rowData[rownames(feature_tab),] # delete this <---------------------------------------------
    # return(rowData)
    
    # Saves column that includes taxonomical information
    MOTHUR_TAX_COL <- "Taxonomy"
    # Reads the file
    data <- read.table(taxonomyFile, check.names=FALSE,
                       header=TRUE, sep="\t", stringsAsFactors=FALSE)
    # Checks that colnames contain information, it is not NULL, and taxonomical information is present 
    if ( !(length(colnames(data)) > 0) || is.null(colnames(data)) || is.null(data[[MOTHUR_TAX_COL]]) ){
      stop("'taxonomy' does not include taxonomical information.",
           call. = FALSE)
    }
    
    # Removes additional characters between taxa
    data [,MOTHUR_TAX_COL] <- gsub("[\"]", "", data [,MOTHUR_TAX_COL])
    data [,MOTHUR_TAX_COL] <- gsub("[(1-100)]", "", data [,MOTHUR_TAX_COL])
    # Splits taxa level into separate columns
    tax <- tidyr::separate(data , MOTHUR_TAX_COL, into=c("Kingdom", "Phylum", "Order", "Class", "Family", "Genus"), sep=";", extra="merge")
    # Removes ";" from the end of genus level names
    tax$Genus <- gsub(";", "", tax$Genus)
    
    # Combines the extracted taxonomical data and unprocessed data and creates a data frame
    # rowData <- cbind(DataFrame(tax),
    #                  DataFrame(data[, !(colnames(data) %in% MOTHUR_TAX_COL)]))
    rowData <- tax
    # If the data includes "OTU" column, that includes taxa
    if( !is.null(rowData$OTU) ){
      # Checks if the rownames are the same as in assay.
      if( !identical(rowData$OTU, rownames(feature_tab)) ){
        stop("taxa in 'taxonomyFile' does not match to taxa in 'sharedFile'.",
             call. = FALSE)
      }
      # Adds rownames
      rownames(rowData) <- rowData$OTU
    }
    return(rowData)

}

.read_mothur_sample_meta <- function(designFile, data_to_colData, ...){
  
    if (.get_mothur_file_type(designFile) != "design") {
      stop("The input '", designFile, "' must be in `design` format.",
           call. = FALSE)
    }
  
    # Reads the file
    colData <- read.table(designFile, check.names=FALSE,
                          header=TRUE, sep="\t", stringsAsFactors=FALSE)
    
    # Combines the extracted colData and data from the assay
    colData <- cbind(colData, data_to_colData)
    
    # If the data includes 'group' column that includes the names of the samples
    if( !is.null(colData$group) ){
      # Checks if the sample names extracted from sample meta data ('group' column) 
      # are the same as the ones that were extracted from assay ('Group' column)
      if( !identical(colData$group, colData$Group) ){
        stop("sample names in 'designFile' does not match to sample names in 'sharedFile'.",
             call. = FALSE)
      }
      # Adds sample names to rownames of the data 
      rownames(colData) <- colData$group
    }
    return(colData)
  
}


#' extract file extension
#' @noRd
.get_mothur_file_type <- function(file) {
    ex <- strsplit(basename(file), split = ".", fixed = TRUE)[[1]]
    ex[length(ex)]
}

#' extract file extension
#' @noRd
.get_mothur_taxonomy_file_type <- function(file) {
    ex <- strsplit(basename(file), split = ".", fixed = TRUE)[[1]]
    paste(c(ex[length(ex)-1], ".", ex[length(ex)]), collapse= "")
}
