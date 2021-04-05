#' Import Mothur results to \code{TreeSummarizedExperiment}
#'
#' @param featureTableFile a single \code{character} value defining the file
#'   path of the feature table to be imported.
#'
#' @param taxonomyTableFile a single \code{character} value defining the file
#'   path of the taxonomy table to be imported. (default:
#'   \code{taxonomyTableFile = NULL}).
#'
#' @param sampleMetaFile a single \code{character} value defining the file path
#'   of the sample metadata to be imported. (default: \code{sampleMetaFile = NULL}).
#'
#' @param featureNamesAsRefSeq \code{TRUE} or \code{FALSE}: Should the feature
#'   names of the feature table be regarded as reference sequences? This setting
#'   will be disregarded, if \code{refSeqFile} is not \code{NULL}. If the
#'   feature names do not contain valid DNA characters only, the reference
#'   sequences will not be set.
#'
#' @param refSeqFile a single \code{character} value defining the file path of
#'   the reference sequences for each feature. (default: \code{refSeqFile =
#'   NULL}).
#'
#' @param phyTreeFile a single \code{character} value defining the file path of
#'   the phylogenetic tree. (default: \code{phyTreeFile = NULL}).
#'
#' @param ... additional arguments:
#'
#' @details
#' Results exported from Mothur can be imported as a
#' \code{TreeSummarizedExperiment} using \code{loadFromMothur}. Except for the
#' \code{featureTableFile}, the other data types, \code{taxonomyTableFile},
#' \code{refSeqFile} and \code{phyTreeFile}, are optional, but are highly
#' encouraged to be provided.
#'
#' @return  A
#' \code{\link[TreeSummarizedExperiment:TreeSummarizedExperiment-class]{TreeSummarizedExperiment}}
#' object
#'
#' @name loadFromMothur
#' @seealso
#' \code{\link[=makeTreeSummarizedExperimentFromphyloseq]{makeTreeSummarizedExperimentFromphyloseq}}
#' \code{\link[=makeTreeSummarizedExperimentFromBiom]{makeTreeSummarizedExperimentFromBiom}}
#' \code{\link[=makeTreeSummarizedExperimentFromDADA2]{makeTreeSummarizedExperimentFromDADA2}}
#' \code{\link[=makeTreeSummarizedExperimentFromDADA2]{loadFromQIIME2}}
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

loadFromMothur <- function(featureTableFile,
                           taxonomyTableFile = NULL,
                           sampleMetaFile = NULL,
                           featureNamesAsRefSeq = TRUE,
                           refSeqFile = NULL,
                           phyTreeFile = NULL,
                           ...) {
  
    # input check
    if(!.is_non_empty_string(featureTableFile)){
      stop("'featureTableFile' must be a single character value.",
           call. = FALSE)
    }
    if(!is.null(taxonomyTableFile) && !.is_non_empty_string(taxonomyTableFile)){
      stop("'taxonomyTableFile' must be a single character value or NULL.",
           call. = FALSE)
    }
    if(!is.null(sampleMetaFile) && !.is_non_empty_string(sampleMetaFile)){
      stop("'sampleMetaFile' must be a single character value or NULL.",
           call. = FALSE)
    }
    if(!.is_a_bool(featureNamesAsRefSeq)){
      stop("'featureNamesAsRefSeq' must be TRUE or FALSE.", call. = FALSE)
    }
    if(!is.null(refSeqFile) && !.is_non_empty_string(refSeqFile)){
      stop("'refSeqFile' must be a single character value or NULL.",
           call. = FALSE)
    }
    if(!is.null(phyTreeFile) && !.is_non_empty_string(phyTreeFile)){
      stop("'phyTreeFile' must be a single character value or NULL.",
           call. = FALSE)
    }
    #
    
    # Reads the featureTablefile 
    feature_tab <- .read_mothur_feature(featureTableFile, ...)
    
    # If rowData information exists, gets that. Otherwise, tax_tab is just data frame without information
    if (!is.null(taxonomyTableFile)) {
      taxa_tab <- .read_mothur_taxonomy(taxonomyTableFile, feature_tab, ...)
    } else {
      taxa_tab <- S4Vectors:::make_zero_col_DataFrame(nrow(feature_tab))
    }
    
    # If colData information exists, gets that. Otherwise, sample_tab is just data frame without information
    if (!is.null(sampleMetaFile)) {
      sample_meta <- .read_mothur_sample_meta(sampleMetaFile, feature_tab)
    } else {
      sample_meta <- S4Vectors:::make_zero_col_DataFrame(ncol(feature_tab))
    }
    
    if (!is.null(phyTreeFile)) {
      tree <- .read_qza(phyTreeFile, ...)
    } else {
      tree <- NULL
    }
    
    # if row.names(feature_tab) is a DNA sequence,  set it as refseq
    if (!is.null(refSeqFile)){
      refseq <- .read_qza(refSeqFile, ...)
    } else if (featureNamesAsRefSeq) {
      refseq <- .rownames_as_dna_seq(rownames(feature_tab))
    } else {
      refseq <- NULL
    }
    
    return(TreeSummarizedExperiment(
      assays = S4Vectors::SimpleList(counts = feature_tab),
      rowData = taxa_tab,
      colData = sample_meta,
      rowTree = tree,
      referenceSeq = refseq
    ))
    
}

.read_mothur_feature <- function(featureTableFile, ...){

    # Reads the file
    assay <- read.table(featureTableFile, check.names=FALSE, header=TRUE,
                      sep="\t", stringsAsFactors=FALSE)
  
    # File contains additional columns "label", "numOtus", and "Group". They are removed
    assay$label <- NULL
    assay$numOtus <- NULL
    # Information of "Group" column is saved, it includes taxa information
    x <- assay$Group
    assay$Group <- NULL
    
    # Initializes rownames
    rownames(assay) <- NULL
    # Saves taxa to rownames
    rownames(assay) <- x
    
    # Transposes and converts assay from dataframe to matrix
    assay <- t(as.matrix(assay))
    
    return(assay)
    
}

.read_mothur_taxonomy <- function(taxonomyTableFile, feature_tab, ...){
    
    # Reads the file
    rowData <- read.table(taxonomyTableFile, check.names=FALSE,
                          header=TRUE, sep="\t", stringsAsFactors=FALSE)
    # Deletes additional information between taxa levels.
    rowData$Taxonomy <- gsub("[\"]", "", rowData$Taxonomy)
    rowData$Taxonomy <- gsub("[(1-100)]", "", rowData$Taxonomy)
    # Separate taxa to own columns
    rowData <- tidyr::separate(rowData, "Taxonomy", into=c("Kingdom",
                                                           "Phylum", "Order", "Class", "Family", "Genus"), sep=";",
                               extra="merge")
    # Deletes ";" from the end of Genus names
    rowData$Genus <- gsub(";", "", rowData$Genus)
    # Deletes addition "Size" column
    rowData$Size <- NULL
    
    # Adds taxa to rownames
    rownames(rowData) <- rowData$OTU
    # Creates a matrix from the table
    rowData <- as.matrix(rowData)
    
    # Gets only those taxa that are included in feature_tab
    rowData <- rowData[rownames(feature_tab),]
    
    return(rowData)

}

.read_mothur_sample_meta <- function(sampleMetaFile, feature_tab, ...){
  
    # Reads the file
    colData <- read.csv(sampleMetaFile, row.names=1, check.names=FALSE)
    # Gets only those samples that are included in feature_tab
    colData <- colData[colnames(feature_tab), ]
    
    return(colData)

  
}
