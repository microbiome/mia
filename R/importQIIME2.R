#' Import QIIME2 results to \code{TreeSummarizedExperiment}
#'
#' Results exported from QIMME2 can be imported as a
#' \code{TreeSummarizedExperiment} using \code{importQIIME2}. Except for the
#' \code{featureTableFile}, the other data types, \code{taxonomyTableFile},
#' \code{refSeqFile} and \code{phyTreeFile}, are optional, but are highly
#' encouraged to be provided.
#'
#' @param featureTableFile a single \code{character} value defining the file
#'   path of the feature table to be imported.
#'
#' @param taxonomyTableFile a single \code{character} value defining the file
#'   path of the taxonomy table to be imported. (default:
#'   \code{taxonomyTableFile = NULL}).
#'
#' @param sampleMetaFile a single \code{character} value defining the file path
#'   of the sample metadata to be imported. The file has to be in tsv format.
#'   (default: \code{sampleMetaFile = NULL}).
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
#' \itemize{
#'   \item{\code{temp}:} {the temporary directory used for decompressing the
#'     data. (default: \code{tempdir()})}
#'   \item{\code{removeTaxaPrefixes}:} {\code{TRUE} or \code{FALSE}: Should
#'     taxonomic prefixes be removed? (default:
#'     \code{removeTaxaPrefixes = FALSE})}
#' }
#'
#' @details
#' Both arguments \code{featureNamesAsRefSeq} and \code{refSeqFile} can be used
#' to define reference sequences of features. \code{featureNamesAsRefSeq} is
#' only taken into account, if \code{refSeqFile} is \code{NULL}. No reference
#' sequences are tried to be created, if \code{featureNameAsRefSeq} is
#' \code{FALSE} and \code{refSeqFile} is \code{NULL}.
#'
#' @return  A
#' \code{\link[TreeSummarizedExperiment:TreeSummarizedExperiment-class]{TreeSummarizedExperiment}}
#' object
#'
#' @name importQIIME2
#' @seealso
#' \code{\link[=makeTreeSEFromPhyloseq]{makeTreeSEFromPhyloseq}}
#' \code{\link[=makeTreeSEFromBiom]{makeTreeSEFromBiom}}
#' \code{\link[=makeTreeSEFromDADA2]{makeTreeSEFromDADA2}}
#' \code{\link[=loadFromMothur]{loadFromMothur}}
#'
#' @export
#' @author Yang Cao
#'
#' @references
#' Bolyen E et al. 2019: Reproducible, interactive, scalable and extensible
#' microbiome data science using QIIME 2. Nature Biotechnology 37: 852–857.
#' \url{https://doi.org/10.1038/s41587-019-0209-9}
#'
#' \url{https://qiime2.org}
#'
#' @examples
#' featureTableFile <- system.file("extdata", "table.qza", package = "mia")
#' taxonomyTableFile <- system.file("extdata", "taxonomy.qza", package = "mia")
#' sampleMetaFile <- system.file("extdata", "sample-metadata.tsv", package = "mia")
#' phyTreeFile <- system.file("extdata", "tree.qza", package = "mia")
#' refSeqFile <- system.file("extdata", "refseq.qza", package = "mia")
#' tse <- importQIIME2(
#'   featureTableFile = featureTableFile,
#'   taxonomyTableFile = taxonomyTableFile,
#'   sampleMetaFile = sampleMetaFile,
#'   refSeqFile = refSeqFile,
#'   phyTreeFile = phyTreeFile
#' )
#'
#' tse

#' @importFrom S4Vectors make_zero_col_DFrame
importQIIME2 <- function(featureTableFile,
                            taxonomyTableFile = NULL,
                            sampleMetaFile = NULL,
                            featureNamesAsRefSeq = TRUE,
                            refSeqFile = NULL,
                            phyTreeFile = NULL,
                           ...) {
    .require_package("yaml")
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

    feature_tab <- importQZA(featureTableFile, ...)

    if (!is.null(taxonomyTableFile)) {
        taxa_tab <- importQZA(taxonomyTableFile, ...)
        taxa_tab <- .subset_taxa_in_feature(taxa_tab, feature_tab)
    } else {
        taxa_tab <- S4Vectors::make_zero_col_DFrame(nrow(feature_tab))
        rownames(taxa_tab) <- rownames(feature_tab)
    }

    if (!is.null(sampleMetaFile)) {
        sample_meta <- .read_q2sample_meta(sampleMetaFile)
    } else {
        sample_meta <- S4Vectors::make_zero_col_DFrame(ncol(feature_tab))
        rownames(sample_meta) <- colnames(feature_tab)
    }

    if (!is.null(phyTreeFile)) {
        tree <- importQZA(phyTreeFile, ...)
    } else {
        tree <- NULL
    }

    # if row.names(feature_tab) is a DNA sequence,  set it as refseq
    if (!is.null(refSeqFile)){
        refseq <- importQZA(refSeqFile, ...)
    } else if (featureNamesAsRefSeq) {
        refseq <- .rownames_as_dna_seq(rownames(feature_tab))
    } else {
        refseq <- NULL
    }
    
    feature_tab <- .set_feature_tab_dimnames(feature_tab, sample_meta, taxa_tab)
    TreeSummarizedExperiment(
        assays = S4Vectors::SimpleList(counts = feature_tab),
        rowData = taxa_tab,
        colData = sample_meta,
        rowTree = tree,
        referenceSeq = refseq
    )
}

#' Read the qza file output from QIIME2
#'
#' Import the QIIME2 artifacts to R.
#'
#' @param file character, path of the input qza file. Only files in format of
#'   `BIOMV210DirFmt` (feature table), `TSVTaxonomyDirectoryFormat` (taxonomic
#'   table), `NewickDirectoryFormat` (phylogenetic tree ) and
#'   `DNASequencesDirectoryFormat` (representative sequences) are supported
#'    right now.
#' @param temp character, a temporary directory in which the qza file will be
#'   decompressed to, default `tempdir()`.
#' @return `matrix` object for feature table, `DataFrame` for taxonomic table,
#'   [`ape::phylo`] object for phylogenetic tree,
#'   [`Biostrings::DNAStringSet-class`] for representative sequences of taxa.
#'   
#' @name importQIIME2
#' @export
#'
#' @examples 
#' # Read individual files
#' featureTableFile <- system.file("extdata", "table.qza", package = "mia")
#' taxonomyTableFile <- system.file("extdata", "taxonomy.qza", package = "mia")
#' sampleMetaFile <- system.file("extdata", "sample-metadata.tsv", package = "mia")
#' 
#' assay <- importQZA(featureTableFile)
#' rowdata <- importQZA(taxonomyTableFile, removeTaxaPrefixes = TRUE)
#' coldata <- read.table(sampleMetaFile, header = TRUE, sep = "\t", comment.char = "")
#' 
#' # Assign rownames 
#' rownames(coldata) <- coldata[, 1]
#' coldata[, 1] <- NULL
#' 
#' # Order coldata based on assay
#' coldata <- coldata[match(colnames(assay), rownames(coldata)), ]
#' 
#' # Create SE from individual files
#' se <- SummarizedExperiment(assays = list(assay), rowData = rowdata, colData = coldata)
#' se
#' 
#' @importFrom utils unzip
#' @importFrom ape read.tree
#' @importFrom Biostrings readDNAStringSet
importQZA <- function(file, temp = tempdir(), ...) {
    if (!file.exists(file)) {
        stop(file, " does not exist", call. = FALSE)
    }
    if (.get_ext(file) != "qza") {
        stop("The input '", file, "' must be in `qza` format (QIIME2 Artifact)",
            call. = FALSE)
    }

    unzipped_file <- unzip(file, exdir = temp)
    on.exit(unlink(c(unzipped_file,unique(dirname(unzipped_file))),
                    recursive = TRUE))
    meta_file <- grep("metadata.yaml", unzipped_file, value = TRUE)
    metadata <- yaml::read_yaml(meta_file[1])
    uuid <- metadata$uuid

    format <- metadata$format
    # support for multiple BIOM formats:  V100, V210
    if (grepl("BIOMV", format)) {
        format <- "BIOMV"
    }

    format_files <- c(
        "feature-table.biom", "taxonomy.tsv",
        "tree.nwk", "dna-sequences.fasta"
    )
    formats <- c(
        "BIOMV", "TSVTaxonomyDirectoryFormat",
        "NewickDirectoryFormat", "DNASequencesDirectoryFormat"
    )
    file <- file.path(temp, uuid, "data", format_files[match(format, formats)])

    res <- switch (
        format,
        BIOMV = .read_q2biom(file),
        TSVTaxonomyDirectoryFormat = .read_q2taxa(file, ...),
        NewickDirectoryFormat = read.tree(file),
        DNASequencesDirectoryFormat = readDNAStringSet(file),
        stop(
            "Only files in format of 'BIOMV210DirFmt', ",
            "'TSVTaxonomyDirectoryFormat', NewickDirectoryFormat' and ",
            "'DNASequencesDirectoryFormat' are supported.",
            call. = FALSE
        )
    )
    res
}

#' Read QIIME2 feature table
#' @param file character, file name of the biom file.
#' @noRd
.read_q2biom <- function(file) {
    .require_package("biomformat")
    biomobj <- suppressWarnings(biomformat::read_biom(file))
    feature_tab <- as(biomformat::biom_data(biomobj),"matrix")

    feature_tab
}

#' Read QIIME2 taxa file
#' @keywords internal
#' @importFrom utils read.table
#' @noRd
.read_q2taxa <- function(file, ...) {
    taxa_tab <- utils::read.table(file, sep = '\t', header = TRUE)
    
    confidence <- NULL
    featureID <- NULL
    # make sure confidence in numeric
    if ("Confidence" %in% colnames(taxa_tab)) {
        confidence <- as.numeric(taxa_tab[,"Confidence"])
    }
    if("Feature.ID" %in% names(taxa_tab)) {
        featureID <- taxa_tab[,"Feature.ID"]
    }
    
    taxa_tab <- .parse_taxonomy(taxa_tab, sep = "; |;", column_name = "Taxon", ...)
    
    rownames(taxa_tab) <- featureID
    taxa_tab$Confidence <- confidence
    
    taxa_tab
}

#' Read QIIME2 sample meta data file
#' @keywords internal
#' @importFrom utils read.table
#' @noRd
.read_q2sample_meta <- function(file) {
    sam <- read.table(file = file, header = TRUE, sep = "\t", comment.char = "")
    rownames(sam) <- as.character(sam[, 1])

    # Find if there is #q2:types row, and store its index
    idx <- which(sam == "#q2:types", arr.ind = TRUE)

    # If the length is over zero, "#q2:types" row was found. Then it is removed.
    if(!(length(idx)==0)){
        sam <- sam[-idx[, "row"], ]
    }

    S4Vectors::DataFrame(sam)
}


#' Subset taxa according to the taxa in feature table
#' @keywords internal
#' @noRd
.subset_taxa_in_feature <- function(taxa_tab, feature_tab) {
    idx <- match( rownames(feature_tab), rownames(taxa_tab) )
    taxa_tab <- taxa_tab[idx, , drop = FALSE]

    taxa_tab
}

#' check the row.names of feature table is DNA sequence or not
#' @keywords internal
#' @importFrom Biostrings DNAStringSet
#' @noRd
.rownames_as_dna_seq <- function(seq){
    names(seq) <- paste0("seq_", seq_along(seq))
    seq <- try({DNAStringSet(seq)}, silent = TRUE)
    if (is(seq, "try-error")) {
        return(NULL)
    }

    seq
}

#' extract file extension
#' @noRd
.get_ext <- function(file) {
    ex <- strsplit(basename(file), split = ".", fixed = TRUE)[[1]]
    ex[length(ex)]
}
