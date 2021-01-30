#' Import QIIME2 results to `TreeSummarizedExperiment`
#'
#' Results exported from QIMME2 can be imported as a `TreeSummarizedExperiment`
#' using `loadFromQIIME2`. Except for the
#' `featureTableFile`, the other data types, `taxonomyTableFile`, `refSeqFile`
#' and `phyTreeFile`, are optional, but are highly encouraged to be provided.
#'
#' @param featureTableFile a single \code{character} value defining the file
#'   path of the feature table to be imported.
#' @param taxonomyTableFile a single \code{character} value defining the file
#'   path of the taxonomy table to be imported. (default: `NULL`).
#' @param sampleMetaFile a single \code{character} value defining the file
#'   path of the sample metadata to be imported. The file has to be in tsv
#'   format. (default: `NULL`).
#' @param featureNamesAsRefseq \code{TRUE} or \code{FALSE}: Should the feature
#'   names of the feature table be regarded as reference sequences? This setting
#'   will be disregarded, if \code{refSeqFile} is not \code{NULL}. If the
#'   feature names do not contain valid DNA characters only, the reference
#'   sequences will not be set.
#' @param refSeqFile a single \code{character} value defining the file
#'   path of the reference sequences for each feature. (default: `NULL`).
#' @param phyTreeFile a single \code{character} value defining the file
#'   path of the phylogenetic tree. (default: `NULL`).
#' @param ... additional arguments:
#' \itemize{
#'   \item{\code{temp}:} {the temporary directory used for decompressing the
#'     data. (default: \code{tempdir()})}
#'   \item{\code{removeTaxaPrefixes}:} {\code{TRUE} or \code{FALSE}: Should
#'     taxonomic prefixes be removed? (default: \code{FALSE}).)}
#' }
#'
#' @details
#' Both arguments `featureNamesAsRefseq` and `refSeqFile` can be used to define
#' reference sequences of features. `featureNamesAsRefseq` is only taken into
#' account, if `refSeqFile` is `NULL`. No reference sequences are tried to be
#' created, if `featureNameAsRefSeq` is `FALSE` and  `refSeqFile` is `NULL`.
#'
#' @return  An object of class
#'   [`TreeSummarizedExperiment::TreeSummarizedExperiment-class`]
#'
#' @name loadFromQIIME2
#' @seealso
#' \code{\link[=makeTreeSummarizedExperimentFromphyloseq]{makeTreeSummarizedExperimentFromphyloseq}}
#' \code{\link[=makeTreeSummarizedExperimentFromBiom]{makeTreeSummarizedExperimentFromBiom}}
#' \code{\link[=makeTreeSummarizedExperimentFromDADA2]{makeTreeSummarizedExperimentFromDADA2}}
#'
#' @export
#' @author Yang Cao
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
#' tse <- loadFromQIIME2(
#'   featureTableFile = featureTableFile,
#'   taxonomyTableFile = taxonomyTableFile,
#'   sampleMetaFile = sampleMetaFile,
#'   refSeqFile = refSeqFile,
#'   phyTreeFile = phyTreeFile
#' )
#'
#' tse
loadFromQIIME2 <- function(featureTableFile,
                           taxonomyTableFile = NULL,
                           sampleMetaFile = NULL,
                           featureNamesAsRefseq = TRUE,
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
    if(!.is_a_bool(featureNamesAsRefseq)){
        stop("'featureNamesAsRefseq' must be TRUE or FALSE.", call. = FALSE)
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

    feature_tab <- .read_qza(featureTableFile, ...)

    if (!is.null(taxonomyTableFile)) {
        taxa_tab <- .read_qza(taxonomyTableFile, ...)
        taxa_tab <- .subset_taxa_in_feature(taxa_tab, feature_tab)
        taxa_tab <- .parse_q2taxonomy(taxa_tab, ...)
    } else {
        taxa_tab <- S4Vectors:::make_zero_col_DataFrame(nrow(feature_tab))
    }

    if (!is.null(sampleMetaFile)) {
        sample_meta <- .read_q2sample_meta(sampleMetaFile)
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
    } else if (featureNamesAsRefseq) {
        refseq <- .rownames_as_dna_seq(rownames(feature_tab))
    } else {
        refseq <- NULL
    }

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
#' @noRd
#' @importFrom yaml read_yaml
#' @importFrom utils unzip
#' @importFrom ape read.tree
#' @importFrom Biostrings readDNAStringSet
#' @keywords internal
.read_qza <- function(file, temp = tempdir(), ...) {
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
    metadata <- read_yaml(meta_file[1])
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
        TSVTaxonomyDirectoryFormat = .read_q2taxa(file),
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
.read_q2taxa <- function(file) {
    taxa_tab <- utils::read.table(file, sep = '\t', header = TRUE)

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
    idx <- match(rownames(feature_tab), taxa_tab[, "Feature.ID"])
    taxa_tab <- taxa_tab[idx, , drop = FALSE]
    rownames(taxa_tab) <- taxa_tab[, "Feature.ID"]

    taxa_tab
}

#' Parse QIIME2 taxa in different taxonomic levels
#' @param taxa_tab `data.frame` object.
#' @param sep character string containing a regular expression, separator
#'  between different taxonomic levels, defaults to one compatible with both
#'  GreenGenes and SILVA `; |;"`.
#' @return  a `data.frame`.
#' @keywords internal
#' @importFrom IRanges CharacterList IntegerList
#' @importFrom S4Vectors DataFrame
#' @noRd
.parse_q2taxonomy <- function(taxa_tab, sep = "; |;",
                              removeTaxaPrefixes = FALSE, ...) {
    if(!.is_a_bool(removeTaxaPrefixes)){
        stop("'removeTaxaPrefixes' must be TRUE or FALSE.", call. = FALSE)
    }

    confidence <- NULL
    featureID <- NULL
    # make sure confidence in numeric
    if ("Confidence" %in% colnames(taxa_tab)) {
        confidence <- as.numeric(taxa_tab[,"Confidence"])
    }
    if("Feature.ID" %in% names(taxa_tab)) {
        featureID <- taxa_tab[,"Feature.ID"]
    }

    #  work with any combination of taxonomic ranks available
    all_ranks <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
    all_prefixes <- c("k__", "p__", "c__", "o__", "f__", "g__", "s__")

    # split the taxa strings
    taxa_split <- CharacterList(strsplit(taxa_tab[,"Taxon"],sep))
    # extract present prefixes
    taxa_prefixes <- lapply(taxa_split, substr, 1L, 3L)
    # match them to the order given by present_prefixes
    taxa_prefixes_match <- lapply(taxa_prefixes, match, x = all_prefixes)
    taxa_prefixes_match <- IntegerList(taxa_prefixes_match)
    # get the taxa values
    if(removeTaxaPrefixes){
        taxa_split <- lapply(taxa_split,
                             gsub,
                             pattern = "([kpcofgs]+)__",
                             replacement = "")
        taxa_split <- CharacterList(taxa_split)
    }
    # extract by order matches
    taxa_split <- taxa_split[taxa_prefixes_match]
    #
    if(length(unique(lengths(taxa_split))) != 1L){
        stop("Internal error. Something went wrong")
    }
    taxa_tab <- DataFrame(as.matrix(taxa_split))
    colnames(taxa_tab) <- all_ranks
    rownames(taxa_tab) <- featureID
    taxa_tab$Confidence <- confidence

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
