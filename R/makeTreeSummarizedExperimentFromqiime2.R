#' Import qiime2 results to `TreeSummarizedExperiment`
#'
#' Results exported from Qiime2 can be imported as a `TreeSummarizedExperiment`
#' using `makeTreeSummarizedExperimentFromqiime2`. Except for the
#' `featureTableFile`, the other data types, `taxonomyTableFile`, `refSeqFile`
#' and `phyTreeFile`, are optional, but are highly encouraged to be provided.
#'
#' @param featureTableFile character, file path of the feature table.
#' @param taxonomyTableFile character, file path of the taxonomy table.
#'   default `NULL`.
#' @param sampleMetaFile character, file path of the sample metadata in tsv
#'   format, default `NULL`.
#' @param featureNamesAsRefseq logical, only worked while argument `refSeqFile`
#'   is `NULL`. The default is `TRUE`, which means if the elments of `row.names`
#'   of feature table is DNA sequences, use it as reference sequences.
#' @param refSeqFile character, file path of the representative sequence for
#'   each feature, default `NULL`.
#' @param phyTreeFile character, file path of the phylogenetic tree, default
#'   `NULL`.
#' @details
#' Both arguments `featureNamesAsRefseq` and `refSeqFile` were used to define
#' reference sequences of features. Firstly, `featureNamesAsRefseq` only worked
#' while `refSeqFile` is `NULL` (if elements of `row.names` of feature table are
#' DNA sequences,  set these sequences as reference sequences, or reference
#' sequences is `NULL`). Otherwise, the reference sequences were created from
#' `refSeqFile` if it is not `NULL`. Finally, no reference sequences is created
#' if `featureNameAsRefSeq` is `FALSE` and  `refSeqFile` is `NULL`.
#'
#' @return  An object of class
#'   [`TreeSummarizedExperiment::TreeSummarizedExperiment-class`]
#' @importFrom Biostrings DNAStringSet
#' @export
#' @author Yang Cao
#' @references \url{https://qiime2.org}
#' @examples
#' featureTableFile <- system.file("extdata", "table.qza", package = "mia")
#' taxonomyTableFile <- system.file("extdata", "taxonomy.qza", package = "mia")
#' sampleMetaFile <- system.file("extdata", "sample-metadata.tsv", package = "mia")
#' phyTreeFile <- system.file("extdata", "tree.qza", package = "mia")
#' refSeqFile <- system.file("extdata", "refseq.qza", package = "mia")
#' tse <- makeTreeSummarizedExperimentFromqiime2(
#'   featureTableFile = featureTableFile,
#'   taxonomyTableFile = taxonomyTableFile,
#'   sampleMetaFile = sampleMetaFile,
#'   refSeqFile = refSeqFile,
#'   phyTreeFile = phyTreeFile
#' )
#'
#' tse
makeTreeSummarizedExperimentFromqiime2 <- function(featureTableFile,
                                                   taxonomyTableFile = NULL,
                                                   sampleMetaFile = NULL,
                                                   featureNamesAsRefseq = TRUE,
                                                   refSeqFile = NULL,
                                                   phyTreeFile = NULL) {
    feature_tab <- .read_qza(featureTableFile)

    if (!is.null(taxonomyTableFile)) {
        taxa_tab <- .read_qza(taxonomyTableFile)
        taxa_tab <- .subset_taxa_in_feature(taxa_tab, feature_tab)
        taxa_tab <- .parse_q2taxonomy(taxa_tab)
    } else {
        taxa_tab <- S4Vectors:::make_zero_col_DataFrame(nrow(feature_tab))
    }

    if (!is.null(sampleMetaFile)) {
        sample_meta <- .read_q2sample_meta(sampleMetaFile)
    } else {
        sample_meta <- S4Vectors:::make_zero_col_DataFrame(ncol(feature_tab))
    }

    if (!is.null(phyTreeFile)) {
        tree <- .read_qza(phyTreeFile)
    } else {
        tree <- NULL
    }

    # if row.names(feature_tab) is a DNA sequence,  set it as refseq
    if (!is.null(refSeqFile)){
        refseq <- .read_qza(refSeqFile)
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


#' Read the qza file output from qiime2
#'
#' Import the qiime2 artifacts to R.
#'
#' @param file character, path of the input qza file. Only files in format of
#'   `BIOMV210DirFmt` (feature table), `TSVTaxonomyDirectoryFormat` (taxonomic
#'   table), `NewickDirectoryFormat` (phylogenetic tree ) and
#'   `DNASequencesDirectoryFormat` (representative sequences) are supported
#'    right now.
#' @param temp character, a temporary directory in which the qza file will be
#'   decompressed to, default `tempdir()`.
#' @return `matirx` object for feature table, `data.frame` for taxonomic table,
#'   [`ape::phylo`] object for phylogenetic tree,
#'   [`Biostrings::DNAStringSet-class`] for representative sequences of taxa.
#' @noRd
#' @importFrom yaml read_yaml
#' @importFrom utils unzip
#' @importFrom ape read.tree
#' @importFrom Biostrings readDNAStringSet
#' @keywords internal
.read_qza <- function(file, temp = tempdir()) {
    if (!file.exists(file)) {
        stop(file, " does not exist", call. = FALSE)
    }
    if (.get_ext(file) != "qza") {
        stop("The input '", file, "' must be in `qza` format (qiime2 Artifact)", call. = FALSE)
    }

    unzipped_file <- unzip(file, exdir = temp)
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

#' Read qiime2

#' Read qiime2 feature table
#'
#' Import qiime2 feature table into otu_table object
#'
#' @param file character, file name of the biom file.
#' @return [`phyloseq::otu_table-class`] object.
#' @noRd
.read_q2biom <- function(file) {
    .require_package("biomformat")
    biomobj <- suppressWarnings(biomformat::read_biom(file))
    feature_tab <- as(biomformat::biom_data(biomobj),"matrix")

    feature_tab
}

#' Read qiime2 taxa file
#' @keywords internal
#' @importFrom utils read.table
#' @noRd
.read_q2taxa <- function(file) {
    taxa <- utils::read.table(file, sep = '\t', header = TRUE)

    taxa
}

#' Read qiime2 sample meta data file
#' @keywords internal
#' @importFrom utils read.table
#' @noRd
.read_q2sample_meta <- function(file) {
    sam <- read.table(file = file, header = TRUE, sep = "\t", comment.char = "")
    rownames(sam) <- as.character(sam[, 1])

    # remove row: #q2:types
    idx <- which(sam == "#q2:types", arr.ind = TRUE)
    sam <- sam[-idx[, "row"], ]

    S4Vectors::DataFrame(sam)
}


#' Subset taxa according to the taxa in feature table
#' @keywords internal
#' @noRd
.subset_taxa_in_feature <- function(taxa_tab, feature_tab) {
    idx <- match(row.names(feature_tab), taxa_tab[, "Feature.ID"])
    taxa_tab <- taxa_tab[idx, , drop = FALSE]
    row.names(taxa_tab) <- taxa_tab[, "Feature.ID"]

    taxa_tab
}

#' Parse qiime2 taxa in different taxonomic levels
#' @param taxa `tax_table` object.
#' @param sep character string containing a regular expression, separator
#'  between different taxonomic levels, defaults to on compatible with both
#'  GreenGenes and SILVA `; |;"`.
#' @return  a `data.frame`.
#' @keywords internal
#' @importFrom rlang .data
#' @importFrom tidyr separate
#' @noRd
.parse_q2taxonomy <- function(taxa, sep = "; |;") {
    taxa <- data.frame(taxa)

    #  work with any combination of taxonomic ranks available
    all_ranks <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
    all_prefixes <- c("k__", "p__", "c__", "o__", "f__", "g__", "s__")
    present_ranks <- strsplit(taxa[, "Taxon"], sep)
    # at certain rank (e.g. genus): some taxa can not be determined which genus
    # it assigned to
    idx <- which.max(lengths(present_ranks))
    present_ranks <- present_ranks[[idx]]
    present_prefixes <- substr(present_ranks, 1, 3)

    taxa <- suppressWarnings(
        separate(
            taxa, .data$Taxon,
            all_ranks[match(present_prefixes, all_prefixes)],
            sep = sep
        )
    )
    taxa <- sapply(taxa, function(x) ifelse(x == "", NA_character_, x))
    taxa <- data.frame(taxa)
    # make sure confidence in numeric
    if ("Confidence" %in% colnames(taxa)) {
        taxa[["Confidence"]] <- as.numeric(taxa[["Confidence"]])
    }
    if("Feature.ID" %in% names(taxa)) {
        rownames(taxa) <- taxa[["Feature.ID"]]
        taxa["Feature.ID"] <- NULL
    }

    taxa
}

#' check the row.names of feature table is DNA sequence or not
#' @keywords internal
#' @importFrom Biostrings DNAStringSet
#' @noRd
.rownames_as_dna_seq <- function(seq){
    names(seq) <- paste0("seq_", seq_along(seq))
    seq <- try({DNAStringSet(seq)}, silent = TRUE)
    if (inherits(seq, "try-error")) {
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
