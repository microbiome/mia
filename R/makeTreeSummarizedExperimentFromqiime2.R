#' Import qiime2 results to `TreeSummarizedExperiment`
#'
#' Import qiime2 sample metadata, feature table, taxonomy table, phylogenetic
#' tree, and representative sequence for each feature to
#' `TreeSummarizedExperiment` object.
#'
#' @param featureTableFile character, file path of the feature table.
#' @param taxonomyTableFile character, file path of the taxonomy table.
#'   default `NULL`.
#' @param sampleMetaFile character, file path of the sample metadata in tsv
#'  format, default `NULL`.
#' @param featureNamesAsRepseq logical, only worked while argument `repSeqFile`
#'   is `NULL`. The default is `TRUE`, which means if the name of feature is a
#'   DNA sequence, use it as representative sequence of the feature.
#' @param repSeqFile character, file path of the representative sequence for
#'   each feature, default `NULL`.
#' @param phyTreeFile character, file path of the phylogenetic tree, default
#'   `NULL`.
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
#' repSeqFile <- system.file("extdata", "refseq.qza", package = "mia")
#' tse <- makeTreeSummarizedExperimentFromqiime2(
#'   featureTableFile = featureTableFile,
#'   taxonomyTableFile = taxonomyTableFile,
#'   sampleMetaFile = sampleMetaFile,
#'   repSeqFile = repSeqFile,
#'   phyTreeFile = phyTreeFile
#' )
#'
#' tse
makeTreeSummarizedExperimentFromqiime2 <- function(featureTableFile,
                                                   taxonomyTableFile = NULL,
                                                   sampleMetaFile = NULL,
                                                   featureNamesAsRepseq = TRUE,
                                                   repSeqFile = NULL,
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

    # check the row.names of feature table: DNA sequence or not
    # if is DNA, row.names(feature_tab) is set as refseq
    if  (.is_dna_seq(row.names(feature_tab))) {
        refseq <- row.names(feature_tab)
        refseq_nm <- paste0("OTU", seq_along(refseq))
        names(refseq) <- refseq_nm
    } else {
        refseq <- NULL
    }

    if (!is.null(repSeqFile)){
        if (!is.null(refseq)) {
            warning(
                "use sequences from `repSeqFile` as representative sequences",
                "argument `featureNamesAsRepseq` does not work"
            )
        }
        refseq <- .read_qza(repSeqFile)
    }


    if(!is.null(refseq)) {
        refseq <- DNAStringSet(refseq)
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
#' @return [`phyloseq::otu_table-class`] object for feature table,
#'   [`phyloseq::taxonomyTable-class`] object for taxonomic table,
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
        stop("The input `file` must be in `qza` format (qiime2 Artifact)")
    }

    unzipped_file <- unzip(file, exdir = temp)
    meta_file <- grep("metadata.yaml", unzipped_file, value = TRUE)
    metadata <- read_yaml(meta_file[1])
    format <- metadata$format
    uuid <- metadata$uuid

    format_files <- c(
        "feature-table.biom", "taxonomy.tsv",
        "tree.nwk", "dna-sequences.fasta"
    )
    formats <- c(
        "BIOMV210DirFmt", "TSVTaxonomyDirectoryFormat",
        "NewickDirectoryFormat", "DNASequencesDirectoryFormat"
    )
    file <- file.path(temp, uuid, "data", format_files[match(format, formats)])

    res <- switch (
        format,
        BIOMV210DirFmt = .read_q2biom(file),
        TSVTaxonomyDirectoryFormat = .read_q2taxa(file),
        NewickDirectoryFormat = read.tree(file),
        DNASequencesDirectoryFormat = readDNAStringSet(file),
        stop(
            "Only files in format of 'BIOMV210DirFmt', ",
            "'TSVTaxonomyDirectoryFormat', NewickDirectoryFormat' and ",
            "'DNASequencesDirectoryFormat' are supported."
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
    as.matrix(taxa)
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
#' @param trim_rank_prefix logical whether remove leading characters from
#'   taxonomic levels, e.g. k__ or D_0__, default `FALSE`.
#' @return [`phyloseq::taxonomyTable-class`] object.
#' @keywords internal
#' @importFrom rlang .data
#' @importFrom tidyr separate
#' @noRd
.parse_q2taxonomy <- function(taxa, sep = "; |;", trim_rank_prefix = FALSE) {
    stopifnot(is.logical(trim_rank_prefix))
    taxa <- data.frame(taxa)
    if(trim_rank_prefix){
        # remove leading characters from GG
        taxa$Taxon <- gsub("[kpcofgs]__", "", taxa$Taxon)
        #remove leading characters from SILVA
        taxa$Taxon <- gsub("D_\\d__", "", taxa$Taxon)
    }

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
    rownames(taxa) <- taxa[, "Feature.ID"]
    taxa <- taxa[, setdiff(colnames(taxa), "Feature.ID")]

    taxa
}

#' check the row.names of feature table is DNA sequence or not
#' @keywords internal
#' @importFrom Biostrings DNAStringSet
#' @noRd
.is_dna_seq <- function(seq){
    seq <- try({DNAStringSet(seq)}, silent = TRUE)
    !inherits(seq, "try-error")
}

#' extract file extension
#' @noRd
.get_ext <- function(file) {
    ex <- strsplit(basename(file), split = ".", fixed = TRUE)[[1]]
    ex[length(ex)]
}
