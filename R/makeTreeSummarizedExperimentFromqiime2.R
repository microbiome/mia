# This function are modified from `microbiomeMarker::import_qiime2`
# https://github.com/yiluheihei/microbiomeMarker/blob/master/R/import_qiime2.R

#' Import qiime2 results to `TreeSummairzedExperiment`
#'
#' Import [qiime2](https://qiime2.org) results to `TreeSummairzedExperiment`
#'
#' @param otu_qza character, file path of the feature table from qiime2.
#' @param taxa_qza character, file path of the taxonomic table from qiime2,
#'   default `NULL`.
#' @param sam_tab character, file path of the sample metadata in tsv format,
#'   default `NULL`.
#' @param refseq_qza character, file path of the representative sequences from
#'   qiime2, default `NULL`.
#' @param tree_qza character, file path of the phylogenetic tree from
#'   qiime2, default `NULL`.
#' @return  An object of class
#'   [`TreeSummarizedExperiment::TreeSummarizedExperiment-class`]
#' @export
#' @examples
#' otuqza_file <- system.file("extdata", "table.qza", package = "mia")
#' taxaqza_file <- system.file("extdata", "taxonomy.qza", package = "mia")
#' sample_file <- system.file("extdata", "sample-metadata.tsv", package = "mia")
#' treeqza_file <- system.file("extdata", "tree.qza", package = "mia")
#' refseqqza_file <- system.file("extdata", "refseq.qza", package = "mia")
#' tse <- makeTreeSummarizedExperimentFromqiime2(
#'   otu_qza = otuqza_file, taxa_qza = taxaqza_file,
#'   sam_tab = sample_file, refseq_qza = refseqqza_file,
#'   tree_qza = treeqza_file
#' )
#'
#' tse
makeTreeSummarizedExperimentFromqiime2 <- function(otu_qza,
                                                   taxa_qza = NULL,
                                                   sam_tab = NULL,
                                                   refseq_qza = NULL,
                                                   tree_qza = NULL) {
    feature_tab <- read_qza(otu_qza)

    if (!is.null(taxa_qza)) {
        taxa_tab <- read_qza(taxa_qza)
        taxa_tab <- subset_taxa_in_feature(taxa_tab, feature_tab)
        taxa_tab <- parse_q2taxonomy(taxa_tab)
    } else {
        taxa_tab <- NULL
    }

    if (!is.null(sam_tab)) {
        sam_tab <- read_q2sample_meta(sam_tab)
    } else {
        sam_tab <- NULL
    }

    if (!is.null(tree_qza)) {
        tree <- read_qza(tree_qza)
    } else {
        tree <- NULL
    }

    # check the row.names of feature table: DNA sequence or not
    # if is DNA, row.names(feature_tab) is set as refseq
    if  (is_dna_seq(row.names(feature_tab))) {
        refseq <- row.names(feature_tab)
        refseq_nm <- paste0("OTU", seq_along(refseq))
        names(refseq) <- refseq_nm

        # set the rownames of feature and taxa tab as OTU1, OTU2,...
        if (!is.null(taxa_tab)) {
            rownames(taxa_tab) <- refseq_nm
        }
        rownames(feature_tab) <- refseq_nm
    } else {
        refseq <- NULL
    }

    if (!is.null(refseq_qza)){
        refseq <- read_qza(refseq_qza)
    }


    if(!is.null(refseq)) {
        refseq <- Biostrings::DNAStringSet(refseq)
    } else {
        refseq <- NULL
    }

    TreeSummarizedExperiment(
        assays = S4Vectors::SimpleList(counts = feature_tab),
        rowData = taxa_tab,
        colData = sam_tab,
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
#' @keywords internal
read_qza <- function(file, temp = tempdir()) {
    unzipped_file <- utils::unzip(file, exdir = temp)
    meta_file <- grep("metadata.yaml", unzipped_file, value = TRUE)
    metadata <- yaml::read_yaml(meta_file[1])
    format <- metadata$format
    uuid <- metadata$uuid

    if (grepl("BIOMV", metadata$format)) {
        biom_file <- file.path(temp, uuid,"data/feature-table.biom")
        res <- read_q2biom(biom_file)
    } else if (format == "TSVTaxonomyDirectoryFormat") {
        taxa_file <- file.path(temp, uuid, "data/taxonomy.tsv")
        res <- read_q2taxa(taxa_file)
    } else if (format == "NewickDirectoryFormat") {
        tree_file <- file.path(temp, uuid, "data/tree.nwk")
        res <- ape::read.tree(tree_file)
    } else if (format == "DNASequencesDirectoryFormat") {
        refseq_file <- file.path(temp, uuid, "data/dna-sequences.fasta")
        res <- Biostrings::readDNAStringSet(refseq_file)
    } else {
        stop(
            "Only files in format of 'BIOMV210DirFmt', ",
            "'TSVTaxonomyDirectoryFormat', NewickDirectoryFormat' and ",
            "'DNASequencesDirectoryFormat' are supported."
        )
    }

    res
}


#' Read qiime2 feature table
#'
#' Import qiime2 feature table into otu_table object
#'
#' @param file character, file name of the biom file.
#' @return [`phyloseq::otu_table-class`] object.
#' @noRd
#' @keywords internal
read_q2biom <- function(file) {
    biomobj <- suppressWarnings(biomformat::read_biom(file))
    feature_tab <- as(biomformat::biom_data(biomobj),"matrix")

    feature_tab
}

#' Read qiime2 taxa file
#' @keywords internal
#' @noRd
read_q2taxa <- function(file) {
    taxa <- utils::read.table(file, sep = '\t', header = TRUE)
    if ("Confidence" %in% names(taxa)) {
        taxa$Confidence <- NULL
    }

    as.matrix(taxa)
}

#' Read qiime2 sample meta data file
#' @keywords internal
#' @noRd
read_q2sample_meta <- function(file) {
    sam <- utils::read.table(file = file, header = TRUE, sep = "\t", comment.char = "")
    rownames(sam) <- as.character(sam[, 1])

    # remove row: #q2:types
    idx <- which(sam == "#q2:types", arr.ind = TRUE)
    sam <- sam[-idx[, "row"], ]

    S4Vectors::DataFrame(sam)
}


#' Subset taxa according to the taxa in feature table
#' @keywords internal
#' @noRd
subset_taxa_in_feature <- function(taxa_tab, feature_tab) {
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
#'   taxonomic levels, e.g. k__ or D_0__, default `TRUE`.
#' @return [`phyloseq::taxonomyTable-class`] object.
#' @keywords internal
#' @importFrom rlang .data
#' @noRd
parse_q2taxonomy <- function(taxa, sep = "; |;", trim_rank_prefix = TRUE) {
    taxa <- data.frame(taxa)
    if(trim_rank_prefix){
        # remove leading characters from GG
        taxa$Taxon <- gsub("[kpcofgs]__", "", taxa$Taxon)
        #remove leading characters from SILVA
        taxa$Taxon <- gsub("D_\\d__", "", taxa$Taxon)
    }

    taxa <- suppressWarnings(
        tidyr::separate(
            taxa, .data$Taxon,
            c("Kingdom","Phylum","Class","Order","Family","Genus","Species"),
            sep = sep
        )
    )
    taxa <- sapply(taxa, function(x) ifelse(x == "", NA_character_, x))
    rownames(taxa) <- taxa[, "Feature.ID"]
    taxa <- taxa[, setdiff(colnames(taxa), "Feature.ID")]

    taxa
}

#' check the row.names of feature table is DNA sequence or not
#' This function is from https://github.com/YuLab-SMU/MicrobiotaProcess/blob/master/R/import_qiime2.R#L169-L177
#' @keywords internal
#' @noRd
is_dna_seq <- function(names){
    x <- unlist(strsplit(toupper(names[1]), split = ""))
    freq <- mean(x %in% c("A", "G","T","C", "N", "X", "-"))

    if (length(x) > 30 & freq > 0.9) {
        return(TRUE)
    } else {
        return(FALSE)
    }
}
