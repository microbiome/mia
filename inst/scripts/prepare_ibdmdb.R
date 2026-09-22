# Build compact demo objects (.rda) for IBDMDB examples/vignettes.
#
# Produces:
#   data/ibdmdb.rda  (ibdmdb)
#
# Source raw inputs from inst/extdata and pre-process for speed/size.
#
# This script is for DEVELOPERS only. It is NOT run during R CMD check.

message("== IBDMDB demo data preparation ==")

# ------------------------------------------------------------------------------
# Config
# ------------------------------------------------------------------------------

# This script is for DEVELOPERS only.
# It downloads large raw files into a local cache (NOT committed to git),
# then creates a compact demo dataset in data/ibdmdb.rda.
# It is NOT run during R CMD check.

cache_dir <- file.path("tools", "cache", "mia_ibdmdb_cache")
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

# Stable download URLs (Zenodo)
url_meta <- "https://zenodo.org/records/18280405/files/hmp2_metadata_2018-08-20.csv"
url_mtx  <- "https://zenodo.org/records/18280535/files/ecs_relab.tsv"
url_mgx  <- "https://zenodo.org/records/18280521/files/taxonomic_profiles_mgx.tsv"

# Local cached file paths (downloaded if missing)
f_mgx  <- file.path(cache_dir, "taxonomic_profiles_mgx.tsv")
f_mtx  <- file.path(cache_dir, "ecs_relab.tsv")
f_meta <- file.path(cache_dir, "hmp2_metadata_2018-08-20.csv")

# Prevalence thresholds (fraction of samples)
prev_mgx_frac <- 0.05
prev_mtx_frac <- 0.02

# Cap feature counts for speed/size
cap_mgx <- 800L
cap_mtx <- 800L
max_samples <- 60L

# ------------------------------------------------------------------------------
# Dependencies
# ------------------------------------------------------------------------------

need <- c("data.table", "matrixStats", "SummarizedExperiment", "MultiAssayExperiment",
          "S4Vectors")
miss <- setdiff(need, rownames(installed.packages()))
if (length(miss)) {
    stop(
        "Missing packages: ",
        paste(miss, collapse = ", "),
        ". Install them, then re-run this script."
    )
}

library(data.table)
library(matrixStats)
library(SummarizedExperiment)
library(MultiAssayExperiment)
library(S4Vectors)

# ------------------------------------------------------------------------------
# Helpers
# ------------------------------------------------------------------------------

download_if_missing <- function(url, dest) {
    if (file.exists(dest)) return(invisible(dest))
    message("Downloading: ", basename(dest))
    utils::download.file(url, dest, mode = "wb", quiet = FALSE)
    invisible(dest)
}

read_ibdmdb_tsv <- function(path) {
    stopifnot(file.exists(path))
    first <- readLines(path, n = 200L, warn = FALSE)
    comment_idx <- which(grepl("^#", first))
    if (length(comment_idx) == 0L) {
        stop("No commented header line found in: ", path)
    }

    header_line <- first[max(comment_idx)]
    header_line <- sub("^#\\s*", "", header_line)
    header_line <- sub("^\ufeff", "", header_line)
    header_vec  <- strsplit(header_line, "\t", fixed = TRUE)[[1]]
    header_vec  <- gsub('^"|"$', "", header_vec)

    dt <- data.table::fread(
        path,
        skip   = length(comment_idx),
        header = FALSE,
        sep    = "\t",
        quote  = "",
        fill   = TRUE
    )
    if (ncol(dt) < length(header_vec)) {
        stop(sprintf(
            "Header columns (%d) > data columns (%d) in %s",
            length(header_vec), ncol(dt), path
        ))
    }
    if (ncol(dt) > length(header_vec)) {
        dt <- dt[, seq_len(length(header_vec)), with = FALSE]
    }
    data.table::setnames(dt, header_vec)
    dt
}

to_matrix <- function(dt) {
    rn  <- dt[[1]]
    mat <- as.matrix(dt[, -1, with = FALSE])
    rownames(mat) <- rn
    storage.mode(mat) <- "numeric"
    mat[is.na(mat)] <- 0
    mat
}

dedup_rownames <- function(mat) {
    stopifnot(!is.null(rownames(mat)))
    rn <- rownames(mat)
    empty_or_na <- rn == "" | is.na(rn)
    if (any(empty_or_na)) {
        rn[empty_or_na] <- paste0("feat_", seq_len(sum(empty_or_na)))
    }
    rownames(mat) <- make.unique(as.character(rn), sep = "_")
    mat
}

sanitize_matrix <- function(mat) {
    mat <- as.matrix(mat)
    storage.mode(mat) <- "numeric"
    mat[!is.finite(mat)] <- 0
    mat
}

cap_by_var <- function(M, cap) {
    if (nrow(M) <= cap) return(M)
    ord <- order(matrixStats::rowVars(M), decreasing = TRUE)
    M[ord[seq_len(cap)], , drop = FALSE]
}

read_metadata <- function(path) {
    stopifnot(file.exists(path))
    read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

make_SE <- function(mat, meta_df = NULL, assay_name = "counts") {
    if (is.null(meta_df)) {
        return(SummarizedExperiment::SummarizedExperiment(
            assays  = setNames(list(mat), assay_name),
            colData = S4Vectors::DataFrame(row.names = colnames(mat))
        ))
    }

    md <- meta_df

    overlaps <- vapply(
        md,
        function(col) sum(as.character(col) %in% colnames(mat)),
        numeric(1)
    )
    best <- names(overlaps)[which.max(overlaps)]

    if (length(best) == 0 || overlaps[[best]] == 0) {
        return(SummarizedExperiment::SummarizedExperiment(
            assays  = setNames(list(mat), assay_name),
            colData = S4Vectors::DataFrame(row.names = colnames(mat))
        ))
    }

    md_sub <- md[md[[best]] %in% colnames(mat), , drop = FALSE]
    md_sub <- md_sub[!duplicated(md_sub[[best]]), , drop = FALSE]

    rownames(md_sub) <- as.character(md_sub[[best]])
    md_sub <- md_sub[colnames(mat), , drop = FALSE]

    if (anyDuplicated(rownames(md_sub))) {
        rownames(md_sub) <- make.unique(rownames(md_sub), sep = "_dup")
    }
    stopifnot(identical(rownames(md_sub), colnames(mat)))

    SummarizedExperiment::SummarizedExperiment(
        assays  = setNames(list(mat), assay_name),
        colData = md_sub
    )
}

# ------------------------------------------------------------------------------
# I/O guards + download raw inputs
# ------------------------------------------------------------------------------

if (!dir.exists("data")) dir.create("data", recursive = TRUE)

download_if_missing(url_mgx,  f_mgx)
download_if_missing(url_mtx,  f_mtx)
download_if_missing(url_meta, f_meta)

missing_files <- c(f_mgx, f_mtx, f_meta)[!file.exists(c(f_mgx, f_mtx, f_meta))]
if (length(missing_files)) {
    stop(
        "Missing cached file(s):\n- ",
        paste(missing_files, collapse = "\n- "),
        "\n\nDownload failed or paths are wrong."
    )
}

# ------------------------------------------------------------------------------
# Metadata
# ------------------------------------------------------------------------------

meta_full <- read_metadata(f_meta)
meta_full <- DataFrame(meta_full, check.names = FALSE)

# ------------------------------------------------------------------------------
# Prepare 2-omic (MGX + MTX)
# ------------------------------------------------------------------------------

message("Preparing 2-omic MGX + MTX demo ...")

dt_mgx <- read_ibdmdb_tsv(f_mgx)
dt_mtx <- read_ibdmdb_tsv(f_mtx)

M_mgx <- sanitize_matrix(dedup_rownames(to_matrix(dt_mgx)))
M_mtx <- sanitize_matrix(dedup_rownames(to_matrix(dt_mtx)))

shared <- intersect(colnames(M_mgx), colnames(M_mtx))
shared <- sort(unique(shared[nchar(shared) > 0]))
if (length(shared) < 20) {
    warning("Few shared samples for 2-omic: ", length(shared), " (keeping anyway).")
}
M_mgx <- M_mgx[, shared, drop = FALSE]
M_mtx <- M_mtx[, shared, drop = FALSE]

n_samp <- ncol(M_mgx)
keep_mgx <- rowSums(M_mgx > 0) >= ceiling(prev_mgx_frac * n_samp)
keep_mtx <- rowSums(M_mtx > 0) >= ceiling(prev_mtx_frac * n_samp)
M_mgx <- M_mgx[keep_mgx, , drop = FALSE]
M_mtx <- M_mtx[keep_mtx, , drop = FALSE]

M_mgx <- M_mgx[, colSums(M_mgx) > 0, drop = FALSE]
M_mtx <- M_mtx[, colSums(M_mtx) > 0, drop = FALSE]

shared2 <- intersect(colnames(M_mgx), colnames(M_mtx))
shared2 <- sort(unique(shared2))
set.seed(1)
if (length(shared2) > max_samples) {
    shared2 <- sort(sample(shared2, max_samples))
}
M_mgx <- M_mgx[, shared2, drop = FALSE]
M_mtx <- M_mtx[, shared2, drop = FALSE]

M_mgx <- cap_by_var(M_mgx, cap_mgx)
M_mtx <- cap_by_var(M_mtx, cap_mtx)

meta_df <- meta_full[meta_full$data_type %in% "metagenomics", ]
se_mgx  <- make_SE(M_mgx, meta_df, assay_name = "mgx")
# Split taxa into taxonomy
rd <- mia:::.parse_taxonomy(data.frame(Taxon = rownames(se_mgx)), sep = "\\|", remove.prefix = TRUE)
rownames(rd) <- rownames(se_mgx)
rowData(se_mgx) <- rd
se_mgx <- agglomerateByRanks(se_mgx)
se_mgx <- swapAltExp(se_mgx, "species", "original")

meta_df <- meta_full[meta_full$data_type %in% "metatranscriptomics", ]
se_mtx  <- make_SE(M_mtx, meta_df, assay_name = "mtx")
# Create a feature metadata table
df <- data.frame(
    feature = rownames(se_mtx),
    stringsAsFactors = FALSE
)
df$gene_function <- sub("\\|.*$", "", df$feature)
df$taxon <- ifelse(
    grepl("\\|", df$feature),
    sub("^.*\\|", "", df$feature),
    NA_character_
)
# Create function + bacteria identifier
df$gene_taxon <- ifelse(
    !is.na(df$taxon),
    paste(df$gene_function, df$taxon, sep = "|"),
    NA_character_
)
df <- DataFrame(df)
rownames(df) <- rownames(se_mtx)
rowData(se_mtx) <- df
se_mtx <- as(se_mtx, "TreeSummarizedExperiment")

# Function level
altExp(se_mtx, "gene_function") <- agglomerateByVariable(
    se_mtx,
    by = 1,
    group = "gene_function"
)
# Function + bacteria level
altExp(se_mtx, "gene_taxon") <- agglomerateByVariable(
    se_mtx,
    by = 1,
    group = "gene_taxon"
)

mae <- MultiAssayExperiment::MultiAssayExperiment(
    experiments = list(MGX = se_mgx, MTX = se_mtx)
)

if (!is.null(SummarizedExperiment::colData(se_mgx))) {
    MultiAssayExperiment::colData(mae) <- SummarizedExperiment::colData(se_mgx)
}

# Save only one dataset object
ibdmdb <- mae
save(
    ibdmdb,
    file     = file.path("data", "ibdmdb.rda"),
    compress = "xz"
)

message("Saved: data/ibdmdb.rda")
message("== Done. Re-run devtools::document(); devtools::check(); BiocCheck::BiocCheck(). ==")
