# Build compact demo objects (.rda) for IBDMDB examples/vignettes.
#
# Produces:
#   data/ibdmdb_2omic_demo.rda  (se_mgx, se_mtx, mae2)
#
# Source raw inputs from inst/extdata and pre-process for speed/size.

message("== IBDMDB demo data preparation ==")

# ------------------------------------------------------------------------------
# Config
# ------------------------------------------------------------------------------

# This script is for DEVELOPERS only.
# It downloads large raw files and prepares a compact demo dataset (.rda).
# It is NOT run during R CMD check.

cache_dir <- file.path("tools", "cache", "mia_ibdmdb_cache")
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

# Stable download URLs (Zenodo)
url_meta <- "https://zenodo.org/records/18280405/files/hmp2_metadata_2018-08-20.csv"
url_mtx  <- "https://zenodo.org/records/18280535/files/ecs_relab.tsv"
url_mgx  <- "https://zenodo.org/records/18280521/files/taxonomic_profiles_mgx.tsv"

f_mgx  <- file.path(cache_dir, "taxonomic_profiles_mgx.tsv")
f_mtx  <- file.path(cache_dir, "ecs_relab.tsv")
f_meta <- file.path(cache_dir, "hmp2_metadata_2018-08-20.csv")

# ------------------------------------------------------------------------------
# Dependencies
# ------------------------------------------------------------------------------

need <- c("data.table", "matrixStats", "SummarizedExperiment", "MultiAssayExperiment")
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
        quote  = ""
    )
    if (ncol(dt) != length(header_vec)) {
        stop(sprintf(
            "Header columns (%d) != data columns (%d) in %s",
            length(header_vec), ncol(dt), path
        ))
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

# Attach metadata rows to a SummarizedExperiment, if possible
make_SE <- function(mat, meta_df = NULL, assay_name = "counts") {
    if (is.null(meta_df)) {
        return(SummarizedExperiment::SummarizedExperiment(
            assays  = setNames(list(mat), assay_name),
            colData = S4Vectors::DataFrame(row.names = colnames(mat))
        ))
    }
    
    md <- as.data.frame(meta_df, stringsAsFactors = FALSE, check.names = FALSE)
    
    # Find best matching ID column between metadata and matrix samples
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
        colData = S4Vectors::DataFrame(md_sub)
    )
}

# ------------------------------------------------------------------------------
# Download raw inputs (if missing)
# ------------------------------------------------------------------------------

download_if_missing(url_mgx,  f_mgx)
download_if_missing(url_mtx,  f_mtx)
download_if_missing(url_meta, f_meta)

# ------------------------------------------------------------------------------
# Metadata (shared with 2-omic demo)
# ------------------------------------------------------------------------------

meta_full <- read_metadata(f_meta)

demo_samples <- character(0)

# ------------------------------------------------------------------------------
# Prepare 2-omic (MGX + MTX)
# ------------------------------------------------------------------------------

message("Preparing 2-omic MGX + MTX demo ...")

message("== Done. Re-run devtools::document(); devtools::check(); BiocCheck::BiocCheck(). ==")