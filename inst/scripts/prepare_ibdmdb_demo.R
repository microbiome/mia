#build compact demo objects (.rda) for IBDMDB examples/vignettes.

message("== IBDMDB demo data preparation ==")

#config
raw_dir <- file.path("inst", "extdata")

#2-omic (used in ibdmdb_benchmarking.qmd & ibdmdb_2omic_jointrpca.qmd)
f_mgx   <- file.path(raw_dir, "taxonomic_profiles_mgx.tsv")     
f_mtx   <- file.path(raw_dir, "ecs_relab.tsv")                   
f_meta  <- file.path(raw_dir, "hmp2_metadata_2018-08-20.csv")

#3-omic (used in ibdmdb_3omic_jointrpca.qmd)
f_16s_n <- file.path(raw_dir, "taxonomic_profiles_16s_new.tsv")
f_mgx_n <- file.path(raw_dir, "taxonomic_profiles_mgx_new.tsv")
f_mtx_n <- file.path(raw_dir, "taxonomic_profiles_mtx_new.tsv")

#prevalence thresholds (in samples)
prev_16s_frac <- 0.05
prev_mgx_frac <- 0.05
prev_mtx_frac <- 0.02

#cap feature counts for speed/size 
cap_16s <- 10000L
cap_mgx <- 10000L
cap_mtx <- 10000L

#dependencies
need <- c("data.table", "matrixStats", "SummarizedExperiment", "MultiAssayExperiment")
miss <- setdiff(need, rownames(installed.packages()))
if (length(miss)) {
  stop("Missing packages: ", paste(miss, collapse = ", "),
       ". Install them, then re-run this script.")
}
library(data.table)
library(matrixStats)
library(SummarizedExperiment)
library(MultiAssayExperiment)

#helpers
read_ibdmdb_tsv <- function(path) {
  stopifnot(file.exists(path))
  first <- readLines(path, n = 200L, warn = FALSE)
  comment_idx <- which(grepl("^#", first))
  if (length(comment_idx) == 0L) stop("No commented header line found in: ", path)
  
  header_line <- first[max(comment_idx)]
  header_line <- sub("^#\\s*", "", header_line)
  header_line <- sub("^\ufeff", "", header_line)
  header_vec  <- strsplit(header_line, "\t", fixed = TRUE)[[1]]
  header_vec  <- gsub('^"|"$', "", header_vec)
  
  dt <- data.table::fread(path, skip = length(comment_idx), header = FALSE, sep = "\t", quote = "")
  if (ncol(dt) != length(header_vec)) {
    stop(sprintf("Header columns (%d) != data columns (%d) in %s",
                 length(header_vec), ncol(dt), path))
  }
  data.table::setnames(dt, header_vec)
  dt
}

to_matrix <- function(dt) {
  rn <- dt[[1]]
  mat <- as.matrix(dt[, -1, with = FALSE])
  rownames(mat) <- rn
  storage.mode(mat) <- "numeric"
  mat[is.na(mat)] <- 0
  mat
}

dedup_rownames <- function(mat) {
  stopifnot(!is.null(rownames(mat)))
  rn <- rownames(mat)
  rn[rn == "" | is.na(rn)] <- paste0("feat_", seq_len(sum(rn == "" | is.na(rn))))
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
  meta <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  meta
}

make_SE <- function(mat, meta_df = NULL, assay_name = "counts") {
  if (is.null(meta_df)) {
    return(SummarizedExperiment::SummarizedExperiment(
      assays  = setNames(list(mat), assay_name),
      colData = S4Vectors::DataFrame(row.names = colnames(mat))
    ))
  }
  
  md <- as.data.frame(meta_df, stringsAsFactors = FALSE, check.names = FALSE)
  
  overlaps <- vapply(md, function(col) sum(as.character(col) %in% colnames(mat)), numeric(1))
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

#I/O guards 
if (!dir.exists(raw_dir)) stop("Raw input dir not found: ", raw_dir)
if (!dir.exists("data")) dir.create("data", recursive = TRUE)

#prepare 2-omic (MGX + MTX) 
has_mgx <- file.exists(f_mgx)
has_mtx <- file.exists(f_mtx)
has_meta <- file.exists(f_meta)

if (has_mgx && has_mtx) {
  message("Preparing 2-omic MGX+MTX demo ...")
  
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
  
  #per-view prevalence
  n_samp <- ncol(M_mgx)
  keep_mgx <- rowSums(M_mgx > 0) >= ceiling(prev_mgx_frac * n_samp)
  keep_mtx <- rowSums(M_mtx > 0) >= ceiling(prev_mtx_frac * n_samp)
  M_mgx <- M_mgx[keep_mgx, , drop = FALSE]
  M_mtx <- M_mtx[keep_mtx, , drop = FALSE]
  
  #drop all-zero samples per-view
  M_mgx <- M_mgx[, colSums(M_mgx) > 0, drop = FALSE]
  M_mtx <- M_mtx[, colSums(M_mtx) > 0, drop = FALSE]
  
  #recompute strict shared
  shared2 <- intersect(colnames(M_mgx), colnames(M_mtx))
  shared2 <- sort(unique(shared2))
  M_mgx <- M_mgx[, shared2, drop = FALSE]
  M_mtx <- M_mtx[, shared2, drop = FALSE]
  
  #cap by variance
  M_mgx <- cap_by_var(M_mgx, cap_mgx)
  M_mtx <- cap_by_var(M_mtx, cap_mtx)
  
  meta_df <- if (has_meta) read_metadata(f_meta) else NULL
  se_mgx  <- make_SE(M_mgx, meta_df, assay_name = "mgx")
  se_mtx  <- make_SE(M_mtx, meta_df, assay_name = "mtx")
  
  mae2 <- MultiAssayExperiment::MultiAssayExperiment(
    experiments = list(MGX = se_mgx, MTX = se_mtx)
  )
  
  #save objects
  save(se_mgx, se_mtx, mae2, file = file.path("data", "ibdmdb_2omic_demo.rda"), compress = "xz")
  if (has_meta) {
    ibdmdb_meta_demo <- meta_df
    save(ibdmdb_meta_demo, file = file.path("data", "ibdmdb_meta_demo.rda"), compress = "xz")
  }
  message("Saved: data/ibdmdb_2omic_demo.rda", if (has_meta) " and data/ibdmdb_meta_demo.rda" else "")
} else {
  message("Skipping 2-omic demo (missing files):",
          "\n  MGX: ", f_mgx, " (", has_mgx, ")",
          "\n  MTX: ", f_mtx, " (", has_mtx, ")")
}

#prepare 3-omic (16S + MGX_new + MTX_new) 
has_16s_n <- file.exists(f_16s_n)
has_mgx_n <- file.exists(f_mgx_n)
has_mtx_n <- file.exists(f_mtx_n)

if (has_16s_n && has_mgx_n && has_mtx_n) {
  message("Preparing 3-omic 16S+MGX_new+MTX_new demo ...")
  
  dt_16s <- read_ibdmdb_tsv(f_16s_n)
  dt_mgx2 <- read_ibdmdb_tsv(f_mgx_n)
  dt_mtx2 <- read_ibdmdb_tsv(f_mtx_n)
  
  X_16s <- sanitize_matrix(dedup_rownames(to_matrix(dt_16s)))
  X_mgx2 <- sanitize_matrix(dedup_rownames(to_matrix(dt_mgx2)))
  X_mtx2 <- sanitize_matrix(dedup_rownames(to_matrix(dt_mtx2)))
  
  shared3 <- Reduce(intersect, list(colnames(X_16s), colnames(X_mgx2), colnames(X_mtx2)))
  shared3 <- sort(unique(shared3[nchar(shared3) > 0]))
  if (length(shared3) < 10) {
    warning("Few shared samples for 3-omic: ", length(shared3), " (keeping anyway).")
  }
  X_16s  <- X_16s[,  shared3, drop = FALSE]
  X_mgx2 <- X_mgx2[, shared3, drop = FALSE]
  X_mtx2 <- X_mtx2[, shared3, drop = FALSE]
  
  #per-view prevalence
  n_samp3 <- length(shared3)
  keep_16s  <- rowSums(X_16s  > 0) >= ceiling(prev_16s_frac * n_samp3)
  keep_mgx2 <- rowSums(X_mgx2 > 0) >= ceiling(prev_mgx_frac * n_samp3)
  keep_mtx2 <- rowSums(X_mtx2 > 0) >= ceiling(prev_mtx_frac * n_samp3)
  
  X_16s  <- X_16s[keep_16s, , drop = FALSE]
  X_mgx2 <- X_mgx2[keep_mgx2, , drop = FALSE]
  X_mtx2 <- X_mtx2[keep_mtx2, , drop = FALSE]
  
  #drop all-zero samples per-view
  X_16s  <- X_16s[,  colSums(X_16s)  > 0, drop = FALSE]
  X_mgx2 <- X_mgx2[, colSums(X_mgx2) > 0, drop = FALSE]
  X_mtx2 <- X_mtx2[, colSums(X_mtx2) > 0, drop = FALSE]
  
  #recompute strict shared
  sharedf <- Reduce(intersect, list(colnames(X_16s), colnames(X_mgx2), colnames(X_mtx2)))
  sharedf <- sort(unique(sharedf))
  X_16s  <- X_16s[,  sharedf, drop = FALSE]
  X_mgx2 <- X_mgx2[, sharedf, drop = FALSE]
  X_mtx2 <- X_mtx2[, sharedf, drop = FALSE]
  
  #cap by variance
  X_16s  <- cap_by_var(X_16s,  cap_16s)
  X_mgx2 <- cap_by_var(X_mgx2, cap_mgx)
  X_mtx2 <- cap_by_var(X_mtx2, cap_mtx)
  
  meta_df <- if (has_meta) read_metadata(f_meta) else NULL
  se_16s  <- make_SE(X_16s,  meta_df, assay_name = "abundance")
  se_mgx2 <- make_SE(X_mgx2, meta_df, assay_name = "abundance")
  se_mtx2 <- make_SE(X_mtx2, meta_df, assay_name = "abundance")
  
  mae3 <- MultiAssayExperiment::MultiAssayExperiment(
    experiments = list(`16S` = se_16s, MGX = se_mgx2, MTX = se_mtx2)
  )
  
  save(se_16s, se_mgx2, se_mtx2, mae3, file = file.path("data", "ibdmdb_3omic_demo.rda"), compress = "xz")
  message("Saved: data/ibdmdb_3omic_demo.rda")
} else {
  message("Skipping 3-omic demo (missing files):",
          "\n  16S_new: ", f_16s_n, " (", has_16s_n, ")",
          "\n  MGX_new: ", f_mgx_n, " (", has_mgx_n, ")",
          "\n  MTX_new: ", f_mtx_n, " (", has_mtx_n, ")")
}

message("== Done. Re-run devtools::document(); devtools::check(); BiocCheck::BiocCheck(). ==")