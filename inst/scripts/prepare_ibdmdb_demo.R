#build compact demo objects (.rda) for IBDMDB examples/vignettes.
#
#produces:
#  data/ibdmdb_2omic_demo.rda  (se_mgx, se_mtx, mae2)
#  data/ibdmdb_meta_demo.rda   (ibdmdb_meta_demo: sample-level metadata subset)
#
#source raw inputs from inst/extdata and pre-process for speed/size.

message("== IBDMDB demo data preparation ==")

# ------------------------------------------------------------------------------
#config
# ------------------------------------------------------------------------------

raw_dir <- file.path("inst", "extdata")

#2-omic (used in ibdmdb_benchmarking.qmd & ibdmdb_2omic_jointrpca.qmd)
f_mgx  <- file.path(raw_dir, "taxonomic_profiles_mgx.tsv")
f_mtx  <- file.path(raw_dir, "ecs_relab.tsv")
f_meta <- file.path(raw_dir, "hmp2_metadata_2018-08-20.csv")

#prevalence thresholds (fraction of samples)
prev_mgx_frac <- 0.05
prev_mtx_frac <- 0.02

#cap feature counts for speed/size
cap_mgx <- 800L
cap_mtx <- 800L
max_samples <- 60L

# ------------------------------------------------------------------------------
#dependencies
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
#helpers
# ------------------------------------------------------------------------------

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

#attach metadata rows to a SummarizedExperiment, if possible
make_SE <- function(mat, meta_df = NULL, assay_name = "counts") {
    if (is.null(meta_df)) {
        return(SummarizedExperiment::SummarizedExperiment(
            assays  = setNames(list(mat), assay_name),
            colData = S4Vectors::DataFrame(row.names = colnames(mat))
        ))
    }
    
    md <- as.data.frame(meta_df, stringsAsFactors = FALSE, check.names = FALSE)
    
    #find best matching ID column between metadata and matrix samples
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
#I/O guards
# ------------------------------------------------------------------------------

if (!dir.exists(raw_dir)) stop("Raw input dir not found: ", raw_dir)
if (!dir.exists("data")) dir.create("data", recursive = TRUE)

# ------------------------------------------------------------------------------
#metadata (shared with 2-omic demo)
# ------------------------------------------------------------------------------

has_meta <- file.exists(f_meta)
meta_full <- NULL
if (has_meta) {
    meta_full <- read_metadata(f_meta)
}

demo_samples <- character(0)

# ------------------------------------------------------------------------------
#prepare 2-omic (MGX + MTX)
# ------------------------------------------------------------------------------

has_mgx <- file.exists(f_mgx)
has_mtx <- file.exists(f_mtx)

if (has_mgx && has_mtx) {
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
    
    #per-view prevalence
    n_samp <- ncol(M_mgx)
    keep_mgx <- rowSums(M_mgx > 0) >= ceiling(prev_mgx_frac * n_samp)
    keep_mtx <- rowSums(M_mtx > 0) >= ceiling(prev_mtx_frac * n_samp)
    M_mgx <- M_mgx[keep_mgx, , drop = FALSE]
    M_mtx <- M_mtx[keep_mtx, , drop = FALSE]
    
    #drop all-zero samples per view
    M_mgx <- M_mgx[, colSums(M_mgx) > 0, drop = FALSE]
    M_mtx <- M_mtx[, colSums(M_mtx) > 0, drop = FALSE]
    
    #recompute strict shared
    shared2 <- intersect(colnames(M_mgx), colnames(M_mtx))
    shared2 <- sort(unique(shared2))
    set.seed(1)
    if (length(shared2) > max_samples) {
        shared2 <- sort(sample(shared2, max_samples))
    }
    M_mgx <- M_mgx[, shared2, drop = FALSE]
    M_mtx <- M_mtx[, shared2, drop = FALSE]
    
    #cap by variance
    M_mgx <- cap_by_var(M_mgx, cap_mgx)
    M_mtx <- cap_by_var(M_mtx, cap_mtx)
    
    #attach metadata if available
    meta_df <- if (has_meta) meta_full else NULL
    se_mgx  <- make_SE(M_mgx, meta_df, assay_name = "mgx")
    se_mtx  <- make_SE(M_mtx, meta_df, assay_name = "mtx")
    
    mae2 <- MultiAssayExperiment::MultiAssayExperiment(
        experiments = list(MGX = se_mgx, MTX = se_mtx)
    )
    
    #attach sample metadata at the MAE level (if consistent across experiments)
    if (!is.null(colData(se_mgx))) {
        MultiAssayExperiment::colData(mae2) <- colData(se_mgx)
    }
    
    #track demo sample IDs
    demo_samples <- union(demo_samples, colnames(M_mgx))
    
    #save demo objects
    save(
        se_mgx, se_mtx, mae2,
        file     = file.path("data", "ibdmdb_2omic_demo.rda"),
        compress = "xz"
    )
    
    message("Saved: data/ibdmdb_2omic_demo.rda")
} else {
    message(
        "Skipping 2-omic demo (missing files):",
        "\n  MGX: ", f_mgx, " (", has_mgx, ")",
        "\n  MTX: ", f_mtx, " (", has_mtx, ")"
    )
}

# ------------------------------------------------------------------------------
# save ibdmdb_meta_demo (metadata subset for demo samples)
# ------------------------------------------------------------------------------

if (has_meta) {
    if (exists("mae2")) {
        ibdmdb_meta_demo <- as.data.frame(MultiAssayExperiment::colData(mae2))
        message("Saving ibdmdb_meta_demo derived from colData(mae2).")
    } else {
        ibdmdb_meta_demo <- meta_full
        message("mae2 not available; saving full metadata as ibdmdb_meta_demo.")
    }
    
    save(
        ibdmdb_meta_demo,
        file     = file.path("data", "ibdmdb_meta_demo.rda"),
        compress = "xz"
    )
    message("Saved: data/ibdmdb_meta_demo.rda")
}

message("== Done. Re-run devtools::document(); devtools::check(); BiocCheck::BiocCheck(). ==")