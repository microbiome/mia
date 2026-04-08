#' Fetch data from the Metalog database as
#' \code{TreeSummarizedExperiment}
#'
#' \code{fetchMetalogTSE} downloads MetaPhlAn4 taxonomic profiles and
#' associated sample metadata from the
#' \href{https://metalog.embl.de/}{Metalog} database and returns them as a
#' \code{TreeSummarizedExperiment} object. An optional sample list can be
#' provided to retrieve only a subset of samples (exported from the Metalog
#' web UI).
#'
#' @param collection \code{Character vector}. One or more Metalog collections
#'   to download. Each value must be one of \code{"human"}, \code{"animal"},
#'   \code{"ocean"}, or \code{"environmental"}. When multiple collections are
#'   given, their assays and metadata are merged before constructing the
#'   final object; a \code{collection} column is added to \code{colData}
#'   identifying the source collection of each sample.
#'
#' @param meta.type \code{Character scalar}. The metadata scope to download.
#'   Must be one of \code{"core"}, \code{"extended"}, \code{"all"}, or
#'   \code{"none"}. When \code{"none"}, no sample metadata is downloaded and
#'   \code{colData} is empty (in the multi-collection case it retains only
#'   the \code{collection} column identifying the source of each sample).
#'   (Default: \code{"core"}).
#'
#' @param samplelist \code{Character scalar} or \code{NULL}. File path to a
#'   sample list exported from the Metalog web UI. Supported formats are
#'   \code{csv}, \code{tsv}, and \code{json}. When provided, the returned
#'   object is filtered to include only the listed samples; taxa with zero
#'   abundance across the retained samples are also removed.
#'   (Default: \code{NULL}).
#'
#' @param use.cache \code{Logical scalar}. Should previously downloaded files
#'   be reused? When \code{TRUE}, cached files in the download directory are
#'   used if available. (Default: \code{TRUE}).
#'
#' @param make.dense \code{Logical scalar}. Should the assay be returned as a
#'   dense base R matrix? Internally the assay is built as a sparse matrix to
#'   conserve memory while loading and merging collections. When \code{TRUE}
#'   (the default), it is converted to a dense matrix as a final step for
#'   broader compatibility with downstream tools. Set to \code{FALSE} to
#'   retain the sparse representation. (Default: \code{TRUE}).
#'
#' @details
#' Data is downloaded from the Metalog database
#' (\url{https://metalog.embl.de/}) and returned as a
#' \code{TreeSummarizedExperiment}. The assay is stored as a sparse matrix
#' named \code{"relabundance"} containing MetaPhlAn4 relative abundances at
#' SGB (species-level genome bin) resolution. Full taxonomic lineages are
#' mapped to standard ranks (Kingdom through SGB) and stored in
#' \code{rowData}. Sample metadata (in long format on the server) is pivoted
#' to wide format and stored in \code{colData}.
#'
#' This function requires an internet connection to download data from the
#' Metalog server.
#'
#' Provenance information is stored in \code{metadata(tse)$metalog} as a
#' list containing the source URL, license, collection and metadata type
#' used, the date stamps parsed from the downloaded file names
#' (\code{date_profile}, \code{date_metadata}), and the date the data was
#' fetched (\code{date_fetched}).
#'
#' @return A
#' \code{\link[TreeSummarizedExperiment:TreeSummarizedExperiment-class]{TreeSummarizedExperiment}}
#' object
#'
#' @name fetchMetalogTSE
#' @seealso
#' \code{\link[=importMetaPhlAn]{importMetaPhlAn}}
#'
#' @export
#'
#' @author Rasmus Hindström
#'
#' @references
#' Metalog database: \url{https://metalog.embl.de/}
#'
#' The data is made available under the Open Database License (ODbL) v1.0.
#'
#' @examples
#' \dontrun{
#' # Fetch the human collection with core metadata
#' tse <- fetchMetalogTSE("human")
#'
#' # Fetch with extended metadata
#' tse <- fetchMetalogTSE("human", meta.type = "extended")
#'
#' # Fetch a subset of samples using a sample list
#' tse <- fetchMetalogTSE("human", samplelist = "my_samples.csv")
#'
#' # Fetch and merge multiple collections
#' tse <- fetchMetalogTSE(c("human", "animal"))
#' }
#'
NULL

#' @rdname fetchMetalogTSE
#' @importFrom TreeSummarizedExperiment TreeSummarizedExperiment
#' @importFrom S4Vectors SimpleList DataFrame metadata metadata<-
#' @importFrom data.table fread setnames setkey dcast tstrsplit
#' @importFrom ape read.tree keep.tip
#' @export
fetchMetalogTSE <- function(
        collection,
        meta.type = "core",
        samplelist = NULL,
        use.cache = TRUE,
        make.dense = TRUE) {
    ################################ Input check ###############################
    allowed_collections <- c("human", "animal", "ocean", "environmental")
    if (!is.character(collection) || length(collection) < 1 ||
            anyNA(collection) || !all(nzchar(collection)) ||
            !all(collection %in% allowed_collections)) {
        stop(
            "'collection' must be a character vector with values from: ",
            paste(dQuote(allowed_collections), collapse = ", "),
            call. = FALSE
        )
    }
    collection <- unique(collection)
    allowed_meta_types <- c("core", "extended", "all", "none")
    if (!.is_non_empty_string(meta.type) ||
            !meta.type %in% allowed_meta_types) {
        stop(
            "'meta.type' must be one of: ",
            paste(dQuote(allowed_meta_types), collapse = ", "),
            call. = FALSE
        )
    }
    if (!is.null(samplelist) && !.is_non_empty_string(samplelist)) {
        stop("'samplelist' must be a single character value or NULL.",
            call. = FALSE)
    }
    if (!is.null(samplelist)) {
        if (!file.exists(samplelist)) {
            stop("'samplelist' file does not exist: ", samplelist,
                call. = FALSE)
        }
        ext <- tolower(tools::file_ext(samplelist))
        allowed_exts <- c("csv", "tsv", "json")
        if (!ext %in% allowed_exts) {
            stop(
                "'samplelist' file type must be one of: ",
                paste(allowed_exts, collapse = ", "),
                call. = FALSE
            )
        }
    }
    if (!.is_a_bool(use.cache)) {
        stop("'use.cache' must be TRUE or FALSE.", call. = FALSE)
    }
    if (!.is_a_bool(make.dense)) {
        stop("'make.dense' must be TRUE or FALSE.", call. = FALSE)
    }
    ############################## Input check end #############################
    # Latest database file for taxonomy mapping (shared across collections)
    mapping_db <- .download_metalog_file(
        "https://metalog.embl.de/static/download/profiles/metaphlan4_clades.tsv.gz",
        use.cache = use.cache
    )
    # Per-collection: download, load assay, optional sample filter, load md
    data_files_list <- lapply(collection, .resolve_metalog_url,
        meta.type = meta.type, use.cache = use.cache)
    names(data_files_list) <- collection
    per_coll <- lapply(collection, function(co) {
        df <- data_files_list[[co]]
        al <- .load_metalog_assay(df[["assay"]])
        if (!is.null(samplelist)) {
            message("Filtering to requested samples in '", co, "'...")
            al <- .filter_metalog_samples(al, samplelist)
        }
        md <- if (identical(meta.type, "none")) {
            data.frame(row.names = al[["samples"]])
        } else {
            .load_metalog_metadata(df[["md"]], al[["samples"]])
        }
        list(assay_list = al, md = md)
    })
    names(per_coll) <- collection
    # Merge assays and metadata across collections
    merged <- .merge_metalog_assays(lapply(per_coll, `[[`, "assay_list"))
    if (make.dense) {
        merged[["assay"]] <- as.matrix(merged[["assay"]])
    }
    md_df <- .merge_metalog_metadata(
        lapply(per_coll, `[[`, "md"), collection, merged[["samples"]],
        add.collection = !(identical(meta.type, "none") &&
            length(collection) == 1))
    # Map SGBs to full taxonomic lineage
    tax <- .construct_metalog_taxmap(mapping_db, merged[["taxa"]])
    # Download and prune the MetaPhlAn4 SGB phylogeny to taxa in the data
    tree_info <- .construct_metalog_tree(merged[["taxa"]], use.cache)

    tse <- TreeSummarizedExperiment(
        assays = SimpleList("relabundance" = merged[["assay"]]),
        colData = DataFrame(md_df),
        rowData = DataFrame(tax),
        rowTree = tree_info[["tree"]],
        rowNodeLab = tree_info[["node_lab"]]
    )

    # Store provenance information
    date_profile <- vapply(data_files_list,
        function(df) .parse_metalog_date(df[["assay"]]), character(1))
    date_metadata <- vapply(data_files_list, function(df) {
        if (is.na(df[["md"]])) NA_character_
        else .parse_metalog_date(df[["md"]])
    }, character(1))
    metadata(tse)$metalog <- list(
        source = "https://metalog.embl.de/",
        license = "Open Database License (ODbL) v1.0",
        collection = collection,
        meta.type = meta.type,
        date_profile = date_profile,
        date_metadata = date_metadata,
        date_fetched = Sys.Date()
    )
    return(tse)
}

# Merge per-collection assay matrices into a single matrix.
# Rows (taxa) are the union; columns (samples) are concatenated. Sample
# alias collisions across collections trigger a hard error.
#' @importFrom Matrix sparseMatrix
.merge_metalog_assays <- function(assay_lists) {
    if (length(assay_lists) == 1) {
        al <- assay_lists[[1]]
        return(list(
            assay = al[["assay"]],
            taxa = al[["taxa"]],
            samples = al[["samples"]]
        ))
    }
    taxa <- sort(unique(unlist(
        lapply(assay_lists, `[[`, "taxa"), use.names = FALSE)))
    samples <- unlist(
        lapply(assay_lists, `[[`, "samples"), use.names = FALSE)
    dups <- unique(samples[duplicated(samples)])
    if (length(dups) > 0) {
        stop(
            "Duplicate sample aliases across collections: ",
            paste(dQuote(utils::head(dups, 5)), collapse = ", "),
            if (length(dups) > 5) ", ..." else "",
            call. = FALSE
        )
    }
    # Collect (i, j, x) triplets from each per-collection sparse matrix and
    # remap them into the merged taxa/sample index space.
    sample_offset <- 0L
    triplets <- lapply(assay_lists, function(al) {
        m <- methods::as(al[["assay"]], "TsparseMatrix")
        i <- match(rownames(m)[m@i + 1L], taxa)
        j <- m@j + 1L + sample_offset
        sample_offset <<- sample_offset + ncol(m)
        list(i = i, j = j, x = m@x)
    })
    X <- Matrix::sparseMatrix(
        i = unlist(lapply(triplets, `[[`, "i"), use.names = FALSE),
        j = unlist(lapply(triplets, `[[`, "j"), use.names = FALSE),
        x = unlist(lapply(triplets, `[[`, "x"), use.names = FALSE),
        dims = c(length(taxa), length(samples)),
        dimnames = list(taxa, samples)
    )
    list(assay = X, taxa = taxa, samples = samples)
}

# Merge per-collection metadata data.frames. Adds a 'collection' column
# identifying the source collection of each sample. Non-overlapping
# columns are filled with NA.
.merge_metalog_metadata <- function(md_list, collections, samples,
        add.collection = TRUE) {
    if (add.collection) {
        for (i in seq_along(md_list)) {
            md_list[[i]][["collection"]] <- collections[[i]]
        }
    }
    combined <- data.table::rbindlist(
        lapply(md_list, function(x) {
            x[["sample_alias"]] <- rownames(x)
            x
        }),
        fill = TRUE, use.names = TRUE
    )
    df <- as.data.frame(combined)
    rownames(df) <- df[["sample_alias"]]
    df[["sample_alias"]] <- NULL
    # Reorder to match merged sample order
    df <- df[samples, , drop = FALSE]
    df
}

################################ HELP FUNCTIONS ################################

# Extract the YYYY-MM-DD date stamp from a Metalog filename.
# Returns NA_character_ if no date is found.
.parse_metalog_date <- function(path) {
    fname <- basename(path)
    m <- regmatches(fname, regexpr("[0-9]{4}-[0-9]{2}-[0-9]{2}", fname))
    if (length(m) == 0) NA_character_ else m
}

# Download a file from Metalog, handling HTTP -> HTTPS redirect issues.
# Adapted from Metalog's example script.
#' @importFrom httr2 request req_user_agent req_options req_error
#'   req_timeout req_retry req_progress req_perform resp_status resp_header
.download_metalog_file <- function(
        target_url,
        download_dir = tools::R_user_dir("mia", "cache"),
        use.cache = TRUE) {
    base_filename <- basename(target_url)
    # Check cache
    if (use.cache) {
        pattern <- sub(
            "latest", "[0-9]{4}-[0-9]{2}-[0-9]{2}", base_filename)
        matching_files <- list.files(
            download_dir, pattern = pattern, full.names = TRUE)
        if (length(matching_files) > 0) {
            latest_file <- max(matching_files)
            message("Loaded cached file: ", latest_file)
            return(latest_file)
        }
    }
    if (!dir.exists(download_dir)) {
        dir.create(download_dir, recursive = TRUE, showWarnings = FALSE)
    }
    message("Fetching file from: ", target_url)
    pkg_version <- utils::packageDescription("mia", fields = "Version")
    ua <- paste0("mia/", pkg_version, " (R; Bioconductor)")
    # Metalog server may downgrade HTTPS to HTTP on redirect. Intercept
    # the redirect and force HTTPS.
    initial_resp <- httr2::request(target_url) |>
        httr2::req_user_agent(ua) |>
        httr2::req_options(followlocation = FALSE) |>
        httr2::req_error(is_error = ~ FALSE) |>
        httr2::req_timeout(300) |>
        httr2::req_perform()
    status <- httr2::resp_status(initial_resp)
    if (status >= 300 && status < 400) {
        final_url <- httr2::resp_header(initial_resp, "location")
        if (is.null(final_url)) {
            stop("Redirect response missing Location header.",
                call. = FALSE)
        }
        final_url <- sub("^http://", "https://", final_url)
        message("Intercepted redirect. Forcing HTTPS: ", final_url)
    } else if (status == 200) {
        final_url <- target_url
    } else {
        stop("Initial request failed with status: ", status,
            call. = FALSE)
    }
    # Download the file with retry and timeout
    filename <- basename(final_url)
    destfile <- file.path(download_dir, filename)
    message("Downloading to: ", destfile)
    tryCatch({
        httr2::request(final_url) |>
            httr2::req_user_agent(ua) |>
            httr2::req_timeout(300) |>
            httr2::req_retry(max_tries = 3) |>
            httr2::req_progress() |>
            httr2::req_perform(path = destfile)
    }, error = function(e) {
        if (file.exists(destfile)) file.remove(destfile)
        stop("Error downloading the file: ", conditionMessage(e),
            call. = FALSE)
    })
    return(destfile)
}

# Construct Metalog download URLs and fetch assay + metadata files
.resolve_metalog_url <- function(collection, meta.type, use.cache) {
    cache_dir <- tools::R_user_dir("mia", "cache")
    base_url <- "https://metalog.embl.de/static/download"
    assay_url <- sprintf(
        "%s/profiles/%s_metaphlan4_latest.tsv.gz", base_url, collection)
    md_url <- sprintf(
        "%s/metadata/%s_%s_long_latest.tsv.gz",
        base_url, collection, meta.type)
    assay_file <- .download_metalog_file(
        target_url = assay_url,
        download_dir = cache_dir,
        use.cache = use.cache
    )
    md_file <- if (identical(meta.type, "none")) {
        NA_character_
    } else {
        .download_metalog_file(
            target_url = md_url,
            download_dir = cache_dir,
            use.cache = use.cache
        )
    }
    list(assay = assay_file, md = md_file)
}

# Load MetaPhlAn4 profiles into a sparse matrix (rows = taxa, cols = samples)
#' @importFrom Matrix sparseMatrix
.load_metalog_assay <- function(path, sep = "\t") {
    # data.table NSE bindings
    clade_name <- rel_abund <- sample_alias <- NULL
    dt <- data.table::fread(path, sep = sep)
    data.table::setnames(
        dt, seq_len(3), c("sample_alias", "clade_name", "rel_abund"))
    dt <- dt[startsWith(clade_name, "t__SGB"), ]
    dt[, rel_abund := as.numeric(rel_abund)]
    dt <- dt[!is.na(rel_abund) & rel_abund != 0]
    # Aggregate duplicates
    data.table::setkey(dt, clade_name, sample_alias)
    dt <- dt[, .(rel_abund = sum(rel_abund)),
        by = .(clade_name, sample_alias)]
    taxa <- sort(unique(dt$clade_name))
    samples <- sort(unique(dt$sample_alias))
    X <- Matrix::sparseMatrix(
        i = match(dt$clade_name, taxa),
        j = match(dt$sample_alias, samples),
        x = dt$rel_abund,
        dims = c(length(taxa), length(samples)),
        dimnames = list(taxa, samples)
    )
    list(assay = X, taxa = taxa, samples = samples)
}

# Load Metalog metadata, pivot to wide, and subset to samples present in assay
.load_metalog_metadata <- function(meta_path, samples, sep = "\t") {
    sample_alias <- metadata_item <- NULL
    dt <- data.table::fread(meta_path, sep = sep, na.strings = c("", "NA"))
    wide <- data.table::dcast(
        dt,
        sample_alias ~ metadata_item,
        value.var = "value",
        fill = NA_character_
    )
    missing <- setdiff(samples, wide$sample_alias)
    if (length(missing) > 0) {
        warning(
            length(missing), " sample(s) present in assay data but missing ",
            "from metadata.", call. = FALSE
        )
    }
    # Subset and reorder via data.table key lookup
    data.table::setkey(wide, sample_alias)
    wide <- wide[.(samples)]
    meta_df <- as.data.frame(wide)
    rownames(meta_df) <- meta_df$sample_alias
    meta_df$sample_alias <- NULL
    meta_df
}

# Map SGB clade names to full taxonomic lineages
.construct_metalog_taxmap <- function(database, taxa) {
    clade_name <- lineage <- NULL
    taxmap <- data.table::fread(database, sep = "\t", header = TRUE)
    taxmap <- taxmap[startsWith(clade_name, "t__SGB")]
    taxmap <- taxmap[!duplicated(clade_name)]
    # Align to taxa order; unmatched SGBs get NA
    idx <- match(taxa, taxmap$clade_name)
    n_missing <- sum(is.na(idx))
    if (n_missing > 0) {
        warning(
            n_missing, " of ", length(taxa),
            " taxa could not be mapped to a full lineage.",
            call. = FALSE
        )
    }
    taxmap <- taxmap[idx]
    # Parse lineage into standard taxonomy ranks
    lineage_cols <- c(
        "Kingdom", "Phylum", "Class", "Order",
        "Family", "Genus", "Species", "SGB"
    )
    taxmap[, (lineage_cols) := data.table::tstrsplit(lineage, "\\|")]
    result <- as.data.frame(taxmap[, lineage_cols, with = FALSE])
    rownames(result) <- taxa
    result
}

# Download the MetaPhlAn4 SGB phylogeny and prune to taxa in the dataset.
# Returns a list with the pruned tree and a per-taxon tip label vector
# (NA for taxa absent from the tree) suitable for rowNodeLab.
.construct_metalog_tree <- function(taxa, use.cache) {
    tree_url <- paste0(
        "http://cmprod1.cibio.unitn.it/biobakery4/metaphlan_databases/",
        "mpa_vJun23_CHOCOPhlAnSGB_202307.nwk"
    )
    tree_file <- .download_metalog_file(tree_url, use.cache = use.cache)
    tree <- ape::read.tree(tree_file)
    # Tree tips are bare numeric SGB ids (e.g. "122987"); taxa are
    # MetaPhlAn clade names like "t__SGB1234" or "t__SGB1234_group".
    # Extract the numeric id from each taxon to match tips.
    taxa_id <- vapply(taxa, function(x) {
        m <- regmatches(x, regexpr("SGB[0-9]+", x))
        if (length(m) == 0) NA_character_ else sub("^SGB", "", m)
    }, character(1))
    node_lab <- ifelse(taxa_id %in% tree$tip.label, taxa_id, NA_character_)
    keep <- node_lab[!is.na(node_lab)]
    n_missing <- sum(is.na(node_lab))
    if (n_missing > 0) {
        warning(
            n_missing, " of ", length(taxa),
            " taxa could not be matched to a tip in the SGB tree.",
            call. = FALSE
        )
    }
    if (length(keep) == 0) {
        stop("No taxa matched any tip in the SGB tree.", call. = FALSE)
    }
    tree <- ape::keep.tip(tree, keep)
    list(tree = tree, node_lab = node_lab)
}

# Filter assay data to samples listed in a sample list file
.filter_metalog_samples <- function(assay_list, samplelist) {
    ext <- tolower(tools::file_ext(samplelist))
    if (ext %in% c("csv", "tsv")) {
        sl_df <- data.table::fread(samplelist)
    } else if (ext == "json") {
        .require_package("jsonlite")
        sl_df <- as.data.frame(jsonlite::fromJSON(samplelist))
    }
    target_samples <- sl_df[["sample_alias"]]
    available_samples <- assay_list[["samples"]]
    keep_samples <- intersect(target_samples, available_samples)
    if (length(keep_samples) == 0) {
        stop(
            "None of the samples in 'samplelist' were found in the dataset.",
            call. = FALSE
        )
    }
    # Subset columns (samples)
    assay_list$assay <- assay_list$assay[, keep_samples, drop = FALSE]
    assay_list$samples <- keep_samples
    # Drop taxa with zero abundance
    row_sums <- rowSums(assay_list$assay)
    keep_taxa <- names(row_sums[row_sums > 0])
    assay_list$assay <- assay_list$assay[keep_taxa, , drop = FALSE]
    assay_list$taxa <- keep_taxa
    assay_list
}
