#' Exporters to common formats for microbiome data outside of R
#'
#' @name export-methods
#' @aliases exportRaw exportQIIME2 exportMothur
#'
#' @description
#' There are a few very popular external tools for microbiome analysis,
#' including QIIME2 and mothur. However, R does not currently provide any class
#' to accommodate those data formats. When exporting data from mia to external
#' tools, the best approach is therefore to break a data container into its
#' building blocks (assays, side information, trees, etc.).
#'
#' Thanks to \code{exportRaw}, \code{exportQIIME2} and \code{exportMothur},
#' it is now possible to export a
#' \code{\link[TreeSummarizedExperiment]{TreeSummarizedExperiment}} object as
#' raw elements or near-ready QIIME2 and mothur formats, respectively. This way,
#' migrating from mia to an external system is still a bad idea, but at least it
#' is fairly straightforward.
#'
#' @param x a \code{\link[TreeSummarizedExperiment]{TreeSummarizedExperiment}}
#'   object.
#'
#' @param dpath \code{Character scalar}. String specifying the directory where
#'   \code{x} should be exported. If non-existent, it is created recursively.
#'
#' @param assay.type \code{Character scalar}. String specifying which assay
#'   should be exported. Only for \code{exportQIIME2} and \code{exportMothur}.
#'   (Default: \code{"counts"})
#'
#' @param tree.name \code{Character scalar}. String specifying which row tree
#'   should be exported. Only for \code{exportQIIME2} and \code{exportMothur}.
#'   (Default: \code{"phylo"})
#'
#' @param group.var \code{Character scalar}. String specifying one variable of
#'   rowData to export as a custom grouping file. Only for \code{exportQIIME2}
#'   and \code{exportMothur}. (Default: \code{NULL})
#'
#' @param ... additional arguments:
#'
#' \itemize{
#'   \item \code{rowdata.file}: \code{Character scalar}. Specifies the
#'   file where rowData should be written. (Default: \code{"rowData"})
#'
#'    \item \code{coldata.file}: \code{Character scalar}. Specifies the
#'   file where colData should be written. (Default: \code{"colData"})
#'
#'    \item \code{assay.dir}: \code{Character scalar}. Specifies the
#'   file where assays should be written. (Default: \code{"assays"})
#'
#'    \item \code{refseq.file}: \code{Character scalar}. Specifies the
#'   file where referenceSeq should be written. (Default: \code{"referenceSeq"})
#'
#'    \item \code{rowtree.dir}: \code{Character scalar}. Specifies the
#'   file where rowTrees should be written. (Default: \code{"rowTrees"})
#'
#'   \item \code{coltree.dir}: \code{Character scalar}. Specifies the
#'   file where colTrees should be written. (Default: \code{"colTrees"})
#'
#'   \item \code{dimred.dir}: \code{Character scalar}. Specifies the
#'   file where reducedDims should be written. (Default: \code{"reducedDims"})
#'
#'   \item \code{altexp.dir}: \code{Character scalar}. Specifies the
#'   file where altExps should be written. (Default: \code{"altExps"})
#' }
#'
#' @returns Directory at \code{dpath} with components of \code{x} each stored
#' as a file in the proper format.
#'
#' @details The output directory contains the elements of \code{x}. For
#' \code{exportQIIME2} and \code{exportMothur}, data will need some more
#' processing using the target tool. For some tips, check the rbiom package
#' vignettes on converting data:
#' \url{https://cmmr.github.io/rbiom/articles/convert.html}
#'
#' @examples
#' # Generate example TreeSE object
#' tse <- makeTSE()
#' # Give name to counts assay
#' assayNames(tse) <- "counts"
#' # Give name to taxonomy variable
#' names(rowData(tse))[1] <- "Genus"
#'
#' # Export raw TreeSE components in custom directory
#' exportRaw(tse, "out")
#'
#' # Export TreeSE components in near-ready QIIME2 format
#' exportQIIME2(tse, "qiime2_dir")
#'
#' # Export TreeSE components in near-ready mothur format
#' exportMothur(tse, "mothur_dir")
NULL

#' @rdname export-methods
#' @importFrom ape write.tree write.FASTA
setMethod("exportRaw", signature = c(x = "TreeSummarizedExperiment"),
    function(x, dpath, ...){
    # Create export directory
    dpath <- .create_export_dir(dpath)
    # Write row data
    .write_rowdata_table(x, dpath, ...)
    # Write col data
    .write_coldata_table(x, dpath, ...)
    # Create directory for assays
    .write_assay_tables(x, dpath, ...)
    # Create directory for row trees
    .write_rowtree_data(x, dpath, ...)
    # Create directory for col trees
    .write_coltree_data(x, dpath, ...)
    # Create directory for reduced dimensions
    .write_dimred_table(x, dpath, ...)
    # Create directory for alternative experiments
    .write_altexp_data(x, dpath, ...)
    # Write reference sequence
    .write_referenceSeq(x, dpath, ...)
    invisible(NULL)
})

#' @rdname export-methods
#' @importFrom ape write.tree write.FASTA
setMethod("exportQIIME2", signature = c(x = "TreeSummarizedExperiment"),
    function(x, dpath, assay.type = "counts", tree.name = "phylo",
    group.var = NULL){
    # Check for forbidden colnames
    if( "sample-id" %in% colnames(colData(x)) ){
        warning(
            "'sample-id' is a reserved name in QIIME2. colData variables with ",
            "that name were made unique.", call. = FALSE
        )
    }
    # Create export directory
    dpath <- .create_export_dir(dpath)
    # Write file with custom grouping variable
    .write_group_file(x, dpath, group.var, file.format = "tsv")
    # Concatenate taxonomic ranks feature-wise
    row_data <- .collapse_taxranks(x, ";_")
    # Prepare row data using QIIME2 format
    row_data <- data.frame(rownames(x), row_data, 1L, row.names = NULL)
    # Name variables using QIIME2 format
    colnames(row_data) <- c("Feature ID", "Taxon", "Confidence")
    # Write row data
    write.table(
        row_data, paste0(dpath, "taxonomy.tsv"),
        sep = "\t", quote = FALSE, row.names = FALSE
    )
    # Convert col data to data frame
    col_data <- as.data.frame(colData(x))
    # Convert factor variables to character to correctly record col types
    col_data[] <- lapply(
        col_data, function(col) if( is.factor(col) ) as.character(col) else col
    )
    # Record col types
    col_types <- apply(
        col_data, 2L, function(col) switch(
        type(col), character = "categorical", integer = , double = "numeric")
    )
    # Add col types to first row of col data
    col_data <- rbind(col_types, col_data)
    # Add sample id to first var of col data
    col_data <- cbind(`sample-id` = c("#q2:types", colnames(x)), col_data)
    # Write col data
    write.table(
        col_data, paste0(dpath, "metadata.tsv"),
        sep = "\t", quote = FALSE, row.names = FALSE
    )
    # Add feature names as first variable of assay
    sel_assay <- data.frame(rownames(x), assay(x, assay.type), row.names = NULL)
    # Name first variable according to QIIME2 format
    colnames(sel_assay)[1L] <- "#OTU ID"
    # Write selected assay
    write.table(
        sel_assay, paste0(dpath, assay.type, ".tsv"),
        sep = "\t", quote = FALSE, row.names = FALSE
    )
    # Retrieve row tree
    row_tree <- rowTree(x, tree.name)
    # If row tree is present
    if( !is.null(row_tree) ){
        # Write row tree
        write.tree(row_tree, paste0(dpath, "tree.nwk"))
    }
    # Write reference sequence
    .write_referenceSeq(x, dpath, "seqs")
    invisible(NULL)
})

#' @rdname export-methods
#' @importFrom ape write.tree write.FASTA
#' @importFrom stringr str_replace_all fixed
setMethod("exportMothur", signature = c(x = "TreeSummarizedExperiment"),
    function(x, dpath, assay.type = "counts", tree.name = "phylo",
    group.var = NULL){
    # Check for forbidden colnames
    if( "group" %in% colnames(colData(x)) ){
        warning(
            "'group' is a reserved name in Mothur. colData variables with ",
            "that name were made unique.", call. = FALSE
        )
    }
    # Create export directory
    dpath <- .create_export_dir(dpath)
    # Make format compatible with Mothur
    rownames(x) <- str_replace_all(fixed(rownames(x)), "-", "_")
    colnames(x) <- str_replace_all(fixed(colnames(x)), "-", "_")
    # Write file with custom grouping variable
    .write_group_file(x, dpath, group.var, file.format = "group")
    # Concatenate taxonomic ranks feature-wise
    row_data <- .collapse_taxranks(x, ";")
    # Make format compatible with Mothur
    row_data <- str_replace_all(fixed(row_data), "-", "_")
    # Retrieve selected assay
    sel_assay <- assay(x, assay.type)
    # Calculate assay row sums
    row_sums <- rowSums(sel_assay)
    # Add row sums to row data
    row_data <- data.frame(
        names(row_data), row_sums, row_data, row.names = NULL
    )
    # Name variables of row data according to Mothur format
    colnames(row_data) <- c("OTU", "Size", "Taxonomy")
    # Write row data
    write.table(
        row_data, paste0(dpath, "taxonomy.tsv"),
        sep = "\t", quote = FALSE, row.names = FALSE
    )
    # Add sample ids to first variable of col data
    col_data <- data.frame(group = colnames(x), colData(x), row.names = NULL)
    # Write col data
    write.table(
        col_data, paste0(dpath, "metadata.tsv"),
        sep = "\t", quote = FALSE, row.names = FALSE
    )
    # Add feature names and row sums to selected assay
    sel_assay <- data.frame(
        rownames(sel_assay), total = row_sums, sel_assay, row.names = NULL
    )
    # Name first variable according to Mothur format
    colnames(sel_assay)[1L] <- "Representative_Sequence"
    # Write selected assay
    write.table(
        sel_assay, paste0(dpath, assay.type, ".tsv"),
        sep = "\t", quote = FALSE, row.names = FALSE
    )
    # Retrieve row tree
    row_tree <- rowTree(x, tree.name)
    # If row tree is present
    if( !is.null(row_tree) ){
        # Make format compatible with Mothur
        row_tree$tip.label <- str_replace_all(
            fixed(row_tree$tip.label), "-", "_")
        # Write row tree
        write.tree(row_tree, paste0(dpath, "tree.nwk"))
    }
    # Write reference sequence
    .write_referenceSeq(x, dpath, "seqs")
    invisible(NULL)
})

################################ HELP FUNCTIONS ################################

# Define function to create main export directory
.create_export_dir <- function(dpath){
    # Add final slash to path if absent
    if( !endsWith(dpath, "/") ) dpath <- paste0(dpath, "/")
    # Create export directory
    if( !dir.exists(dpath) ) dir.create(dpath, recursive = TRUE)
    return(dpath)
}

# Define function to create directories for TreeSE slots
.create_slot_dir <- function(x, FUN, slot.path, main.path = ""){
    # If slot contains elements
    if( length(FUN(x)) != 0L ){
        # Build path to slot
        slot.path <- paste0(main.path, slot.path)
        # Add final slash if absent
        if( !endsWith(slot.path, "/") ) slot.path <- paste0(slot.path, "/")
        # Create directory for slot
        dir.create(slot.path)
    }
    return(slot.path)
}

.write_rowdata_table <- function(x, dpath, rowdata.file = "rowData", ...){
    write.table(rowData(x), paste0(dpath, rowdata.file, ".tsv"), sep = "\t")
    invisible(NULL)
}

.write_coldata_table <- function(x, dpath, coldata.file = "colData", ...){
    write.table(colData(x), paste0(dpath, coldata.file, ".tsv"), sep = "\t")
    invisible(NULL)
}

.write_assay_tables <- function(x, dpath, assay.dir = "assays", ...){
    assay_dir <- .create_slot_dir(x, assays, assay.dir, dpath)
    # For each assay
    for( assay_name in assayNames(x) ){
        # Write assay
        write.table(
            assay(x, assay_name),
            paste0(assay_dir, assay_name, ".tsv"),
            sep = "\t"
        )
    }
    invisible(NULL)
}

.write_rowtree_data <- function(x, dpath, rowtree.dir = "rowTrees", ...){
    rowtree_dir <- .create_slot_dir(x, rowTreeNames, rowtree.dir, dpath)
    # For each row tree
    for( tree_name in rowTreeNames(x) ){
        # Write row tree
        write.tree(
            rowTree(x, tree_name), paste0(rowtree_dir, tree_name, ".nwk")
        )
    }
    invisible(NULL)
}

.write_coltree_data <- function(x, dpath, coltree.dir = "colTrees", ...){
    coltree_dir <- .create_slot_dir(x, colTreeNames, coltree.dir, dpath)
    # For each col tree
    for( tree_name in colTreeNames(x) ){
        # Write col tree
        write.tree(
            colTree(x, tree_name), paste0(coltree_dir, tree_name, ".nwk")
        )
    }
    invisible(NULL)
}

.write_dimred_table <- function(x, dpath, dimred.dir = "reducedDims", ...){
    dimred_dir <- .create_slot_dir(x, reducedDims, dimred.dir, dpath)
    # For each reduced dimension
    for( dimred_name in reducedDimNames(x) ){
        # Write reduced dimension
        write.table(
            reducedDim(x, dimred_name),
            paste0(dimred_dir, dimred_name, ".tsv"),
            sep = "\t"
        )
    }
    invisible(NULL)
}

.write_altexp_data <- function(x, dpath, altexp.dir = "altExps", ...){
    altexp_dir <- .create_slot_dir(x, altExps, altexp.dir, dpath)
    # For each alternative experiment
    for( altexp_name in altExpNames(x) ){
        # Build path to experiment
        altexp_path <- paste0(altexp_dir, altexp_name)
        # If instance of SE
        if( is(altExp(x, altexp_name), "SummarizedExperiment") ){
            # Apply exporter recursively
            exportToRaw(altExp(x, altexp_name), altexp_path)
        }else{
            warning("altExps skipped.", call. = FALSE)
        }
    }
    invisible(NULL)
}

.write_referenceSeq <- function(x, dpath, refseq.file = "referenceSeq", ...){
    # If reference sequence is present
    if( !is.null(referenceSeq(x)) ){
        # Write reference sequence
        write.FASTA(referenceSeq(x), paste0(dpath, refseq.file, ".fna"))
    }
    invisible(NULL)
}

# Define function to write file for custom grouping variable
#' @importFrom stringr str_replace_all fixed
.write_group_file <- function(x, dpath, group.var, file.format){
    if( !is.null(group.var) ){
        # Extract grouping variable from row data
        group <- rowData(x)[[group.var]]
        # Make format compatible with QIIME2 and Mothur
        group <- str_replace_all(fixed(group), "-", "_")
        # Add rownames as first variable
        group <- cbind(rownames(x), group)
        # Name variables according to QIIME2 format
        colnames(group) <- c("Feature ID", group.var)
        # Write file with custom grouping variable
        write.table(
            group, paste0(dpath, group.var, ".", file.format),
            sep = "\t", quote = FALSE, row.names = FALSE
        )
    }
}

# Define function to concatenate taxonomic ranks feature-wise
#' @importFrom stringr str_remove
.collapse_taxranks <- function(x, sep){
    # Concatenate taxonomic ranks feature-wise
    row_data <- apply(rowData(x)[taxonomyRanks(x)], 1L, paste, collapse = sep)
    # Remove empty ranks
    pattern <- paste0("(", sep, "|", sep, "NA)+$")
    row_data <- str_remove(row_data, pattern)
}
