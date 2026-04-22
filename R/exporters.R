#' Exporters to common formats for microbiome data outside of R
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
#' @param x a \code{\link[TreeSummarizedExperiment]{TreeSummarizedExperiment}} object.
#' 
#' @param dpath \code{Character scalar}.
#' 
#' @param assay.type \code{Character scalar}. (Default: \code{"counts"})
#' 
#' @param rowdata \code{Character scalar}. (Default: \code{"rowdata"})
#' 
#' @param coldata \code{Character scalar}. (Default: \code{"coldata"})
#' 
#' @param assay.dir \code{Character scalar}. (Default: \code{"assays"})
#' 
#' @param rowtree.dir \code{Character scalar}. (Default: \code{"row_trees"})
#' 
#' @param coltree.dir \code{Character scalar}. (Default: \code{"col_trees"})
#' 
#' @param dimred.dir \code{Character scalar}. (Default: \code{"dim_reds"})
#' 
#' @param altexp.dir \code{Character scalar}. (Default: \code{"alt_exps"})
#' 
#' @param ... Unused.
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
#' library(TreeSummarizedExperiment)
#' 
#' tse <- makeTSE()
#' assayNames(tse) <- "counts"
#' 
#' # Export raw TreeSE components in custom directory
#' exportRaw(tse, "out")
#' 
#' # Export TreeSE components in near-ready QIIME2 format
#' exportQIIME2(tse, "qiime2_dir")
#' 
#' # Export TreeSE components in near-ready mothur format
#' exportMothur(tse, "mothur_dir")
#' 
#' @name export-methods
#' @aliases exportRaw exportQIIME2 exportMothur
NULL


#' @rdname export-methods
#' @importFrom ape write.tree
setMethod("exportRaw", signature = c(x = "TreeSummarizedExperiment"),
    function(x, dpath, rowdata.file = "rowdata", coldata.file = "coldata",
    assay.dir = "assays", rowtree.dir = "row_trees", coltree.dir = "col_trees",
    dimred.dir = "dim_reds", altexp.dir = "alt_exps"){
    
    if( !endsWith(dpath, "/") ) dpath <- paste0(dpath, "/")
    if( !dir.exists(dpath) ) dir.create(dpath)
    
    write.table(rowData(x), paste0(dpath, rowdata.file, ".tsv"), sep = "\t")
    
    write.table(colData(x), paste0(dpath, coldata.file, ".tsv"), sep = "\t")
    
    assay.dir <- .create_slot_dir(x, assays, assay.dir, dpath)
    
    for( assay_name in assayNames(x) ){
        write.table(
            assay(x, assay_name),
            paste0(assay.dir, assay_name, ".tsv"),
            sep = "\t"
        )
    }
    
    rowtree.dir <- .create_slot_dir(x, rowTreeNames, rowtree.dir, dpath)
    
    for( tree_name in rowTreeNames(x) ){
        write.tree(
            rowTree(x, tree_name), paste0(rowtree.dir, tree_name, ".nwk")
        )
    }
    
    coltree.dir <- .create_slot_dir(x, colTreeNames, coltree.dir, dpath)
    
    for( tree_name in colTreeNames(x) ){
        write.tree(
            colTree(x, tree_name), paste0(coltree.dir, tree_name, ".nwk")
        )
    }
    
    dimred.dir <- .create_slot_dir(x, reducedDims, dimred.dir, dpath)
    
    for( dimred_name in reducedDimNames(x) ){
        write.table(
            reducedDim(x, dimred_name),
            paste0(dimred.dir, dimred_name, ".tsv"),
            sep = "\t"
        )
    }
    
    altexp.dir <- .create_slot_dir(x, altExps, altexp.dir, dpath)
    
    for( altexp_name in altExpNames(x) ){
        
        altexp_path <- paste0(altexp.dir, altexp_name)
        
        if( is(altExp(x, altexp_name), "SummarizedExperiment") ){
        
            exportToRaw(altExp(x, altexp_name), altexp_path)
        
        }else{
            
            write.table(altExp(x, altexp_name), altexp_path, sep = "\t")
            
        }
    }
    invisible(NULL)
})


.create_slot_dir <- function(x, FUN, slot.path, main.path = ""){
  
    if( length(FUN(x)) != 0L ){
        slot.path <- paste0(main.path, slot.path)
        if( !endsWith(slot.path, "/") ) slot.path <- paste0(slot.path, "/")
        dir.create(slot.path)
    }
    return(slot.path)
}


#' @rdname export-methods
#' @importFrom ape write.tree write.FASTA
setMethod("exportQIIME2", signature = c(x = "TreeSummarizedExperiment"),
    function(x, dpath, assay.type = "counts", tree.name = "phylo"){
    
    if( !endsWith(dpath, "/") ) dpath <- paste0(dpath, "/")
    if( !dir.exists(dpath) ) dir.create(dpath)
    
    row_data <- apply(rowData(x), 1L, paste, collapse = ";_")
    row_data <- gsub(";_$", "", row_data)
    
    row_data <- cbind(names(row_data), row_data, 1L)
    colnames(row_data) <- c("Feature ID", "Taxon", "Confidence")
    
    write.table(
        row_data, paste0(dpath, "taxonomy.tsv"), sep = "\t", row.names = FALSE
    )
    
    col_data <- as.data.frame(colData(x))
    
    col_types <- apply(
        col_data, 2L, function(col) switch(
        type(col), character = "categorical", integer = , double = "numeric")
    )
    
    col_data <- rbind(col_types, col_data)
    col_data <- cbind(`sample-id` = rownames(col_data), col_data)
    col_data[1L, "sample-id"] <- "#q2:types"
    
    write.table(
        col_data, paste0(dpath, "metadata.tsv"), sep = "\t", row.names = FALSE
    )
    
    sel_assay <- assay(x, assay.type)
    sel_assay <- cbind(`#OTU ID` = rownames(sel_assay), sel_assay)
    
    write.table(
        sel_assay, paste0(dpath, assay.type, ".tsv"),
        sep = "\t", row.names = FALSE
    )
    
    row_tree <- rowTree(x, tree.name)
    
    if( !is.null(row_tree) ){
        write.tree(row_tree, paste0(dpath, "tree.nwk"))
    }
    
    if( !is.null(referenceSeq(x)) ){
        write.FASTA(referenceSeq(x), paste0(dpath, "seqs.fna"))
    }
    invisible(NULL)
})


#' @rdname export-methods
#' @importFrom ape write.tree write.FASTA
setMethod("exportMothur", signature = c(x = "TreeSummarizedExperiment"),
    function(x, dpath, assay.type = "counts", tree.name = "phylo"){
    
    if( !endsWith(dpath, "/") ) dpath <- paste0(dpath, "/")
    if( !dir.exists(dpath) ) dir.create(dpath)
    
    row_data <- apply(rowData(x), 1L, paste, collapse = ";")
    row_data <- gsub(";$", "", row_data)
    
    sel_assay <- assay(x, assay.type)
    row_sums <- rowSums(sel_assay)
    
    row_data <- cbind(names(row_data), row_sums, row_data)
    colnames(row_data) <- c("OTU", "Size", "Taxonomy")
    
    write.table(
        row_data, paste0(dpath, "taxonomy.tsv"), sep = "\t", row.names = FALSE
    )
    
    col_data <- as.data.frame(colData(x))
    col_data <- cbind(group = rownames(col_data), col_data)

    write.table(
        col_data, paste0(dpath, "metadata.tsv"), sep = "\t", row.names = FALSE
    )
    
    sel_assay <- cbind(
        `Representative Sequence` = rownames(sel_assay),
        total = row_sums, sel_assay
    )
    
    write.table(
        sel_assay, paste0(dpath, assay.type, ".tsv"),
        sep = "\t", row.names = FALSE
    )
    
    row_tree <- rowTree(x, tree.name)
    
    if( !is.null(row_tree) ){
        write.tree(row_tree, paste0(dpath, "tree.nwk"))
    }
    
    if( !is.null(referenceSeq(x)) ){
        write.FASTA(referenceSeq(x), paste0(dpath, "seqs.fna"))
    }
    invisible(NULL)
})
