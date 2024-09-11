.norm_f <- function(i, group, dim.type = c("rows","columns"), na.rm = FALSE, ...){
    if(!.is_a_bool(na.rm)){
        stop("'na.rm' must be TRUE or FALSE.", call. = FALSE)
    }
    dim.type <- match.arg(dim.type)
    if(!is.character(group) && !is.factor(group)){
        stop("'group' must be a factor or character vector coercible to a ",
            "meaningful factor.",
            call. = FALSE)
    }
    if(i != length(group)){
        stop("'group' must have the same number of ",dim.type," as 'x'",
            call. = FALSE)
    }
    # This is done otherwise we lose NA values
    if( !na.rm && any(is.na(group)) ){
        group <- as.character(group)
        group[ is.na(group) ] <- "NA"
    }
    if(is.character(group)){
        group <- factor(group)
    }
    group
}

.norm_archetype <- function(group, archetype){
    if(length(archetype) > 1L){
        if(length(levels(group)) != length(archetype)){
            stop("length of 'archetype' must have the same length as ",
                "levels('group')",
                call. = FALSE)
        }
    }
    f_table <- table(group)
    if(!is.null(names(archetype))){
        if(anyNA(names(archetype)) || anyDuplicated(names(archetype))){
            stop("If 'archetype' is named, names must be non-NA and unqiue.",
                call. = FALSE)
        }
        archetype <- archetype[names(f_table)]
    }
    if(any(f_table < archetype)){
        stop("'archetype' out of bounds for some levels of 'group'. The maximum of",
            " 'archetype' is defined as table('group')", call. = FALSE)
    }
    if(length(archetype) == 1L){
        archetype <- rep(archetype,length(levels(group)))
    }
    archetype
}

#' @importFrom S4Vectors splitAsList
.get_element_pos <- function(group, archetype){
    archetype <- as.list(archetype)
    f_pos <- seq.int(1L, length(group))
    f_pos_split <- S4Vectors::splitAsList(f_pos, group)
    f_pos <- unlist(f_pos_split[archetype])
    f_pos
}

#' @importFrom S4Vectors SimpleList
#' @importFrom scuttle sumCountsAcrossFeatures
.merge_rows <- function(x, group, archetype = 1L,
                        average = FALSE,
                        BPPARAM = SerialParam(),
                        check.assays = TRUE,
                        ...){
    # input check
    if( !.is_a_bool(average) ){
        stop("'average' must be TRUE or FALSE.", call. = FALSE)
    }
    if( !.is_a_bool(check.assays) ){
        stop("'check.assays' must be TRUE or FALSE.", call. = FALSE)
    }
    if( .is_a_string(group) && group %in% colnames(rowData(x)) ){
        group <- rowData(x)[[ group ]]
    }
    group <- .norm_f(nrow(x), group, ...)
    if(length(levels(group)) == nrow(x)){
        return(x)
    }

    archetype <- .norm_archetype(group, archetype)
    # merge assays
    assays <- assays(x)
    if( check.assays ){
        mapply(.check_assays_for_merge, names(assays), assays)
    }
    assays <- S4Vectors::SimpleList(lapply(assays,
                                            scuttle::sumCountsAcrossFeatures,
                                            ids = group,
                                            subset.row = NULL,
                                            subset.col = NULL,
                                            average = average,
                                            BPPARAM = BPPARAM))
    names(assays) <- names(assays(x))
    # merge to result
    x <- x[.get_element_pos(group, archetype = archetype),]
    assays(x, withDimnames = FALSE) <- assays
    # Change rownames to group names
    rownames(x) <- rownames(assays[[1]])
    x
}

#' @importFrom scuttle sumCountsAcrossFeatures
.check_assays_for_merge <- function(assay.type, assay){
    # Check if assays include binary or negative values
    if( all(assay == 0 | assay == 1) ){
        warning("'",assay.type,"'", " includes binary values.",
                "\nAgglomeration of it might lead to meaningless values.",
                "\nCheck the assay, and consider doing transformation again",
                "manually with agglomerated data.",
                call. = FALSE)
    }
    if( !all( assay >= 0 | is.na(assay) ) ){
        warning("'",assay.type,"'", " includes negative values.",
                "\nAgglomeration of it might lead to meaningless values.",
                "\nCheck the assay, and consider doing transformation again",
                "manually with agglomerated data.",
                call. = FALSE)
    }
}

#' @importFrom S4Vectors SimpleList
#' @importFrom scuttle summarizeAssayByGroup
.merge_cols <- function(x, group, archetype = 1L, ...){
    # input check
    if( .is_a_string(group) && group %in% colnames(colData(x)) ){
      group <- colData(x)[[ group ]]
    }
    group <- .norm_f(ncol(x), group, "columns", ...)
    
    if(length(levels(group)) == ncol(x)){
        return(x)
    }
    archetype <- .norm_archetype(group, archetype)
    # merge col data
    element_pos <- .get_element_pos(group, archetype = archetype)
    col_data <- colData(x)[element_pos,,drop=FALSE]
    # merge assays
    assays <- assays(x)
    mapply(.check_assays_for_merge, names(assays), assays)
    FUN <- function(mat, ...){
        temp <- scuttle::summarizeAssayByGroup(mat,
                                                statistics = "sum",
                                                ...)
        # "sum" includes agglomerated (summed up) data
        mat <- assay(temp, "sum")
        return(mat)
    }
    assays <- S4Vectors::SimpleList(lapply(assays,
                                            FUN = FUN,
                                            ids = group,
                                            subset.row = NULL,
                                            subset.col = NULL,
                                            ...))
    names(assays) <- names(assays(x))
    # merge to result
    x <- x[,.get_element_pos(group, archetype = archetype)]
    assays(x, withDimnames = FALSE) <- assays
    # Change colnames to group names
    colnames(x) <- colnames(assays[[1]])
    x
}

#' @importFrom Biostrings DNAStringSetList
.merge_refseq_list <- function(sequences_list, group, names, ...){
    threshold <- list(...)[["threshold"]]
    if(is.null(threshold)){
        threshold <- 0.05
    }
    if(!is(sequences_list,"DNAStringSetList")){
        return(.merge_refseq(sequences_list, group, names, threshold))
    }
    names <- names(sequences_list)
    seqs <- DNAStringSetList(lapply(sequences_list, .merge_refseq, group, names,
                                    threshold))
    names(seqs) <- names
    seqs
}

#' @importFrom Biostrings DNAStringSetList
#' @importFrom DECIPHER ConsensusSequence
.merge_refseq <- function(sequences, group, names, threshold){
    sequences <- split(sequences,group)
    seq <- unlist(DNAStringSetList(lapply(sequences, ConsensusSequence,
                                            threshold = threshold)))
    seq
}

.merge_rows_TSE <- function(x, group, archetype = 1L, update.tree = FALSE,
    update.refseq = mergeRefSeq, mergeRefSeq = FALSE, ...){
    # input check
    if(!.is_a_bool(update.tree)){
        stop("'update.tree' must be TRUE or FALSE.", call. = FALSE)
    }
    if(!.is_a_bool(update.refseq)){
        stop("'update.refseq' must be TRUE or FALSE.", call. = FALSE)
    }
    # for optionally merging referenceSeq
    refSeq <- NULL
    if(update.refseq){
        refSeq <- referenceSeq(x)
    }
    #
    x <- .merge_rows(x, group, archetype = 1L, ...)
    # optionally merge rowTree
    if( update.tree ){
        x <- .agglomerate_trees(x, 1, ...)
    }
    # optionally merge referenceSeq
    if(!is.null(refSeq)){
        referenceSeq(x) <- .merge_refseq_list(refSeq, group, rownames(x), ...)
    }
    x
}

.merge_cols_TSE <- function(x, group, archetype = 1L, update.tree = FALSE, ...){
    # input check
    if(!.is_a_bool(update.tree)){
        stop("'update.tree' must be TRUE or FALSE.", call. = FALSE)
    }
    #
    x <- .merge_cols(x, group, archetype = 1L, ...)
    # optionally merge colTree
    if( update.tree ){
        x <- .agglomerate_trees(x, 2, ...)
    }
    return(x)
}
