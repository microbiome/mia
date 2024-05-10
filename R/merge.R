.norm_f <- function(i, f, dim.type = c("rows","columns"), na.rm = FALSE, ...){
    dim.type <- match.arg(dim.type)
    if(!.is_a_bool(na.rm)){
      stop("'na.rm' must be TRUE or FALSE.", call. = FALSE)
    }
    if(!is.character(f) && !is.factor(f)){
        stop("'f' must be a factor or character vector coercible to a ",
            "meaningful factor.",
            call. = FALSE)
    }
    if(i != length(f)){
        stop("'f' must have the same number of ",dim.type," as 'x'",
            call. = FALSE)
    }
    # This is done otherwise we lose NA values
    if ( na.rm ) {
      f[is.na(f)] <- "NA"
    }
    if(is.character(f)){
        f <- factor(f)
    }
    f
}

.norm_archetype <- function(f, archetype){
    if(length(archetype) > 1L){
        if(length(levels(f)) != length(archetype)){
            stop("length of 'archetype' must have the same length as ",
                "levels('f')",
                call. = FALSE)
        }
    }
    f_table <- table(f)
    if(!is.null(names(archetype))){
        if(anyNA(names(archetype)) || anyDuplicated(names(archetype))){
            stop("If 'archetype' is named, names must be non-NA and unqiue.",
                call. = FALSE)
        }
        archetype <- archetype[names(f_table)]
    }
    if(any(f_table < archetype)){
        stop("'archetype' out of bounds for some levels of 'f'. The maximum of",
            " 'archetype' is defined as table('f')", call. = FALSE)
    }
    if(length(archetype) == 1L){
        archetype <- rep(archetype,length(levels(f)))
    }
    archetype
}

#' @importFrom S4Vectors splitAsList
.get_element_pos <- function(f, archetype){
    archetype <- as.list(archetype)
    f_pos <- seq.int(1L, length(f))
    f_pos_split <- S4Vectors::splitAsList(f_pos, f)
    f_pos <- unlist(f_pos_split[archetype])
    f_pos
}

#' @importFrom S4Vectors SimpleList
#' @importFrom scuttle sumCountsAcrossFeatures
.merge_rows <- function(x, f, archetype = 1L,
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
    if( .is_a_string(f) && f %in% colnames(rowData(x)) ){
        f <- rowData(x)[[ f ]]
    }
    f <- .norm_f(nrow(x), f)
    if(length(levels(f)) == nrow(x)){
        return(x)
    }

    archetype <- .norm_archetype(f, archetype)
    # merge assays
    assays <- assays(x)
    if( check.assays ){
        mapply(.check_assays_for_merge, names(assays), assays)
    }
    assays <- S4Vectors::SimpleList(lapply(assays,
                                            scuttle::sumCountsAcrossFeatures,
                                            ids = f,
                                            subset.row = NULL,
                                            subset.col = NULL,
                                            average = average,
                                            BPPARAM = BPPARAM))
    names(assays) <- names(assays(x))
    # merge to result
    x <- x[.get_element_pos(f, archetype = archetype),]
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
.merge_cols <- function(x, f, archetype = 1L, ...){
    # input check
    if( .is_a_string(f) && f %in% colnames(rowData(x)) ){
      f <- rowData(x)[[ f ]]
    }
    f <- .norm_f(ncol(x), f, "columns")
    
    if(length(levels(f)) == ncol(x)){
        return(x)
    }
    archetype <- .norm_archetype(f, archetype)
    # merge col data
    element_pos <- .get_element_pos(f, archetype = archetype)
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
                                            ids = f,
                                            subset.row = NULL,
                                            subset.col = NULL,
                                            ...))
    names(assays) <- names(assays(x))
    # merge to result
    x <- x[,.get_element_pos(f, archetype = archetype)]
    assays(x, withDimnames = FALSE) <- assays
    # Change colnames to group names
    colnames(x) <- colnames(assays[[1]])
    x
}

.merge_rows_SE <- function(x, f, archetype = 1L, ...){
    .merge_rows(x, f, archetype = archetype, ...)
}

.merge_cols_SE <- function(x, f, archetype = 1L, ...){
    .merge_cols(x, f, archetype = archetype, ...)
}

#' @importFrom Biostrings DNAStringSetList
.merge_refseq_list <- function(sequences_list, f, names, ...){
    threshold <- list(...)[["threshold"]]
    if(is.null(threshold)){
        threshold <- 0.05
    }
    if(!is(sequences_list,"DNAStringSetList")){
        return(.merge_refseq(sequences_list, f, names, threshold))
    }
    names <- names(sequences_list)
    seqs <- DNAStringSetList(lapply(sequences_list, .merge_refseq, f, names,
                                    threshold))
    names(seqs) <- names
    seqs
}

#' @importFrom Biostrings DNAStringSetList
#' @importFrom DECIPHER ConsensusSequence
.merge_refseq <- function(sequences, f, names, threshold){
    sequences <- split(sequences,f)
    seq <- unlist(DNAStringSetList(lapply(sequences, ConsensusSequence,
                                            threshold = threshold)))
    seq
}

.merge_rows_TSE <- function(x, f, archetype = 1L, mergeTree = FALSE,
                            mergeRefSeq = FALSE, ...){
    # input check
    if(!.is_a_bool(mergeTree)){
        stop("'mergeTree' must be TRUE or FALSE.", call. = FALSE)
    }
    if(!.is_a_bool(mergeRefSeq)){
        stop("'mergeRefSeq' must be TRUE or FALSE.", call. = FALSE)
    }
    # for optionally merging referenceSeq
    refSeq <- NULL
    if(mergeRefSeq){
        refSeq <- referenceSeq(x)
    }
    #
    x <- .merge_rows_SE(x, f, archetype = 1L, ...)
    # optionally merge rowTree
    if( mergeTree ){
        x <- .agglomerate_trees(x, 1)
    }
    # optionally merge referenceSeq
    if(!is.null(refSeq)){
        referenceSeq(x) <- .merge_refseq_list(refSeq, f, rownames(x), ...)
    }
    x
}

.merge_cols_TSE <- function(x, f, archetype = 1L, mergeTree = FALSE, ...){
    # input check
    if(!.is_a_bool(mergeTree)){
        stop("'mergeTree' must be TRUE or FALSE.", call. = FALSE)
    }
    #
    x <- .merge_cols_SE(x, f, archetype = 1L, ...)
    # optionally merge colTree
    if( mergeTree ){
        x <- .agglomerate_trees(x, 2)
    }
    return(x)
}
