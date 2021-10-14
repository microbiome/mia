#' Merge a subset of the rows or columns of a \code{SummarizedExperiment}
#'
#' \code{mergeRows}/\code{mergeCols} merge data on rows or columns of a
#' \code{SummarizedExperiment} as defined by a \code{factor} alongside the
#' chosen dimension. Metadata from the \code{rowData} or \code{colData} are
#' retained as defined by \code{archetype}.
#' 
#' \code{\link[SummarizedExperiment:SummarizedExperiment-class]{assay}} are 
#' agglomerated, i.e.. summed up. Other than counts / absolute values might lead
#' to meaningless values. 
#'
#' @param x a \code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}} or
#'   a \code{\link[TreeSummarizedExperiment:TreeSummarizedExperiment-class]{TreeSummarizedExperiment}}
#'
#' @param f A factor for merging. Must be the same length as
#'   \code{nrow(x)/ncol(x)}. Rows/Cols corresponding to the same level will be
#'   merged. If \code{length(levels(f)) == nrow(x)/ncol(x)}, \code{x} will be
#'   returned unchanged.
#'
#' @param archetype Of each level of \code{f}, which element should be regarded
#'   as the archetype and metadata in the columns or rows kept, while merging?
#'   This can be single integer value or an integer vector of the same length
#'   as \code{levels(f)}. (Default: \code{archetype = 1L}, which means the first
#'   element encountered per factor level will be kept)
#'   
#' @param mergeTree \code{TRUE} or \code{FALSE}: should to
#'   \code{rowTree()} also be merged? (Default: \code{mergeTree = FALSE})
#'
#' @param mergeRefSeq \code{TRUE} or \code{FALSE}: should a consensus sequence
#'   calculate? If set to \code{FALSE}, the result from \code{archetype} is
#'   returned; If set to \code{TRUE} the result from
#'   \code{\link[DECIPHER:ConsensusSequence]{DECIPHER::ConsensusSequence}} is
#'   returned. (Default: \code{mergeRefSeq = FALSE})
#'
#' @param ... optional arguments:
#' \itemize{
#'   \item{passed onto \code{\link[scuttle:sumCountsAcrossFeatures]{sumCountsAcrossFeatures}}, except \code{subset_row}, \code{subset_col}}
#' }
#'
#' @details
#' These functions are similar to
#' \code{\link[scuttle:sumCountsAcrossFeatures]{sumCountsAcrossFeatures}}.
#' However, additional support for \code{TreeSummarizedExperiment} was added and
#' science field agnostic names were used. In addition the \code{archetype}
#' argument lets the user select how to preserve row or column data.
#'
#' For merge data of assays the function from \code{scuttle} are used.
#'
#' @name merge-methods
#' @aliases mergeRows mergeCols
#'
#' @return an object with the same class \code{x} with the specified entries
#'   merged into one entry in all relevant components.
#'
#' @seealso
#' \code{\link[scuttle:sumCountsAcrossFeatures]{sumCountsAcrossFeatures}}
#'
#' @examples
#' data(esophagus)
#' esophagus
#' plot(rowTree(esophagus))
#' # get a factor for merging
#' f <- factor(regmatches(rownames(esophagus),
#'                        regexpr("^[0-9]*_[0-9]*",rownames(esophagus))))
#' merged <- mergeRows(esophagus,f)
#' plot(rowTree(merged))
#' #
#' data(GlobalPatterns)
#' GlobalPatterns
#' merged <- mergeCols(GlobalPatterns,colData(GlobalPatterns)$SampleType)
#' merged
NULL

#' @rdname merge-methods
#' @export
setGeneric("mergeRows",
           signature = "x",
           function(x, f, archetype = 1L, ...)
               standardGeneric("mergeRows"))

#' @rdname merge-methods
#' @export
setGeneric("mergeCols",
           signature = "x",
           function(x, f, archetype = 1L, ...)
               standardGeneric("mergeCols"))

.norm_f <- function(i, f, dim_type = c("rows","columns")){
    dim_type <- match.arg(dim_type)
    if(!is.character(f) && !is.factor(f)){
        stop("'f' must be a factor or character vector coercible to a ",
             "meaningful factor.",
             call. = FALSE)
    }
    if(i != length(f)){
        stop("'f' must have the same number of ",dim_type," as 'x'",
             call. = FALSE)
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
.merge_rows <- function(x, f, archetype = 1L, ...){
    # input check
    f <- .norm_f(nrow(x), f)
    if(length(levels(f)) == nrow(x)){
        return(x)
    }
    archetype <- .norm_archetype(f, archetype)
    # merge assays
    assays <- assays(x)
    assays <- S4Vectors::SimpleList(lapply(names(assays), .calculate_merge_row_values, 
                                           x = x, f = f, ...))
    names(assays) <- names(assays(x))
    # merge to result
    x <- x[.get_element_pos(f, archetype = archetype),]
    assays(x, withDimnames = FALSE) <- assays
    x
}

#' @importFrom scuttle sumCountsAcrossFeatures
.calculate_merge_row_values <- function(assay_name, x, f, ...){
    assay <- assay(x, assay_name)
    # Check if assays include binary or negative values
    if( all(assay == 0 | assay == 1) ){
        warning(paste0("'",assay_name,"'", " includes binary values."),
        "\nAgglomeration of it might lead to meaningless values.", 
        "\nCheck the assay, and consider doing transformation again manually", 
        " with agglomerated data.",
                call. = FALSE)
    }
    if( any(assay < 0) ){
        warning(paste0("'",assay_name,"'", " includes negative values."),
                "\nAgglomeration of it might lead to meaningless values.",
                "\nCheck the assay, and consider doing transformation again manually", 
                " with agglomerated data.",
                call. = FALSE)
    }
    assay <- scuttle::sumCountsAcrossFeatures(assay, 
                                              ids = f, 
                                              subset_row = NULL, 
                                              subset_col = NULL, 
                                              ...)
    return(assay)
}

#' @importFrom S4Vectors SimpleList
#' @importFrom scuttle summarizeAssayByGroup
.merge_cols <- function(x, f, archetype = 1L, ...){
    # input check
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
    assays <- S4Vectors::SimpleList(lapply(assays, FUN = function(mat, ...){
        temp <- scuttle::summarizeAssayByGroup(mat,
                                               ids = f,
                                               subset.row = NULL,
                                               subset.col = NULL,
                                               ...)
        # "sum" includes agglomerated (summed up) data
        mat <- assay(temp, "sum")
        return(mat)
    }))
    names(assays) <- names(assays(x))
    # merge to result
    x <- x[,.get_element_pos(f, archetype = archetype)]
    assays(x, withDimnames = FALSE) <- assays
    x
}

#' @rdname merge-methods
#' @export
setMethod("mergeRows", signature = c(x = "SummarizedExperiment"),
    function(x, f, archetype = 1L, ...){
        .merge_rows(x, f, archetype = archetype, ...)
    }
)

#' @rdname merge-methods
#' @export
setMethod("mergeCols", signature = c(x = "SummarizedExperiment"),
    function(x, f, archetype = 1L, ...){
        .merge_cols(x, f, archetype = archetype, ...)
    }
)

.merge_tree <- function(tree, links){
    tips <- sort(setdiff(tree$edge[, 2], tree$edge[, 1]))
    drop_tip <- tips[!(tips %in% unique(links$nodeNum[links$isLeaf]))]
    oldTree <- tree
    newTree <- ape::drop.tip(oldTree, tip = drop_tip)
    track <- trackNode(oldTree)
    track <- ape::drop.tip(track, tip = drop_tip)
    #
    oldAlias <- links$nodeLab_alias
    newNode <- convertNode(tree = track, node = oldAlias)
    newAlias <- convertNode(tree = newTree, node = newNode)
    #
    list(newTree = newTree, newAlias = newAlias)
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

#' @rdname merge-methods
#' @importFrom ape keep.tip
#' @export
setMethod("mergeRows", signature = c(x = "TreeSummarizedExperiment"),
    function(x, f, archetype = 1L, mergeTree = FALSE, mergeRefSeq = FALSE, ...){
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
        x <- callNextMethod(x, f, archetype = 1L, ...)
        # optionally merge rowTree
        tree <- rowTree(x)
        if(!is.null(tree) && mergeTree){
            tmp <- .merge_tree(tree, rowLinks(x))
            #
            x <- changeTree(x = x,
                            rowTree = tmp$newTree,
                            rowNodeLab = tmp$newAlias)
        }
        # optionally merge referenceSeq
        if(!is.null(refSeq)){
            referenceSeq(x) <- .merge_refseq_list(refSeq, f, rownames(x), ...)
        }
        x
    }
)

#' @rdname merge-methods
#' @importFrom ape keep.tip
#' @export
setMethod("mergeCols", signature = c(x = "TreeSummarizedExperiment"),
    function(x, f, archetype = 1L, mergeTree = FALSE, ...){
        # input check
        if(!.is_a_bool(mergeTree)){
            stop("'mergeTree' must be TRUE or FALSE.", call. = FALSE)
        }
        #
        x <- callNextMethod(x, f, archetype = 1L, ...)
        # optionally merge colTree
        tree <- colTree(x)
        if(!is.null(tree) && mergeTree){
            tmp <- .merge_tree(tree, colLinks(x))
            #
            x <- changeTree(x = x,
                            colTree = tmp$newTree,
                            colNodeLab = tmp$newAlias)
        }
        x
    }
)
