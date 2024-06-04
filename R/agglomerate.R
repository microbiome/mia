#' Agglomerate or merge data using taxonomic information
#'
#' Agglomeration functions can be used to sum-up data based on specific criteria
#' such as taxonomic ranks, variables or prevalence.
#'
#' \code{agglomerateByRank} can be used to sum up data based on associations
#' with certain taxonomic ranks, as defined in \code{rowData}. Only available
#' \code{\link{taxonomyRanks}} can be used.
#'
#' \code{agglomerateByVariable} merges data on rows or columns of a
#' \code{SummarizedExperiment} as defined by a \code{factor} alongside the
#' chosen dimension. This function allows agglomeration of data based on other
#' variables than taxonomy ranks.
#' Metadata from the \code{rowData} or \code{colData} are
#' retained as defined by \code{archetype}.
#' \code{\link[SummarizedExperiment:SummarizedExperiment-class]{assay}} are
#' agglomerated, i.e. summed up. If the assay contains values other than counts
#' or absolute values, this can lead to meaningless values being produced.
#'
#' @param x a \code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}} or
#'   a \code{\link[TreeSummarizedExperiment:TreeSummarizedExperiment-class]{TreeSummarizedExperiment}}
#'
#' @param rank a single character defining a taxonomic rank. Must be a value of
#'   \code{taxonomyRanks()} function.
#'
#' @param onRankOnly \code{TRUE} or \code{FALSE}: Should information only from
#'   the specified rank be used or from ranks equal and above? See details.
#'   (default: \code{onRankOnly = FALSE})
#'
#' @param na.rm \code{TRUE} or \code{FALSE}: Should taxa with an empty rank be
#'   removed? Use it with caution, since empty entries on the selected rank
#'   will be dropped. This setting can be tweaked by defining
#'   \code{empty.fields} to your needs. (default: \code{na.rm = TRUE})
#'
#' @param empty.fields a \code{character} value defining, which values should be
#'   regarded as empty. (Default: \code{c(NA, "", " ", "\t")}). They will be
#'   removed if \code{na.rm = TRUE} before agglomeration.
#'
#' @param agglomerate.tree \code{TRUE} or \code{FALSE}: should
#'   \code{rowTree()} also be agglomerated? (Default:
#'   \code{agglomerate.tree = FALSE})
#'
#' @param agglomerateTree alias for \code{agglomerate.tree}.
#'
#' @param ... arguments passed to \code{agglomerateByRank} function for
#'   \code{SummarizedExperiment} objects,
#'   to \code{\link[=agglomerate-methods]{agglomerateByVariable}} and
#'   \code{\link[scuttle:sumCountsAcrossFeatures]{sumCountsAcrossFeatures}},
#'   to \code{getPrevalence} and \code{getPrevalentTaxa} and used in
#'   \code{agglomeratebyPrevalence}
#'   \itemize{
#'        \item \code{remove_empty_ranks}: A single boolean value for selecting 
#'        whether to remove those columns of rowData that include only NAs after
#'        agglomeration. (By default: \code{remove_empty_ranks = FALSE})
#'        \item \code{make_unique}: A single boolean value for selecting 
#'        whether to make rownames unique. (By default: \code{make_unique = TRUE})
#'        \item \code{detection}: Detection threshold for absence/presence. 
#'        Either an absolute value compared directly to the values of \code{x} 
#'        or a relative value between 0 and 1, if \code{as_relative = FALSE}.
#'        \item \code{prevalence}: Prevalence threshold (in 0 to 1). The 
#'        required prevalence is strictly greater by default. To include the 
#'        limit, set \code{include_lowest} to \code{TRUE}.
#'        \item \code{as.relative}: Logical scalar: Should the detection 
#'        threshold be applied on compositional (relative) abundances? 
#'        (default: \code{FALSE})
#'    }
#'
#' @param altexp String or integer scalar specifying an alternative experiment
#'   containing the input data.
#'
#' @param strip_altexp \code{TRUE} or \code{FALSE}: Should alternative
#'   experiments be removed prior to agglomeration? This prevents to many
#'   nested alternative experiments by default (default:
#'   \code{strip_altexp = TRUE})
#'
#' @param MARGIN A character value for selecting if data is merged
#'   row-wise / for features ('rows') or column-wise / for samples ('cols').
#'   Must be \code{'rows'} or \code{'cols'}.
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
#' @param mergeTree \code{TRUE} or \code{FALSE}: Should
#'   \code{rowTree()} also be merged? (Default: \code{mergeTree = FALSE})
#'
#' @param mergeRefSeq \code{TRUE} or \code{FALSE}: Should a consensus sequence
#'   be calculated? If set to \code{FALSE}, the result from \code{archetype} is
#'   returned; If set to \code{TRUE} the result from
#'   \code{\link[DECIPHER:ConsensusSequence]{DECIPHER::ConsensusSequence}} is
#'   returned. (Default: \code{mergeRefSeq = FALSE})
#'
#' @details
#' When using \code{agglomerateByRank}, please note that depending on the 
#' available taxonomic data and its structure, setting\code{onRankOnly = TRUE} 
#' has certain implications on the interpretability of your results. If no loops
#' exist (loops meaning two higher ranks containing the same lower rank), the 
#' results should be comparable. You can check for loops using 
#' \code{\link[TreeSummarizedExperiment:detectLoop]{detectLoop}}.
#' 
#' Also, agglomeration sums up the values of assays at the specified taxonomic level. With
#' certain assays, e.g. those that include binary or negative values, this summing
#' can produce meaningless values. In those cases, consider performing agglomeration
#' first, and then applying the transformation afterwards.
#' 
#' \code{agglomerateByVariable} works similarly to
#' \code{\link[scuttle:sumCountsAcrossFeatures]{sumCountsAcrossFeatures}}.
#' However, additional support for \code{TreeSummarizedExperiment} was added and
#' science field agnostic names were used. In addition the \code{archetype}
#' argument lets the user select how to preserve row or column data.
#'
#' For merge data of assays the function from \code{scuttle} are used.
#'
#' @return 
#' \code{agglomerateByRank} returns a taxonomically-agglomerated, 
#' optionally-pruned object of the same class as \code{x}.
#' \code{agglomerateByVariable} returns an object of the same class as \code{x} 
#' with the specified entries merged into one entry in all relevant components.
#' \code{agglomerateByRank} returns a taxonomically-agglomerated, 
#' optionally-pruned object of the same class as \code{x}.
#'
#' @name agglomerate-methods
#'
#' @seealso
#' \code{\link[=splitOn]{splitOn}}
#' \code{\link[=unsplitOn]{unsplitOn}}
#' \code{\link[=agglomerate-methods]{agglomerateByVariable}},
#' \code{\link[scuttle:sumCountsAcrossFeatures]{sumCountsAcrossFeatures}},
#' \code{\link[=agglomerate-methods]{agglomerateByRank}},
#' \code{\link[SingleCellExperiment:altExps]{altExps}},
#' \code{\link[SingleCellExperiment:splitAltExps]{splitAltExps}}
#'
#' @examples
#'
#' ### Agglomerate data based on taxonomic information
#'
#' data(GlobalPatterns)
#' # print the available taxonomic ranks
#' colnames(rowData(GlobalPatterns))
#' taxonomyRanks(GlobalPatterns)
#'
#' # agglomerate at the Family taxonomic rank
#' x1 <- agglomerateByRank(GlobalPatterns, rank="Family")
#' ## How many taxa before/after agglomeration?
#' nrow(GlobalPatterns)
#' nrow(x1)
#'
#' # agglomerate the tree as well
#' x2 <- agglomerateByRank(GlobalPatterns, rank="Family",
#'                        agglomerate.tree = TRUE)
#' nrow(x2) # same number of rows, but
#' rowTree(x1) # ... different
#' rowTree(x2) # ... tree
#'
#' # If assay contains binary or negative values, summing might lead to
#' # meaningless values, and you will get a warning. In these cases, you might
#' # want to do agglomeration again at chosen taxonomic level.
#' tse <- transformAssay(GlobalPatterns, method = "pa")
#' tse <- agglomerateByRank(tse, rank = "Genus")
#' tse <- transformAssay(tse, method = "pa")
#'
#' # removing empty labels by setting na.rm = TRUE
#' sum(is.na(rowData(GlobalPatterns)$Family))
#' x3 <- agglomerateByRank(GlobalPatterns, rank="Family", na.rm = TRUE)
#' nrow(x3) # different from x2
#'
#' # Because all the rownames are from the same rank, rownames do not include
#' # prefixes, in this case "Family:".
#' print(rownames(x3[1:3,]))
#'
#' # To add them, use getTaxonomyLabels function.
#' rownames(x3) <- getTaxonomyLabels(x3, with_rank = TRUE)
#' print(rownames(x3[1:3,]))
#'
#' # use 'remove_empty_ranks' to remove columns that include only NAs
#' x4 <- agglomerateByRank(GlobalPatterns, rank="Phylum",
#'                         remove_empty_ranks = TRUE)
#' head(rowData(x4))
#'
#' # If the assay contains NAs, you might want to consider replacing them,
#' # since summing-up NAs lead to NA
#' x5 <- GlobalPatterns
#' # Replace first value with NA
#' assay(x5)[1,1] <- NA
#' x6 <- agglomerateByRank(x5, "Kingdom")
#' head( assay(x6) )
#' # Replace NAs with 0. This is justified when we are summing-up counts.
#' assay(x5)[ is.na(assay(x5)) ] <- 0
#' x6 <- agglomerateByRank(x5, "Kingdom")
#' head( assay(x6) )
#'
#' ## Look at enterotype dataset...
#' data(enterotype)
#' ## Print the available taxonomic ranks. Shows only 1 available rank,
#' ## not useful for agglomerateByRank
#' taxonomyRanks(enterotype)
#'
#' ### Merge TreeSummarizedExperiments on rows and columns
#'
#' data(esophagus)
#' esophagus
#' plot(rowTree(esophagus))
#' # get a factor for merging
#' f <- factor(regmatches(rownames(esophagus),
#'                        regexpr("^[0-9]*_[0-9]*",rownames(esophagus))))
#' merged <- agglomerateByVariable(esophagus, MARGIN = "rows", f,
#'                                 mergeTree = TRUE)
#' plot(rowTree(merged))
#' #
#' data(GlobalPatterns)
#' GlobalPatterns
#' merged <- agglomerateByVariable(GlobalPatterns, MARGIN = "cols",
#'                                 colData(GlobalPatterns)$SampleType)
#' merged
NULL

#' @rdname agglomerate-methods
#' @export
setGeneric("agglomerateByRank",
            signature = "x",
            function(x, ...)
                standardGeneric("agglomerateByRank"))

#' @rdname agglomerate-methods
#' @aliases agglomerateByVariable
#' @export
setGeneric("agglomerateByVariable",
            signature = "x",
            function(x, ...)
                standardGeneric("agglomerateByVariable"))

#' @rdname agglomerate-methods
#'
#' @importFrom SummarizedExperiment rowData rowData<-
#'
#' @export
setMethod("agglomerateByRank", signature = c(x = "SummarizedExperiment"),
    function(x, rank = taxonomyRanks(x)[1], onRankOnly = TRUE, na.rm = FALSE,
        empty.fields = c(NA, "", " ", "\t", "-", "_"), ...){
        # input check
        if(nrow(x) == 0L){
            stop("No data available in `x` ('x' has nrow(x) == 0L.)",
                call. = FALSE)
        }
        if(!.is_non_empty_string(rank)){
            stop("'rank' must be a non-empty single character value",
                call. = FALSE)
        }
        if(!.is_a_bool(onRankOnly)){
            stop("'onRankOnly' must be TRUE or FALSE.", call. = FALSE)
        }
        if(!.is_a_bool(na.rm)){
            stop("'na.rm' must be TRUE or FALSE.", call. = FALSE)
        }
        if(ncol(rowData(x)) == 0L){
            stop("taxonomyData needs to be populated.", call. = FALSE)
        }
        .check_taxonomic_rank(rank, x)
        .check_for_taxonomic_data_order(x)
        #

        # Make a vector from the taxonomic data.
        col <- which( taxonomyRanks(x) %in% rank )
        tax_cols <- .get_tax_cols_from_se(x)

        # if na.rm is TRUE, remove the empty, white-space, NA values from
        # tree will be pruned later, if agglomerate.tree = TRUE
        if( na.rm ){
            x <- .remove_with_empty_taxonomic_info(x, tax_cols[col],
                                                    empty.fields)
        }
        # If rank is the only rank that is available and this data is unique,
        # then the data is already 'aggregated' and no further operations
        # are needed.
        if (length(taxonomyRanks(x)) == 1L &&
            !anyDuplicated(rowData(x)[,taxonomyRanks(x)])) {
            return(x)
        }

        # get groups of taxonomy entries
        tax_factors <- .get_tax_groups(x, col = col, onRankOnly = onRankOnly)

        # merge taxa
        x <- agglomerateByVariable(x, MARGIN = "rows", f = tax_factors, ...)

        # "Empty" the values to the right of the rank, using NA_character_.
        if( col < length(taxonomyRanks(x)) ){
            badcolumns <- tax_cols[seq_along(tax_cols) > col]
            if(length(badcolumns) > 0L){
                row_data <- rowData(x)
                row_data[, badcolumns] <- NA_character_
                rowData(x) <- row_data
            }
        }
        # adjust rownames
        rownames(x) <- getTaxonomyLabels(x, empty.fields, ...,
                                        with_rank = FALSE,
                                        resolve_loops = FALSE)
        # Remove those columns from rowData that include only NAs
        x <- .remove_NA_cols_from_rowdata(x, ...)
        x <- .add_values_to_metadata(x, "agglomerated_by_rank", rank)

        # Order the data in alphabetical order
        x <- x[ order(rownames(x)), ]
    }
)

#' @rdname agglomerate-methods
#' @aliases agglomerateByVariable
#' @export
setMethod("agglomerateByVariable", signature = c(x = "SummarizedExperiment"),
            function(x, MARGIN, f, archetype = 1L, ...){
                MARGIN <- .check_MARGIN(MARGIN)
                FUN <- switch(MARGIN, .merge_rows_SE, .merge_cols_SE)
                FUN(x, f, archetype = archetype, ...)
            }
)

#' @rdname agglomerate-methods
#' @aliases agglomerateByVariable
#' @export
setMethod("agglomerateByVariable",
            signature = c(x = "TreeSummarizedExperiment"),
            function(x, MARGIN, f, archetype = 1L, mergeTree = FALSE,
                     mergeRefSeq = FALSE, ...){
                MARGIN <- .check_MARGIN(MARGIN)
                if ( MARGIN == 1L ){
                    .merge_rows_TSE(x, f, archetype = 1L, mergeTree = mergeTree,
                                   mergeRefSeq = mergeRefSeq, ...)
                }
                else{
                    .merge_cols_TSE(x, f, archetype = 1L, mergeTree = mergeTree,
                                   ...)
                }
            }
)

#' @rdname agglomerate-methods
#' @importFrom SingleCellExperiment altExp altExp<- altExps<-
#' @export
setMethod("agglomerateByRank", signature = c(x = "SingleCellExperiment"),
    function(x, ..., altexp = NULL, strip_altexp = TRUE){
        # input check
        if(!.is_a_bool(strip_altexp)){
            stop("'strip_altexp' mus be TRUE or FALSE.", call. = FALSE)
        }
        #
        if (!is.null(altexp)) {
            x <- altExp(x, altexp)
        }
        if(strip_altexp && is(x, "SingleCellExperiment")){
            altExps(x) <- NULL
        }
        callNextMethod(x, ...)
    }
)

#' @rdname agglomerate-methods
#' @export
setMethod(
    "agglomerateByRank", signature = c(x = "TreeSummarizedExperiment"),
    function(
        x, ..., agglomerate.tree = agglomerateTree, agglomerateTree = FALSE){
                # input check
                if(!.is_a_bool(agglomerate.tree)){
                    stop("'agglomerate.tree' must be TRUE or FALSE.",
                        call. = FALSE)
                }
                # If there are multipe rowTrees, it might be that multiple
                # trees are preserved after agglomeration even though the
                # dataset could be presented with one tree.
                # --> order the data so that the taxa are searched from one tree
                # first.
                if( length(rowTreeNames(x)) > 1 ){
                    x <- .order_based_on_trees(x)
                }
                # Agglomerate data
                x <- callNextMethod(x, mergeTree = agglomerate.tree, ...)
                return(x)
            }
)

################################ HELP FUNCTIONS ################################

.remove_with_empty_taxonomic_info <-
    function(x, column, empty.fields = c(NA,""," ","\t","-","_"))
        {
        tax <- as.character(rowData(x)[,column])
        f <- !(tax %in% empty.fields)
        if(any(!f)){
            x <- x[f, , drop=FALSE]
        }
        x
    }

# This function removes empty columns from rowdata. (Those that include only
# NA values)
.remove_NA_cols_from_rowdata <- function(x, remove_empty_ranks = FALSE, ...){
    # Check remove_empty_ranks
    if( !.is_a_bool(remove_empty_ranks) ){
        stop("'remove_empty_ranks' must be a boolean value.",
            call. = FALSE)
    }
    # If user wants to remove those columns
    if( remove_empty_ranks ){
        # Get rowData
        rd <- rowData(x)
        # Does teh column include data?
        columns_including_data <- apply(rd, 2, function(x){!all(is.na(x))})
        # Subset data so that it includes only columns that include data
        rd <- rd[, columns_including_data]
        # Assign it back to SE
        rowData(x) <- rd
    }
    return(x)
}

# Order the data so that taxa from tree1 comes first, then taxa
# from tree2...
.order_based_on_trees <- function(x){
    # Get rowlinks and unique trees
    links <- DataFrame(rowLinks(x))
    uniq_trees <- sort(unique(links$whichTree))
    # Get row index to the data
    links$row_i <- seq_len(nrow(x))
    # Calculate, how many rows each tree has, and add it to data
    freq <- as.data.frame(table(links$whichTree))
    links <- merge(links, freq, all.x = TRUE, all.y = FALSE,
                    by.x = "whichTree", by.y = "Var1")
    # Factorize the names of trees
    links$whichTree <- factor(links$whichTree, levels = uniq_trees)
    # Order the data back to its original order based on row indices
    links <- links[order(links$row_i), ]
    # Get the order based on size of tree and name
    order <- order(links$whichTree)
    # Order the data
    x <- x[order, ]
    return(x)
}

# Agglomerate all rowTrees found in TreeSE object. Get tips that represent
# rows and remove all others.
.agglomerate_trees <- function(x, MARGIN = 1){
    # Get right functions based on direction
    tree_names_FUN <- switch(
        MARGIN, "1" = rowTreeNames, "2" = colTreeNames, stop("."))
    links_FUN <- switch(MARGIN, "1" = rowLinks, "2" = colLinks, stop("."))
    tree_FUN <- switch(MARGIN, "1" = rowTree, "2" = colTree, stop("."))
    # Get right argument names for changeTree call
    args_names <- switch(
        MARGIN, "1" = c("x", "rowTree", "rowNodeLab", "whichRowTree"),
        "2" = c("x", "colTree", "colNodeLab", "whichColTree"),
        stop("."))
    # Get names of trees and links between trees and rows
    tree_names <- tree_names_FUN(x)
    row_links <- links_FUN(x)
    # Loop through tree names
    for( name in tree_names ){
        # Get the tree that is being agglomerated
        tree <- tree_FUN(x, name)
        # Get row links that corresponds this specific tree
        links_temp <- row_links[ row_links[["whichTree"]] == name, ]
        # If the tree represents the data, agglomerate it
        if( nrow(links_temp) > 0 ){
            # Get names of nodes that are preserved
            links_temp <- links_temp[["nodeLab"]]
            # Agglomerate the tree
            tree <- .prune_tree(tree, links_temp)
            # Change the tree with agglomerated version
            args <- list(x, tree, links_temp, name)
            names(args) <- args_names
            x <- do.call(changeTree, args)
        }
    }
    return(x)
}

# This function trims tips until all tips can be found from provided set of
# nodes
#' @importFrom ape drop.tip has.singles collapse.singles
.prune_tree <- function(tree, nodes){
    # Get those tips that can not be found from provided nodes
    remove_tips <- tree$tip.label[!tree$tip.label %in% nodes]
    # As long as there are tips to be dropped, run the loop
    while( length(remove_tips) > 0 ){
        # Drop tips that cannot be found. Drop only one layer at the time. Some
        # dataset might have taxa that are not in tip layer but they are in
        # higher rank. If we delete more than one layer at the time, we might
        # loose the node for those taxa. --> The result of pruning is a tree
        # whose all tips can be found provided nodes i.e., rows of TreeSE. Some
        # taxa might be higher rank meaning that all rows might not be in tips
        # even after pruning; these rows have still child-nodes that represent
        # other rows.
        # Suppress warning: drop all tips of the tree: returning NULL
        suppressWarnings(
            tree <- drop.tip(
                tree, remove_tips,
                trim.internal = FALSE,
                collapse.singles = FALSE)
        )
        # If all tips were dropped, the result is NULL --> stop loop
        if( is.null(tree) ){
            warning("Pruning resulted to empty tree.", call. = FALSE)
            break
        }
        # Again, get those tips of updated tree that cannot be found from
        # provided nodes
        remove_tips <- tree$tip.label[!tree$tip.label %in% nodes]
    }
    # Simplify the tree structure. Remove nodes that have only single
    # descendant.
    if( !is.null(tree) && length(tree$tip.label) > 1 && has.singles(tree) ){
        tree <- collapse.singles(tree)
    }
    return(tree)
}