#' Ag           glomerate data using taxonomic information
#'
#' \code{agglomerateByRank} can be used to sum up data based on the association
#' to certain taxonomic ranks given as \code{rowData}. Only available
#' \code{\link{taxonomyRanks}} can be used.
#'
#' @param x a
#'   \code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#'   object
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
#' @param agglomerateTree \code{TRUE} or \code{FALSE}: should
#'   \code{rowTree()} also be agglomerated? (Default:
#'   \code{agglomerateTree = FALSE})
#'
#' @param ... arguments passed to \code{agglomerateByRank} function for
#'   \code{SummarizedExperiment} objects,
#'   \code{\link[=merge-methods]{mergeRows}} and
#'   \code{\link[scuttle:sumCountsAcrossFeatures]{sumCountsAcrossFeatures}}.
#'   \itemize{
#'        \item{\code{remove_empty_ranks}}{A single boolean value for selecting 
#'        whether to remove those columns of rowData that include only NAs after
#'        agglomeration. (By default: \code{remove_empty_ranks = FALSE})}
#'        \item{\code{make_unique}}{A single boolean value for selecting 
#'        whether to make rownames unique. (By default: \code{make_unique = TRUE})}
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
#' @details
#' Based on the available taxonomic data and its structure setting
#' \code{onRankOnly = TRUE} has certain implications on the interpretability of
#' your results. If no loops exist (loops meaning two higher ranks containing
#' the same lower rank), the results should be comparable. you can check for
#' loops using \code{\link[TreeSummarizedExperiment:detectLoop]{detectLoop}}.
#' 
#' Agglomeration sum up values of assays at specified taxonomic level. Certain assays,
#' e.g. those that include binary or negative values, can lead to meaningless values, 
#' when values are summed. In those cases, consider doing agglomeration first and then 
#' transformation.
#'
#' @return A taxonomically-agglomerated, optionally-pruned object of the same
#'   class as \code{x}.
#'
#' @name agglomerate-methods
#' @seealso
#' \code{\link[=merge-methods]{mergeRows}},
#' \code{\link[scuttle:sumCountsAcrossFeatures]{sumCountsAcrossFeatures}}
#'
#' @examples
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
#' # with agglomeration of the tree
#' x2 <- agglomerateByRank(GlobalPatterns, rank="Family",
#'                        agglomerateTree = TRUE)
#' nrow(x2) # same number of rows, but
#' rowTree(x1) # ... different
#' rowTree(x2) # ... tree
#' 
#'  # If assay contains binary or negative values, summing might lead to meaningless
#'  # values, and you will get a warning. In these cases, you might want to do 
#'  # agglomeration again at chosen taxonomic level.
#'  tse <- transformAssay(GlobalPatterns, method = "pa")
#'  tse <- agglomerateByRank(tse, rank = "Genus")
#'  tse <- transformAssay(tse, method = "pa")
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
#' x4 <- agglomerateByRank(GlobalPatterns, rank="Phylum", remove_empty_ranks = TRUE)
#' head(rowData(x4))
#' 
#' # If assay contains NAs, you might want to consider replacing them since summing-up
#' # NAs lead to NA
#' x5 <- GlobalPatterns
#' # Replace first value with NA
#' assay(x5)[1,1] <- NA
#' x6 <- agglomerateByRank(x5, "Kingdom")
#' head( assay(x6) )
#' # Replace NAs with 0. It is justified when we are summing-up counts
#' assay(x5)[ is.na(assay(x5)) ] <- 0
#' x6 <- agglomerateByRank(x5, "Kingdom")
#' head( assay(x6) )
#' 
#' ## Look at enterotype dataset...
#' data(enterotype)
#' ## print the available taxonomic ranks. Shows only 1 rank available
#' ## not useful for agglomerateByRank
#' taxonomyRanks(enterotype)
NULL

setGeneric("agglomerateByRank",
            signature = "x",
            function(x, ...)
                standardGeneric("agglomerateByRank"))

setGeneric("mergeFeaturesByRank",
           signature = "x",
           function(x, ...)
               standardGeneric("mergeFeaturesByRank"))

#' @rdname agglomerate-methods
#' @aliases agglomerateByRank
#'
#' @importFrom SummarizedExperiment rowData rowData<-
#'
#' @export
setMethod("agglomerateByRank", signature = c(x = "SummarizedExperiment"),
    function(x, rank = taxonomyRanks(x)[1], onRankOnly = FALSE, na.rm = FALSE,
        empty.fields = c(NA, "", " ", "\t", "-", "_"), ...){
        # input check
        if(nrow(x) == 0L){
            stop("No data available in `x` ('x' has nrow(x) == 0L.)",
                 call. = FALSE)
        }
        if(!.is_non_empty_string(rank)){
            stop("'rank' must be an non empty single character value.",
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
        # tree will be pruned later, if agglomerateTree = TRUE
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
        x <- mergeRows(x, f = tax_factors, ...)

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
                                        with_rank = FALSE, resolve_loops = FALSE)
        # Remove those columns from rowData that include only NAs
        x <- .remove_NA_cols_from_rowdata(x, ...)
        x <- .add_values_to_metadata(x, "agglomerated_by_rank", rank)
        x
    }
)

#' @rdname agglomerate-methods
#' @aliases agglomerateByRank
#'
#' @importFrom SummarizedExperiment rowData rowData<-
#'
#' @export
setMethod("mergeFeaturesByRank", signature = c(x = "SummarizedExperiment"),
          function(x, rank = taxonomyRanks(x)[1], onRankOnly = FALSE, na.rm = FALSE,
                   empty.fields = c(NA, "", " ", "\t", "-", "_"), ...){
              .Deprecated(old="agglomerateByRank", new="mergeFeaturesByRank", "Now agglomerateByRank is deprecated. Use mergeFeaturesByRank instead.")
              x <- agglomerateByRank(x)
              x
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
#' @aliases agglomerateByRank
#' @importFrom SingleCellExperiment altExp altExp<- altExps<-
#' @export
setMethod("mergeFeaturesByRank", signature = c(x = "SingleCellExperiment"),
          function(x, ..., altexp = NULL, strip_altexp = TRUE){
              .Deprecated(old="agglomerateByRank", new="mergeFeaturesByRank", "Now agglomerateByRank is deprecated. Use mergeFeaturesByRank instead.")
              x <- agglomerateByRank(x)
              x
          }
)


#' @rdname agglomerate-methods
#' @export
setMethod("agglomerateByRank", signature = c(x = "TreeSummarizedExperiment"),
          function(x, ..., agglomerateTree = FALSE){
              # input check
              if(!.is_a_bool(agglomerateTree)){
                  stop("'agglomerateTree' must be TRUE or FALSE.", call. = FALSE)
              }
              # If there are multipe rowTrees, it might be that multiple
              # trees are preserved after agglomeration even though the dataset
              # could be presented with one tree. --> order the data so that
              # the taxa are searched from one tree first.
              if( length(x@rowTree) > 1 ){
                  x <- .order_based_on_trees(x)
              }
              # Agglomerate data
              x <- callNextMethod(x, ...)
              # Agglomerate also tree, if the data includes only one
              # rowTree --> otherwise it is not possible to agglomerate
              # since all rownames are not found from individual tree.
              if(agglomerateTree){
                  if( length(x@rowTree) > 1 ){
                      warning("The dataset includes multiple tree after ",
                              "agglomeration. Agglomeration of tree is not ",
                              "possible.", call. = FALSE)
                  } else{
                      x <- addTaxonomyTree(x)
                  }
              }
              x
          }
)

#' @rdname agglomerate-methods
#' @aliases agglomerateByRank
#' @export
setMethod("mergeFeaturesByRank", signature = c(x = "TreeSummarizedExperiment"),
          function(x, ..., agglomerateTree = FALSE){
              .Deprecated(old="agglomerateByRank", new="mergeFeaturesByRank", "Now agglomerateByRank is deprecated. Use mergeFeaturesByRank instead.")
              x <- agglomerateByRank(x)
              x
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
