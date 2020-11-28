#' Agglomerate taxa of the same type.
#'
#' \code{agglomerateByRank} can be used to sum up data based on the association
#' to certain taxonomic ranks given as \code{rowData}. Only available
#' \code{\link{taxonomicRanks}} can be used.
#'
#' @param x \code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#'
#' @param rank a single character defining a taxonomic rank. Must be a value of
#'   \code{taxonomicRanks()} function.
#'
#' @param onRankOnly \code{TRUE} or \code{FALSE}: Should information only from
#'   the specified rank used or from ranks equal and above?.
#'   (default: \code{onRankOnly = FALSE})
#'
#' @param na.rm \code{TRUE} or \code{FALSE}: Should taxa with an empty rank be
#'   removed? Use it with caution, since result with NA on the selected rank
#'   will be dropped. This setting can be tweaked by defining
#'   \code{empty.fields} to your needs. (default: \code{na.rm = TRUE})
#'
#' @param empty.fields a \code{character} value defining, which values should be
#'   regarded as empty. (Default: \code{c(NA, "", " ", "\t")}). They will be
#'   removed if \code{na.rm = TRUE} before agglomeration.
#'
#' @param agglomerateTree \code{TRUE} or \code{FALSE}: should to
#'   \code{rowTree()} also be agglomerated? (Default:
#'   \code{agglomerateTree = FALSE})
#'
#' @param ... arguments passed to \code{agglomerateByRank} function for
#'   \code{SummarizedExperiment} objects.
#'
#' @param altexp String or integer scalar specifying an alternative experiment
#'   containing the input data.
#'
#' @param strip_altexp \code{TRUE} or \code{FALSE}: Should alternative
#'   experiments be removed prior to agglomeration? This prohibits to many
#'   nested alternative experiments by default (default:
#'   \code{strip_altexp = TRUE})
#'
#' @return A taxonomically-agglomerated, optionally-pruned object of the same
#'   class \code{x}.
#'
#' @seealso
#'   \code{\link[=merge-methods]{mergeRows}}
#'
#' @name agglomerate-methods
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
#'                         agglomerateTree = TRUE)
#' nrow(x2) # same number of rows, but
#' rowTree(x1) # ... different
#' rowTree(x2) # ... tree
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


#' @rdname agglomerate-methods
#' @aliases agglomerateByRank
#'
#' @importFrom SummarizedExperiment rowData rowData<-
#'
#' @export
setMethod("agglomerateByRank", signature = c(x = "SummarizedExperiment"),
    function(x, rank = taxonomyRanks(x)[1], onRankOnly = FALSE, na.rm = FALSE,
       empty.fields = c(NA, "", " ", "\t", "-"), agglomerateTree = FALSE, ...){
        # input check
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
        # If x is aggregated on the only rank that is available,
        # then the data is already 'aggregated' and no further operations
        # are needed.
        if (length(taxonomyRanks(x)) == 1L) {
            return(x)
        }
        if(!.is_a_bool(agglomerateTree)){
            stop("'agglomerateTree' must be TRUE or FALSE.", call. = FALSE)
        }
        .check_for_taxonomic_data_order(x)
        #

        # Make a vector from the taxonomic data.
        col <- which( taxonomyRanks(x) %in% rank )
        tax_cols <- .get_tax_cols_from_se(x)

        # if na.rm is TRUE, remove the empty, white-space, NA values from
        # tree will be pruned later, if agglomerateTree = TRUE
        if( na.rm ){
            tax <- as.character(rowData(x)[,tax_cols[col]])
            f <- !(tax %in% empty.fields)
            x <- x[f, , drop=FALSE]
        }

        # get groups of taxonomy entries
        tax_factors <- .get_tax_groups(x, col = col, onRankOnly = onRankOnly)

        # merge taxa
        x <- mergeRows(x, f = tax_factors, mergeTree = agglomerateTree)

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
        rownames(x) <- .get_taxonomic_label(x, empty.fields)
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
