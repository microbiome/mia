#' @name taxonomy-methods
#'
#' @title Taxonomy related functions
#'
#' @description
#' These function work on optional data present in \code{rowData}.
#'
#' \code{taxonomyRanks} returns, which columns of \code{rowData(x)} are regarded
#' as columns containing taxonomic information.
#'
#' \code{taxonomyRankEmpty} checks, if a selected rank is empty of information.
#'
#' \code{checkTaxonomy} checks, if taxonomy information is valid and whether
#'   it contains any problems. This is a soft test, which reports some
#'   diagnostic and might mature into a data validator used upon object
#'   creation.
#'
#' \code{getTaxonomyLabels} generates a character vector per row consisting of
#'   the lowest taxonomic information possible. If data from different levels,
#'   is to be mixed, the taxonomic level is prepended by default.
#'
#'
#' \code{taxonomyTree} generates a \code{phylo} tree object from the available
#'   taxonomic information. Internally it uses
#'   \code{\link[TreeSummarizedExperiment:toTree]{toTree}} and
#'   \code{\link[TreeSummarizedExperiment:resolveLoop]{resolveLoop}} to sanitize
#'   data if needed.
#'
#' @param x a
#'   \code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#'   object
#'
#' @param from a \code{Taxa} object as returned by
#'   \code{\link[DECIPHER:IdTaxa]{IdTaxa}}
#'
#' @param rank a single character defining a taxonomic rank. Must be a value of
#'   \code{taxonomicRanks()} function.
#'
#' @param empty.fields a \code{character} value defining, which values should be
#'   regarded as empty. (Default: \code{c(NA, "", " ", "\t")}). They will be
#'   removed if \code{na.rm = TRUE} before agglomeration.
#'
#' @param with_rank \code{TRUE} or \code{FALSE}: Should the level be add as a
#'   suffix? For example: "Phylum:Crenarchaeota" (default:
#'   \code{with_rank = FALSE})
#'
#' @param make_unique \code{TRUE} or \code{FALSE}: Should the labels be made
#'   unique, if there are any duplicates? (default: \code{make_unique = TRUE})
#'
#' @param resolve_loops \code{TRUE} or \code{FALSE}: Should \code{resolveLooops}
#'   be applied to the taxonomic data? Please note that has only an effect,
#'   if the data is unique. (default: \code{resolve_loops = TRUE})
#'
#' @param ... optional arguments not used currently.
#'
#' @details
#' Taxonomic information from the \code{IdTaxa} function of \code{DECIPHER}
#' package are returned as a special class. With \code{as(taxa,"DataFrame")}
#' the information can be easily converted to a \code{DataFrame} compatible
#' with storing the taxonomic information a \code{rowData}. Please not that the
#' assigned confidence information are returned as \code{metatdata} and can
#' be accessed using \code{metadata(df)$confidence}.
#'
#' @return
#' \itemize{
#'   \item{\code{taxonomyRanks}:} {a \code{character} vector with all the
#'     taxonomic ranks found in \code{colnames(rowData(x))}}
#'   \item{\code{taxonomyRankEmpty}:} {a \code{logical} value}
#' }
#'
#' @seealso \code{\link[=agglomerate-methods]{agglomerateByRank}},
#' \code{\link[TreeSummarizedExperiment:toTree]{toTree}},
#' \code{\link[TreeSummarizedExperiment:resolveLoop]{resolveLoop}}
#'
#' @examples
#' data(GlobalPatterns)
#' GlobalPatterns
#' taxonomyRanks(GlobalPatterns)
#'
#' checkTaxonomy(GlobalPatterns)
#'
#' table(taxonomyRankEmpty(GlobalPatterns,"Kingdom"))
#' table(taxonomyRankEmpty(GlobalPatterns,"Species"))
#'
#' getTaxonomyLabels(GlobalPatterns[1:20,])
#'
#' # adding a rowTree() based on the available taxonomic information. Please
#' # note that any tree already stored in rowTree() will be overwritten.
#' x <- GlobalPatterns
#' x <- addTaxonomyTree(x)
#' x
NULL

#' @rdname taxonomy-methods
#' @format a \code{character} vector of length 8 containing the taxonomy ranks
#'   recognized. In functions this is used case insensitive.
#' @export
TAXONOMY_RANKS <- c("domain","kingdom","phylum","class","order","family",
                    "genus","species")

#' @rdname taxonomy-methods
setGeneric("taxonomyRanks", signature = c("x"),
           function(x)
             standardGeneric("taxonomyRanks"))

#' @rdname taxonomy-methods
#' @aliases taxonomicRanks
#'
#' @importFrom SummarizedExperiment rowData
#'
#' @export
setMethod("taxonomyRanks", signature = c(x = "SummarizedExperiment"),
    function(x){
        ranks <- colnames(rowData(x))
        ranks[.get_tax_cols(ranks)]
    }
)

#' @rdname taxonomy-methods
setGeneric("taxonomyRankEmpty",
           signature = "x",
           function(x, rank = taxonomyRanks(x)[1L],
                    empty.fields = c(NA, "", " ", "\t", "-", "_"))
             standardGeneric("taxonomyRankEmpty"))

#' @rdname taxonomy-methods
#' @aliases taxonomyRankEmpty
#'
#' @importFrom SummarizedExperiment rowData
#'
#' @export
setMethod("taxonomyRankEmpty", signature = c(x = "SummarizedExperiment"),
    function(x, rank = taxonomyRanks(x)[1],
           empty.fields = c(NA, "", " ", "\t", "-", "_")){
    # input check
    if(ncol(rowData(x)) == 0L){
        stop("rowData needs to be populated.", call. = FALSE)
    }
    if(!.is_non_empty_string(rank)){
        stop("'rank' must be an non empty single character value.",
             call. = FALSE)
    }
    if(!is.character(empty.fields) || length(empty.fields) == 0L){
        stop("'empty.fields' must be a character vector with one or ",
             "more value", call. = FALSE)
    }
    .check_taxonomic_rank(rank, x)
    .check_for_taxonomic_data_order(x)
    #
    rowData(x)[,rank] %in% empty.fields
    }
)

#' @rdname taxonomy-methods
setGeneric("checkTaxonomy",
           signature = "x",
           function(x, ...)
             standardGeneric("checkTaxonomy"))

#' @rdname taxonomy-methods
#' @aliases checkTaxonomy
#' @export
setMethod("checkTaxonomy", signature = c(x = "SummarizedExperiment"),
    function(x){
        tmp <- try(.check_for_taxonomic_data_order(x), silent = TRUE)
        if(is(tmp,"try-error")){
            FALSE
        }
        TRUE
    }
)

.check_taxonomic_rank <- function(rank, x){
    if( !(rank %in% taxonomyRanks(x) ) ){
        stop("'rank' must be a value from 'taxonomyRanks()'")
    }
}
.check_taxonomic_ranks <- function(ranks, x){
    if( !all(ranks %in% taxonomyRanks(x) ) ){
        stop("'ranks' must contain values from 'taxonomyRanks()'")
    }
}

#' @importFrom SummarizedExperiment rowData
.check_for_taxonomic_data_order <- function(x){
    ranks <- colnames(rowData(x))
    f <- tolower(ranks) %in% TAXONOMY_RANKS
    if(!any(f)){
        stop("no taxonomic ranks detected in rowData(). Columns with one of ",
             "the following names can be used: '",
             paste(TAXONOMY_RANKS, collapse = "', '"), "'", call. = FALSE)
    }
    m <- match(TAXONOMY_RANKS, tolower(ranks[f]))
    m <- m[!is.na(m)]
    # check that taxonomic ranks are in order. If they are all value in check
    # should be 1 or 0
    check <- unique(c(m[-1], m[length(m)]) - m )
    if(!all(check %in% c(1L,0L))){
        stop("Taxonomic ranks are not in order. Please reorder columns, which ",
             "correspond to taxonomic ranks like this:\n'",
             paste(TAXONOMY_RANKS, collapse = "', '"), "'.",
             call. = FALSE)
    }
}


#' @rdname taxonomy-methods
setGeneric("getTaxonomyLabels",
           signature = "x",
           function(x, ...)
               standardGeneric("getTaxonomyLabels"))

#' @rdname taxonomy-methods
#' @aliases checkTaxonomy
#' @export
setMethod("getTaxonomyLabels", signature = c(x = "SummarizedExperiment"),
    function(x, empty.fields = c(NA, "", " ", "\t", "-", "_"),
             with_rank = FALSE, make_unique = TRUE, resolve_loops = FALSE){
        # input check
        if(ncol(rowData(x)) == 0L){
            stop("rowData needs to be populated.", call. = FALSE)
        }
        .check_for_taxonomic_data_order(x)
        if(!is.character(empty.fields) || length(empty.fields) == 0L){
            stop("'empty.fields' must be a character vector with one or ",
                 "more values.", call. = FALSE)
        }
        if(!.is_a_bool(with_rank)){
            stop("'with_rank' must be TRUE or FALSE.", call. = FALSE)
        }
        if(!.is_a_bool(make_unique)){
            stop("'make_unique' must be TRUE or FALSE.", call. = FALSE)
        }
        if(!.is_a_bool(resolve_loops)){
            stop("'resolve_loops' must be TRUE or FALSE.", call. = FALSE)
        }
        #
        dup <- duplicated(rowData(x)[,taxonomyRanks(x)])
        if(any(dup)){
            td <- apply(rowData(x)[,taxonomyRanks(x)],1L,paste,collapse = "___")
            td_non_dup <- td[!dup]
            m <- match(td, td_non_dup)
        }
        ans <- .get_taxonomic_label(x[!dup,],
                                    empty.fields = empty.fields,
                                    with_rank = with_rank,
                                    make_unique = make_unique,
                                    resolve_loops = resolve_loops)
        if(any(dup)){
            ans <- ans[m]
        }
        # last resort - this happens, if annotation data contains ambiguous data
        # sometimes labeled as "circles"
        if(make_unique && anyDuplicated(ans)){
            dup <- which(ans %in% ans[which(duplicated(ans))])
            ans[dup] <- make.unique(ans[dup], sep = "_")
        }
        ans
    }
)

#' @importFrom IRanges CharacterList LogicalList
.get_tax_ranks_selected <- function(x,rd, tax_cols, empty.fields){
    # We need DataFrame here to handle cases with a single entry in tax_cols
    charlist <- CharacterList(t(rd[,tax_cols, drop=FALSE]))
    tax_ranks_non_empty <- !is.na(charlist) &
        !LogicalList(lapply(charlist,"%in%",empty.fields))

    tax_ranks_non_empty <- t(as(tax_ranks_non_empty,"matrix"))
    tax_ranks_selected <- apply(tax_ranks_non_empty,1L,which)
    if(any(lengths(tax_ranks_selected) == 0L)){
        if(!anyDuplicated(rownames(x))){
            return(NULL)
        }
        stop("Only empty taxonomic information detected. Some rows contain ",
             "only entries selected by 'empty.fields'. Cannot generated ",
             "labels. Try option na.rm = TRUE in the function call.",
             call. = FALSE)
    }
    #
    if(is.matrix(tax_ranks_selected)){
        tax_ranks_selected <- apply(tax_ranks_selected,2L,max)
    } else if(is.list(tax_ranks_selected)) {
        tax_ranks_selected <- lapply(tax_ranks_selected,max)
        tax_ranks_selected <- unlist(tax_ranks_selected)
    } else if(is.vector(tax_ranks_selected)){
        tax_ranks_selected <- max(tax_ranks_selected)
    } else {
        stop(".")
    }
    tax_ranks_selected
}

.add_taxonomic_type <- function(rd, ans, tax_cols_selected){
    sep <- rep(":", length(ans))
    tax_cols_selected <- unlist(tax_cols_selected)
    # sep[tax_cols_selected != max(tax_cols_selected)] <- "::"
    types <- colnames(rd)[tax_cols_selected]
    ans <- paste0(types, sep, ans)
    ans
}

.get_taxonomic_label <- function(x,
                                 empty.fields = c(NA, "", " ", "\t", "-", "_"),
                                 with_rank = FALSE, make_unique = TRUE,
                                 resolve_loops = FALSE){
    rd <- rowData(x)
    tax_cols <- .get_tax_cols_from_se(x)
    tax_ranks_selected <- .get_tax_ranks_selected(x, rd, tax_cols, empty.fields)
    if(is.null(tax_ranks_selected)){
        return(rownames(x))
    }
    tax_cols_selected <- tax_cols[tax_ranks_selected]
    # resolve loops
    if(resolve_loops){
        td <- as.data.frame(rd[,tax_cols])
        td <- suppressWarnings(resolveLoop(td))
        rd[,tax_cols] <- as(td,"DataFrame")
        rm(td)
    }
    #
    all_same_rank <- length(unique(tax_cols_selected)) == 1L
    ans <- mapply("[",
                  as.data.frame(t(as.data.frame(rd))),
                  tax_cols_selected,
                  SIMPLIFY = FALSE)
    ans <- unlist(ans, use.names = FALSE)
    if(with_rank || !all_same_rank){
        ans <- .add_taxonomic_type(rd, ans, tax_cols_selected)
    }
    ans
}

#' @rdname taxonomy-methods
setGeneric("taxonomyTree",
           signature = "x",
           function(x, ...)
               standardGeneric("taxonomyTree"))

#' @rdname taxonomy-methods
#' @export
setMethod("taxonomyTree", signature = c(x = "SummarizedExperiment"),
    function(x){
        td <- rowData(x)[,taxonomyRanks(x)]
        # Remove empty taxonomic levels
        td <- td[,!vapply(td,function(tl){all(is.na(tl))},logical(1))]
        # Make information unique
        td_NA <- DataFrame(lapply(td,is.na))
        td <- as.data.frame(td)
        td <- as(suppressWarnings(resolveLoop(td)),"DataFrame")
        # Build tree
        tree <- toTree(td)
        tree$tip.label <- paste0(colnames(td)[ncol(td)],":",tree$tip.label)
        # remove empty nodes
        for(i in rev(seq_len(ncol(td)))){
            if(any(td_NA[,i])){
                to_drop <- paste0(colnames(td)[i],":",td[,i][td_NA[,i]])
                tree <- ape::drop.tip(tree,
                                      to_drop,
                                      trim.internal = FALSE,
                                      collapse.singles = FALSE)
            }
        }
        tree
    }
)

#' @rdname taxonomy-methods
setGeneric("addTaxonomyTree",
           signature = "x",
           function(x, ...)
               standardGeneric("addTaxonomyTree"))

#' @rdname taxonomy-methods
#' @export
setMethod("addTaxonomyTree", signature = c(x = "SummarizedExperiment"),
    function(x){
        #
        tree <- taxonomyTree(x)
        x <- as(x,"TreeSummarizedExperiment")
        rownames(x) <- getTaxonomyLabels(x, with_rank = TRUE,
                                         resolve_loops = TRUE,
                                         make_unique = FALSE)
        x <- changeTree(x, tree, rownames(x))
        x
    }
)


################################################################################
# helper functions

.get_tax_cols_logical <- function(x){
    tolower(x) %in% TAXONOMY_RANKS
}

.get_tax_cols <- function(x){
    which(.get_tax_cols_logical(x))
}

#' @importFrom SummarizedExperiment rowData
.get_tax_cols_from_se <- function(x){
    .get_tax_cols(colnames(rowData(x)))
}

#' @importFrom SummarizedExperiment rowData
.get_tax_groups <- function(x, col, onRankOnly = FALSE){
    tax_cols <- .get_tax_cols_from_se(x)
    tax_col_n <- seq_along(tax_cols)
    if(length(tax_col_n) < col){
        stop(".")
    }
    if(onRankOnly){
        groups <- rowData(x)[,tax_cols[tax_col_n == col],drop=TRUE]
    } else {
        groups <- rowData(x)[,tax_cols[tax_col_n <= col],drop=FALSE]
        groups <- apply(groups,1L,paste,collapse="_")
    }
    factor(groups, unique(groups))
}

################################################################################
# IDTAXA to DataFrame conversion

.idtaxa_to_DataFrame <- function(from){
    ranks <- CharacterList(lapply(from,"[[","rank"))
    conf <- NumericList(lapply(from,"[[","confidence"))
    taxa <- CharacterList(lapply(from,"[[","taxon"))
    # even out the lengths
    l <- lengths(ranks)
    ml <- max(l)
    diff <- ml - l
    add <- CharacterList(lapply(diff,rep,x=NA))
    ranks <- pc(ranks,add)
    conf <- pc(conf,as(add,"NumericList"))
    taxa <- pc(taxa,add)
    # convert to DataFrame
    names <- unique(unlist(ranks))
    names <- names[!is.na(names)]
    taxa <- DataFrame(as.matrix(taxa))
    colnames(taxa) <- names
    conf <- DataFrame(as.matrix(conf))
    colnames(conf) <- names
    # subset to valid taxonomic information
    f <- tolower(names) %in% TAXONOMY_RANKS
    taxa <- taxa[,f]
    conf <- conf[,f]
    # combine with confidence data
    metadata(taxa)$confidence <- conf
    #
    taxa
}

#' @rdname taxonomy-methods
#' @export
IdTaxaToDataFrame <- .idtaxa_to_DataFrame
