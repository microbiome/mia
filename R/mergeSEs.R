#' Merge SE objects into single SE object.
#' 
#' @param x a \code{\link{SummarizedExperiment}} object or a list of 
#' \code{\link{SummarizedExperiment}} objects.
#' 
#' @param y a \code{\link{SummarizedExperiment}} object when \code{x} is a
#' \code{\link{SummarizedExperiment}} object. Disabled when \code{x} is a list.
#' 
#' @param assay_name A single character value for selecting the
#' \code{\link[SummarizedExperiment:SummarizedExperiment-class]{assay}}
#' to be merged. (By default: \code{assay_name = "counts"})
#' 
#' @param join A single character value for selecting the joining method.
#' Must be 'full', 'inner', 'left', or 'right'. 'left' and 'right' are disabled
#' when more than two objects are being merged.  (By default: \code{join = "full"})
#' 
#' @param missing_values NA, 0, or a single character values specifying the notation
#' of missing values. (By default: \code{missing_values = NA})
#' 
#' @param collapse_samples A boolean value for selecting whether to collapse identically
#' named samples to one. (By default: \code{collapse_samples = FALSE})
#' 
#' @param verbose A single boolean value to choose whether to show messages. 
#' (By default: \code{verbose = TRUE})
#'
#' @param ... optional arguments (not used).
#'
#' @return A single \code{SummarizedExperiment} object.
#'
#' @details
#' This function merges multiple \code{SummarizedExperiment} objects. It combines
#' \code{rowData}, \code{assays}, and \code{colData} so that the output includes
#' each unique row and column ones. The merging is done based on \code{rownames} and
#' \code{colnames}. \code{rowTree} and \code{colTree} are preserved if linkage
#' between rows/cols and the tree is found.
#' 
#' Equally named rows are interpreted as equal. Further
#' matching based on \code{rowData} is not done. For samples, collapsing 
#' is disabled by default meaning that equally named samples that are stored 
#' in different objects are interpreted as unique. Collapsing can be enabled 
#' with \code{collapse_samples = TRUE} when equally named samples describe the same
#' sample. 
#' 
#' If, for example, all rows are not shared with
#' individual objects, there are missing values in \code{assays}. The notation of missing
#' can be specified with the \code{missing_values} argument. If input consists of
#' \code{TreeSummarizedExperiment} objects, also \code{rowTree}, \code{colTree}, and
#' \code{referenceSeq} are preserved if possible.
#' 
#' Compared to \code{cbind} and \code{rbind} \code{mergeSEs} 
#' allows more freely merging since \code{cbind} and \code{rbind} expect 
#' that rows and columns are matching, respectively.
#' 
#' You can choose joining methods from \code{'full'}, \code{'inner'},
#'  \code{'left'}, and  \code{'right'}. In all the methods, all the samples are 
#'  included in the result object. However, with different methods, it is possible 
#'  to choose which rows are included.
#' 
#' \itemize{
#'   \item{\code{full} -- all unique features}
#'   \item{\code{inner} -- all shared features}
#'   \item{\code{left} -- all the features of the first object}
#'   \item{\code{right} -- all the features of the second object}
#' }
#' 
#' You can also doe e.g., a full join by using a function \code{full_join} which is 
#' an alias for \code{mergeSEs}. Also other joining methods have dplyr-like aliases.
#' 
#' The output depends on the input. If the input contains \code{SummarizedExperiment}
#' object, then the output will be \code{SummarizedExperiment}. When all the input
#' objects belong to \code{TreeSummarizedExperiment}, the output will be 
#' \code{TreeSummarizedExperiment}.
#'
#' @seealso
#' \itemize{
#'   \item{\code{TreeSummarizedExperiment::cbind}}
#'   \item{\code{TreeSummarizedExperiment::rbind}}
#'   \item{\code{\link[dplyr:full_join]{full_join}}}
#'   \item{\code{\link[dplyr:inner_join]{inner_join}}}
#'   \item{\code{\link[dplyr:left_join]{left_join}}}
#'   \item{\code{\link[dplyr:right_join]{right_join}}}
#' }
#'
#' @name mergeSEs
#' @export
#'
#' @author Leo Lahti and Tuomas Borman. Contact: \url{microbiome.github.io}
#' 
#' @examples
#' data(GlobalPatterns)
#' data(esophagus)
#' data(enterotype)
#' 
#' # Take only subsets so that it wont take so long
#' tse1 <- GlobalPatterns[1:100, ]
#' tse2 <- esophagus
#' tse3 <- enterotype[1:100, ]
#' 
#' # Merge two TreeSEs
#' tse <- mergeSEs(tse1, tse2)
#' 
#' # Merge a list of TreeSEs
#' list <- SimpleList(tse1, tse2, tse3)
#' tse <- mergeSEs(list, assay_name = "counts", missing_values = 0)
#' tse
#' 
#' # With 'join', it is possible to specify the merging method. Subsets are used
#' # here just to show the functionality
#' tse_temp <- mergeSEs(tse[1:10, 1:10], tse[5:100, 11:20], join = "left")
#' tse_temp
#' 
#' # You can also do a left_join by using alias "left_join"
#' tse_temp <- left_join(tse[1:10, 1:10], tse[5:100, 11:20])
#' 
#' # If your objects contain samples that describe one and same sample,
#' # you can collapse equally named samples to one by specifying 'collapse_samples'
#' tse_temp <- inner_join(list(tse[1:10, 1], tse[1:20, 1], tse[1:5, 1]), 
#'                        collapse_samples = TRUE)
#' tse_temp
#' 
NULL

################################### Generic ####################################

#' @rdname mergeSEs
#' @export
setGeneric("mergeSEs", signature = c("x"),
        function(x, ... )
            standardGeneric("mergeSEs"))

###################### Function for SimpleList of TreeSEs ######################

#' @rdname mergeSEs
#' @export
setMethod("mergeSEs", signature = c(x = "SimpleList"),
        function(x, assay_name = "counts", join = "full", 
                 missing_values = NA, collapse_samples = FALSE, verbose = TRUE, 
                 ... ){
            ################## Input check ##################
            # Check the objects 
            class <- .check_objects_and_give_class(x)
            # Can the assay_name the found form all the objects
            assay_name_bool <- .assays_cannot_be_found(assay_name = assay_name, x)
            if( any(assay_name_bool) ){
                stop("'assay_name' must specify an assay from assays. 'assay_name' ",
                     "cannot be found at least in one SE object.",
                     call. = FALSE)
            }
            # Check join
            if( !(.is_a_string(join) &&
                join %in% c("full", "inner", "left", "right") ) ){
                stop("'join' must be 'full', 'inner', 'left', or 'right'.",
                     call. = FALSE)
            }
            # Check if join is not available
            if( length(x) > 2 &&
                join %in% c("left", "right") ){
                stop("Joining method 'left' and 'right' are not available ",
                     "when more than two objects are being merged.",
                     call. = FALSE)
            }
            # Is missing_values one of the allowed ones
            missing_values_bool <- length(missing_values) == 1L &&
                (is.numeric(missing_values) && missing_values == 0) ||
                .is_a_string(missing_values) || is.na(missing_values)
            # If not then give error
            if(  !missing_values_bool ){
                stop("'missing_values' must be 0, NA, or a single character value.",
                     call. = FALSE)
            }
            # Check collapse_samples
            if( !.is_a_bool(collapse_samples) ){
                stop("'collapse_samples' must be TRUE or FALSE.",
                     call. = FALSE)
            }
            # Check verbose
            if( !.is_a_bool(verbose) ){
                stop("'verbose' must be TRUE or FALSE.",
                     call. = FALSE)
            }
            ################ Input check end ################
            # Give message if TRUE
            if( verbose ){
                message("Merging with ", join, " join...")
                message("1/", length(x), appendLF = FALSE)
            }
            # Merge objects
            tse <- .merge_SEs(x, class, join, assay_name, 
                             missing_values, collapse_samples, verbose)
            return(tse)
        }
)

########################### Function for two TreeSEs ###########################

#' @rdname mergeSEs
#' @export
setMethod("mergeSEs", signature = c(x = "SummarizedExperiment"),
        function(x, y = NULL, ...){
            ################## Input check ##################
            # Check y
            if( !(is(y, "SummarizedExperiment")) ){
                stop("'y' must be a 'SummarizedExperiment' object.",
                     call. = FALSE)
            } 
            ################ Input check end ################
            # Create a list based on TreeSEs
            list <- SimpleList(x, y)
            # Call the function for list
            mergeSEs(list, ...)
        }
)

########################### Function for list TreeSEs ##########################

#' @rdname mergeSEs
#' @export
setMethod("mergeSEs", signature = c(x = "list"),
          function(x, ...){
              # Convert into a list
              x <- SimpleList(x)
              # Call the function for list
              mergeSEs(x, ...)
          }
)

################################# full_join ####################################

#' @rdname mergeSEs
#' @export
setGeneric("full_join", signature = c("x"),
    function(x, ...)
        standardGeneric("full_join"))

#' @rdname mergeSEs
#' @export
setMethod("full_join", signature = c(x = "ANY"),
    function(x, ...){
        mergeSEs(x, join = "full", ...)
    }
)

################################# inner_join ###################################

#' @rdname mergeSEs
#' @export
setGeneric("inner_join", signature = c("x"),
    function(x, ...)
        standardGeneric("inner_join"))

#' @rdname mergeSEs
#' @export
setMethod("inner_join", signature = c(x = "ANY"),
    function(x, ...){
        mergeSEs(x, join = "inner", ...)
    }
)

################################# left_join ####################################

#' @rdname mergeSEs
#' @export
setGeneric("left_join", signature = c("x"),
    function(x, ...)
        standardGeneric("left_join"))

#' @rdname mergeSEs
#' @export
setMethod("left_join", signature = c(x = "ANY"),
    function(x, ...){
        mergeSEs(x, join = "left", ...)
    }
)

################################# right_join ###################################

#' @rdname mergeSEs
#' @export
setGeneric("right_join", signature = c("x"),
    function(x, ...)
        standardGeneric("right_join"))

#' @rdname mergeSEs
#' @export
setMethod("right_join", signature = c(x = "ANY"),
    function(x, ...){
        mergeSEs(x, join = "right", ...)
    }
)

################################ HELP FUNCTIONS ################################

################################## .merge_SEs ##################################
# This function merges SE objects into one SE

# Input: A list of SEs
# Output: SE

#' @importFrom SingleCellExperiment SingleCellExperiment
.merge_SEs <- function(x, class, join, assay_name, 
                      missing_values, collapse_samples, verbose){
    # Take first element and remove it from the list
    tse <- x[[1]]
    x[[1]] <- NULL
    
    # If class is TreeSE, get trees and links
    tse_args <- list(
        rowTrees = NULL,
        colTrees = NULL,
        refSeqs = NULL
    )
    if( class == "TreeSummarizedExperiment" ){
        tse_args <- .get_TreeSE_args(tse, tse_args)
    }
    
    # Get the function based on class
    FUN <- switch(class,
                  TreeSummarizedExperiment = .get_TreeSummarizedExperiment_data,
                  SingleCellExperiment = .get_SingleCellExperiment_data,
                  SummarizedExperiment = .get_SummarizedExperiment_data,
                  )
    
    # Get the data in a list
    args <- do.call(FUN, args = list(tse = tse, assay_name = assay_name))
    
    # Get the function based on class
    FUN_constructor <- switch(class,
                              TreeSummarizedExperiment = TreeSummarizedExperiment,
                              SingleCellExperiment = SingleCellExperiment,
                              SummarizedExperiment = SummarizedExperiment
    )
    tse <- do.call(FUN_constructor, args = args)
    
    # Get the function based on class
    FUN <- switch(class,
                  TreeSummarizedExperiment = .merge_TreeSummarizedExperiments,
                  SingleCellExperiment = .merge_SingleCellExperiments,
                  SummarizedExperiment = .merge_SummarizedExperiments,
    )
    
    # Lopp through individual TreeSEs and add them to tse
    if( length(x) > 0 ){
        for( i in 1:length(x) ){
            # Give message if TRUE
            if( verbose ){
                message("\r", i+1, "/", length(x)+1, appendLF = FALSE)
            }
            
            # Get the ith object
            temp <- x[[i]]

            # Modify names if specified
            if( !collapse_samples ){
                temp <- .get_unique_sample_names(tse, temp, i+1)
            }
            # Merge data
            args <- do.call(FUN, args = list(
                tse_original = tse,
                tse = temp,
                join = join,
                assay_name = assay_name,
                missing_values = missing_values
                ))
            # Create an object
            tse <- do.call(FUN_constructor, args = args)
            
            # If class is TreeSE, get trees and links
            # If class is TreeSE, get trees and links
            if( class == "TreeSummarizedExperiment" ){
                tse_args <- .get_TreeSE_args(temp, tse_args)
            }
        }
        rowTrees <- tse_args$rowTrees
        colTrees <- tse_args$colTrees
        refSeqs <- tse_args$refSeqs
        # If class is TreeSE, get trees and links
        if( !is.null(rowTrees) ){
            tse <- .check_and_add_trees(tse, rowTrees, "row")
        }
        # If class is TreeSE, get trees and links
        if( !is.null(colTrees) ){
            tse <- .check_and_add_trees(tse, colTrees, "col")
        }
        # If class is TreeSE, get trees and links
        if( !is.null(refSeqs) ){
            tse <- .check_and_add_trees(tse, colTrees, "col")
        }
        # Add new line to, so that possible warning message has new line
        if( verbose ){
            message("")
        }
    }
    return(tse)
}

.check_and_add_trees <- function(tse, trees_and_links, MARGIN){
    trees <- trees_and_links$trees
    links <- trees_and_links$links
    names <- ifelse(MARGIN == "row", rownames(tse), colnames(tse) )
    
    if( !all(names %in% links[["names"]]) ){
        warning(MARGIN, "Tree(s) does not match with the data so it is discarded.",
                call. = FALSE)
        return(tse)
    }
    
    tree_labs <- lapply(trees, FUN = function(x){
        labs <- c( x$tip.label, x$node.label )
        return(labs)
    })
    # Make names unique
    names(tree_labs) <- names(trees) <- make.unique(names(trees))
    
    # Get links
    # Remove duplicates
    links <- links[ !duplicated(links[["names"]]), ]
    node_labs <- links[["nodeLab"]]
    # Get corresponding row/col
    names(node_labs) <- links[["names"]]
    
    # Check whether nodes can be found from trees
    result <- lapply(tree_labs, FUN = function(x){
        c( node_labs %in% x )
    })
    # Create a data.frame
    result <- as.data.frame(result)
    rownames(result) <- node_labs
    
    # If there are rows/cols that cannot be linked with trees
    if( any( rowSums(result) == 0 ) ){
        warning(MARGIN, "Tree(s) does not match with the data so it is discarded.",
                call. = FALSE)
        return(tse)
    }
    
    tree_names <- NULL
    for( i in seq_len(ncol(result)) ){
        combinations <- combn(result, i, simplify = FALSE)
        res <- lapply(combinations, FUN = function(x){
            all( rowSums(x) > 0 )
        })
        res <- unlist(res)
        if( any(res) ){
            # Take the first combination that have all the rows
            combinations <- combinations[[which(res)[[1]]]]
            # Take the tree name
            tree_names <- colnames(combinations)
            break
        }
    }
    trees <- trees[tree_names]
    result <- result[ , tree_names, drop = FALSE ]
    whichTree <- apply(result, 1, FUN = function(x){
        names(result)[x == TRUE]
    }
    )
    
    rowLinks <- data.frame(
        nodeLab = node_labs,
        whichTree = whichTree
    )
    
    temp <- apply(rowLinks, 1, FUN = function(x){
        temp <- c()
        nodeNum <- convertNode(trees[[x[["whichTree"]]]], x[["nodeLab"]])
        nodeLab_alias <- paste0("alias_", nodeNum)
        isLeaf <- isLeaf(trees[[x[["whichTree"]]]], x[["nodeLab"]])
        
        return(c(nodeNum = nodeNum, nodeLab_alias = nodeLab_alias, isLeaf = isLeaf))
    })
    temp <- t(temp)
    colnames(temp) = c("nodeNum", "nodeLab_alias", "isLeaf")
    rowLinks <- cbind(rowLinks, temp)
    rowLinks <- rowLinks[ , c("nodeLab", "nodeNum", "nodeLab_alias", "isLeaf", "whichTree")]
    rowLinks <- LinkDataFrame(
        nodeLab = rowLinks[["nodeLab"]],
        nodeNum = as.integer(rowLinks[["nodeNum"]]),
        nodeLab_alias = rowLinks[["nodeLab_alias"]],
        isLeaf = as.logical(rowLinks[["isLeaf"]]),
        whichTree = rowLinks[["whichTree"]]
    )
    
    if(MARGIN == "row" ){
        tse@rowTree <- trees
        tse@rowLinks <- rowLinks
    } else{
        tse@colTree <- trees
        tse@colLinks <- rowLinks
    }
    
    tse
}

.get_TreeSE_args <- function(tse, tse_args){
    
    if( !is.null(tse@rowTree) ){
        rowLinks <- DataFrame(rowLinks(tse))
        rowLinks$names <- rownames(tse)
        rowTrees <- list(
            trees = tse@rowTree,
            links = rowLinks
        )
        if( is.null(tse_args$rowTrees) ){
            tse_args$rowTrees <- rowTrees
        }
        else{
            tse_args$rowTrees <- list( 
                trees = c(tse_args$rowTrees$trees, rowTrees$trees),
                links = rbind(tse_args$rowTrees$links, rowTrees$links)
                ) 
        }
    }
    if( !is.null(tse@colTree) ){
        colLinks <- DataFrame(colLinks(tse))
        colLinks$names <- colnames(tse)
        
        colTrees <- list(
            trees = tse@colTree,
            links = colLinks
        )
        if( is.null(tse_args$colTrees) ){
            tse_args$colTrees <- colTrees
        }
        else{
            tse_args$colTrees <- list( 
                trees = c(tse_args$colTrees$trees, colTrees$trees),
                links = rbind(tse_args$colTrees$links, colTrees$links)
            ) 
        }
    }
    if( !is.null(referenceSeq(tse)) ){
        refSeqs <- list(
            referenceSeq(tse)
        )
        if( !is.null(tse_args$refSeqs) ){
            tse_args$refSeqs <- refSeqs
        }
        else{
            tse_args$refSeqs <- c( tse_args$refSeqs, refSeqs ) 
        }
    }
    return(tse_args)
}
###################### .get_TreeSummarizedExperiment_data ######################
# This function gets the desired data from one TreeSE and creates a list of 
# arguments containing the data

# Input; TreeSE
# Output: A list of arguments
.get_TreeSummarizedExperiment_data <- function(tse, assay_name){
    # # Get rowTree and colTree
    # row_tree <- rowTree(tse)
    # row_links <- rowLinks(tse)
    # col_tree <- colTree(tse)
    # col_links <- colLinks(tse)
    # Get a list of arguments of SCE object
    args <- .get_SingleCellExperiment_data(tse, assay_name)
    # # Add TreeSE-specific slots
    # args$rowTree <- row_tree
    # args$rowNodeLab <- row_links[ , "nodeLab"]
    # args$colTree <- col_tree
    # args$colNodeLab <- col_links[ , "nodeLab"]
    return(args)
}

######################## .get_SingleCellExperiment_data ########################
# This function gets the desired data from one SCE object and creates a list of 
# arguments containing the data

# Input; SCE
# Output: A list of arguments
.get_SingleCellExperiment_data <- function(tse, assay_name){
    # reducedDim is additional slot for SCE compared to SE. 
    # However, merging reducedDims leads to non-meaningful data
    # Get the arguments of SE object
    args <- .get_SummarizedExperiment_data(tse, assay_name)
    return(args)
}

######################## .get_SummarizedExperiment_data ########################
# This function gets the desired data from one SE object and creates a list of 
# arguments containing the data

# Input: SE
# Output: A list of arguments
.get_SummarizedExperiment_data <- function(tse, assay_name){
    # Remove all information but rowData, colData, metadata and assay
    row_data <- rowData(tse)
    col_data <- colData(tse)
    assay <- assay(tse, assay_name)
    assays <- SimpleList(name = assay)
    names(assays) <- assay_name
    metadata <- metadata(tse)
    # Create a list of arguments
    args <- list(assays = assays,
                rowData = row_data,
                colData = col_data,
                metadata = metadata
    )
    return(args)
    
}

######################## .check_objects_and_give_class #########################
# This function checks that the object are in correct format

# Input: a list of objects
# Output: A shared class of objects
.check_objects_and_give_class <- function(x){
    # Allowed classes
    allowed_classes <- c("TreeSummarizedExperiment", "SingleCellExperiment", "SummarizedExperiment")
    
    # Get the class based on hierarchy TreeSE --> SCE --> SE
    if( all( unlist( lapply(x, is, class2 = allowed_classes[[1]]) ) ) ){
        class <- allowed_classes[1]
    } else if( all( unlist( lapply(x, is, class2 = allowed_classes[[2]]) ) ) ){
        class <- allowed_classes[2]
    } else if( all( unlist( lapply(x, is, class2 = allowed_classes[[3]]) ) ) ){
        class <- allowed_classes[3]
    # If there is an object that does not belong to these classes give an error
    } else{
        stop("Input includes an object that is not 'SummarizedExperiment'.",
             call. = FALSE)
    }
    # If there are multiple classes, give a warning
    if( length(unique( unlist(lapply(x, function(y){ class(y)})) )) > 1 ){
        warning("The Input consist of multiple classes. ",
                "The output is '", class, "'.",
                call. = FALSE)
    }
    # Check that there are no object with no dimensions
    if( any(unlist(lapply(x, FUN = function(y){ nrow(y) == 0 || ncol(y) == 0}))) ){
        stop("Input includes an object that has either no columns or/and no rows.",
             call. = FALSE)
    }
    return(class)
}

########################### .assays_cannot_be_found ############################
# This function checks that the assay can be found from TreeSE objects of a list.

# Input: the name of the assay and a list of TreeSE objects
# Output: A list of boolean values
.assays_cannot_be_found <- function(assay_name, x){
    # Check if the assay_name can be found. If yes, then FALSE. If not, then TRUE
    list <- lapply(x, .assay_cannot_be_found, assay_name = assay_name)
    # Unlist the list
    result <- unlist(list)
    return(result)
}

############################ .assay_cannot_be_found #############################
# This function checks that the assay can be found from TreeSE. If it cannot be found
# --> TRUE, if it can be found --> FALSE

# Input: the name of the assay and TreSE object
# Output: TRUE or FALSE
.assay_cannot_be_found <- function(assay_name, tse){
    # Check if the assay_name can be found. If yes, then FALSE. If not, then TRUE
    tryCatch(
        {
            .check_assay_present(assay_name, tse)
            return(FALSE)
            
        },
        error = function(cond) {
            return(TRUE)
        }
    )
}

########################### .get_unique_sample_names ###########################
# This function convert colnames unique

# Input: TreeSEs
# Output: One TreeSE with unique sample names compared to other TreeSE
.get_unique_sample_names <- function(tse1, tse2, iteration){
    # Get indices of those sample names that match
    ind <-  colnames(tse2) %in% colnames(tse1)
    # Get duplicated sample names
    duplicated_colnames <-  colnames(tse2)[ind]
    if( length(duplicated_colnames) > 0 ) {
        # Add the number of object to duplicated sample names
        duplicated_colnames <- paste0(duplicated_colnames, "_", iteration)
        # Add new sample names to the tse object
        colnames(tse2)[ind] <- duplicated_colnames
    }
    return(tse2)
}

###################### .merge_TreeSummarizedExperiments ########################
# This function merges the data of two TreeSE objects into one set of arguments that
# can be feed to create a single object.

# Input: Two TreeSEs, the name of the assay, joining method, and the value to
# denote missing values that might occur when object do not share same features, e.g.
# Output: A list of arguments
.merge_TreeSummarizedExperiments <- function(tse_original, tse, join,  
                                             assay_name, missing_values){
    # Merge data to get a list of arguments
    args <- .merge_SingleCellExperiments(tse_original, tse, join,
                                        assay_name, missing_values)
    
    # # Check if 1st rowTree matches with data
    # tree_args <- .check_if_tree_matches_with_data(
    #     tse_original@rowTree,
    #     tse@rowTree,
    #     rowLinks(tse_original),
    #     rowLinks(tse),
    #     rownames(tse_original),
    #     rownames(tse),
    #     rownames(args$rowData)
    #     )
    # # If it is not NULL, then add arguments
    # if( !is.null(tree_args) ){
    #     args$rowTree <- tree_args$tree
    #     args$rowNodeLab <- tree_args$node_labs
    # # Otherwise test the 2nd tree
    # } else{
    #     # Get tree args
    #     tree_args <- .check_if_tree_matches_with_data(
    #         tse_original@rowTree,
    #         tse@rowTree,
    #         rowLinks(tse_original),
    #         rowLinks(tse),
    #         rownames(tse_original),
    #         rownames(tse),
    #         rownames(args$rowData)
    #     )
    #     # If not NULL, then add
    #     if( !is.null(tree_args) ){
    #         args$rowTree <- tree_args$tree
    #         args$rowNodeLab <- tree_args$node_labs
    #     }
    # }
    # 
    # # Check if 1st colTree matches with data
    # tree_args <- .check_if_tree_matches_with_data(
    #     colTree(tse_original),
    #     colLinks(tse_original),
    #     colLinks(tse),
    #     colnames(tse_original),
    #     colnames(tse),
    #     rownames(args$colData)
    # )
    # # If it is not NULL, then add arguments
    # if( !is.null(tree_args) ){
    #     args$colTree <- tree_args$tree
    #     args$colNodeLab <- tree_args$node_labs
    #     # Otherwise test the 2nd tree
    # } else{
    #     # Get tree args
    #     tree_args <- .check_if_tree_matches_with_data(
    #         colTree(tse),
    #         colLinks(tse_original),
    #         colLinks(tse),
    #         colnames(tse_original),
    #         colnames(tse),
    #         rownames(args$colData)
    #     )
    #     # If not NULL, then add
    #     if( !is.null(tree_args) ){
    #         args$colTree <- tree_args$tree
    #         args$colNodeLab <- tree_args$node_labs
    #     }
    # }
    # 
    # # Reference sequences
    # ref_seqs1 <- referenceSeq(tse_original)
    # ref_seqs2 <- referenceSeq(tse)
    # ####################################### BIND refeqs "#######################
    # # If names of 1st sequences match with data, add refseq1
    # if( length(rownames(args$rowData)) > 0 && !is.null( ref_seqs1 ) &&
    #     all( rownames(args$rowData) %in% names(ref_seqs1) ) ){
    #     # Put sequences into correct order
    #     ref_seqs1 <- ref_seqs1[ match(rownames(args$rowData), names(ref_seqs1)), ]
    #     args$referenceSeq <- ref_seqs1
    # # If names of 2nd sequences match with data, add refseq2
    # } else if( length(rownames(args$rowData)) > 0 && !is.null( ref_seqs2 ) &&
    #            all( rownames(args$rowData) %in% names(ref_seqs2) ) ){
    #     # Put sequences into correct order
    #     ref_seqs2 <- ref_seqs2[ match(rownames(args$rowData), names(ref_seqs2)), ]
    #     args$referenceSeq <- ref_seqs2
    # }
    return(args)
}

####################### .check_if_tree_matches_with_data #######################
# This function check if rowTree or colTree matches with the data

# Input: rowTree, rowLinks of both objects, rownames of both objects, 
# rownames of final object
# Output: rowTree and rowLinks as a list if tree matches. Otherwise NULL.
# .check_if_tree_matches_with_data <- function(tree1, tree2, 
#                                              links1, links2, 
#                                              names_original1, names_original2,
#                                              names_final){
#     # Get labels of tree1
#     tree_labels1 <- c( tree1$tip.label, tree1$node.label )
#     tree_labels2 <- c( tree2$tip.label, tree2$node.label )
#     # Get links1
#     node_labs1 <- links1[ , "nodeLab" ]
#     # Otherwise, if links are NULL, Links are names
#     if( is.null(node_labs1) ){
#         node_labs1 <- names_original1
#     }
#     # Add names
#     names(node_labs1) <- names_original1
#     # Get links2
#     node_labs2 <- links2[ , "nodeLab" ]
#     # Otherwise, if links are NULL, Links are names
#     if( is.null(node_labs2) ){
#         node_labs2 <- names_original2
#     }
#     # Add names
#     names(node_labs2) <- names_original2
#     
#     # Combine node labels
#     node_labs <- c(node_labs1, node_labs2)
#     
#     # Check if tree1 matches with data
#     if( length(names_final) > 0 && !is.null( tree1 ) && !is.null(node_labs) &&
#         all( node_labs1 %in% tree_labels1 ) ){
#         # Create a list that combines arguments
#         tree_arguments <- list()
#         # Add tree
#         tree_arguments$tree <- tree1
#         # Put links into correct order
#         node_labs <- node_labs[ match(names_final, names(node_labs)) ]
#         # Add links
#         tree_arguments$node_labs <- node_labs
#     # Check tree2
#     } else if( length(names_final) > 0 && !is.null( tree2 ) && !is.null(node_labs) &&
#                all( node_labs1 %in% tree_labels2 ) ){
#         # Create a list that combines arguments
#         tree_arguments <- list()
#         # Add tree
#         tree_arguments$tree <- tree2
#         # Put links into correct order
#         node_labs <- node_labs[ match(names_final, names(node_labs)) ]
#         # Add links
#         tree_arguments$node_labs <- node_labs
#     # If trees are not NULL
#     } else if( !is.null(tree1) && !is.null(tree2) ){
#         # Create a list that combines arguments
#         tree_arguments <- list()
#         tree_arguments$tree <- list(tree1, tree2)
#         names(tree_arguments$tree) <- c(row_tree_name1, row_tree_name2)
#     # If the trees do not match, give NULL
#     } else{
#         tree_arguments <- NULL
#     }
#     return(tree_arguments)
# }

.check_if_tree_matches_with_data <- function(trees1, trees2,
                                             links1, links2,
                                             names_original1, names_original2,
                                             names_final){
    
    # If there are no trees, labels are NULL
    if( is.null(trees1) ){
        tree_labs1 <- NULL
    } else{
        # Otherwise get labels from tree
        tree_labs1 <- lapply(trees1, FUN = function(x){
            labs <- c( x$tip.label, x$node.label )
            return(labs)
        })
    }
    # If there are no trees, labels are NULL
    if( is.null(trees2) ){
        tree_labs2 <- NULL
    } else{
        # Otherwise get labels from tree
        tree_labs2 <- lapply(trees2, FUN = function(x){
            labs <- c( x$tip.label, x$node.label )
            return(labs)
        })
    }
    # Combine tree labs
    tree_labs <- c(tree_labs1, tree_labs2)
    # Combine trees
    trees <- c(trees1, trees2)
    # Make names uniwue
    names(tree_labs) <- names(trees) <- make.unique(names(trees))
    
    
    # Get links1
    node_labs1 <- links1[ , "nodeLab" ]
    # Otherwise, if links are NULL, Links are names
    if( is.null(node_labs1) ){
        node_labs1 <- names_original1
    }
    # Add names
    names(node_labs1) <- names_original1
    # Get links2
    node_labs2 <- links2[ , "nodeLab" ]
    # Otherwise, if links are NULL, Links are names
    if( is.null(node_labs2) ){
        node_labs2 <- names_original2
    }
    # Add names
    names(node_labs2) <- names_original2
    
    # Unlist node_labs to have one large vector
    node_labs <- unlist(node_labs)
    # Check whether nodes can be found from trees
    result <- lapply(tree_labs, FUN = function(x){
        c( node_labs %in% x )
    })
    # Create a data.frame
    result <- as.data.frame(result)
    rownames(result) <- node_labs
    
    # If there are rows that cannot be linked with trees
    if( any( rowSums(result) == 0 ) ){
        args <- list()
        args$tree <- NULL
        args$links <- NULL
        return(args)
    }
    
    how_many_rows_can_be_found <- colSums(result)
    
    for( i in seq_len(ncol(result)) ){
        combinations <- combn(result, i, simplify = FALSE)
        res <- lapply(combinations, FUN = function(x){
            all( rowSums(x) == ncol(x) )
        })
        res <- unlist(res)
        if( any(res) ){
            # Take the first combination that have all the rows
            combinations <- combinations[[which(res)[[1]]]]
            # Take the tree name
            tree_names <- colnames(combinations)
            break
        }
    }
    
    result <- result[ , tree_names, drop = FALSE ]
    whichTree <- apply(result, 1, FUN = function(x){
        names(result)[x == TRUE]
    }
        )
    
    rowLinks <- data.frame(
        nodeLab = node_labs,
        whichTree = whichTree
        )
    
    temp <- apply(rowLinks, 1, FUN = function(x){
        temp <- c()
        nodeNum <- convertNode(trees[[x[["whichTree"]]]], x[["nodeLab"]])
        nodeLab_alias <- paste0("alias_", nodeNum)
        isLeaf <- isLeaf(trees[[x[["whichTree"]]]], x[["nodeLab"]])
        
        return(c(nodeNum = nodeNum, nodeLab_alias = nodeLab_alias, isLeaf = isLeaf))
    })
    temp <- t(temp)
    colnames(temp) = c("nodeNum", "nodeLab_alias", "isLeaf")
    rowLinks <- cbind(rowLinks, temp)
    rowLinks <- rowLinks[ , c("nodeLab", "nodeNum", "nodeLab_alias", "isLeaf", "whichTree")]
    rowLinks <- LinkDataFrame(
        nodeLab = rowLinks[["nodeLab"]],
        nodeNum = as.integer(rowLinks[["nodeNum"]]),
        nodeLab_alias = rowLinks[["nodeLab_alias"]],
        isLeaf = as.logical(rowLinks[["isLeaf"]]),
        whichTree = rowLinks[["whichTree"]]
    )
    # If there are rows that cannot be linked with any of 
}

######################## .merge_SingleCellExperiments ##########################
# This function merges the data of two SCE objects into one set of arguments that
# can be feed to create a single object.
# Input: Two SCEs
# Output: A single SCE
.merge_SingleCellExperiments <- function(tse_original, tse, join,  
                                         assay_name, missing_values){
    # reducedDim is additional slot for SCE compared to SE. 
    # However, merging reducedDims leads to non-meaningful data
    # Merge data to get a list of arguments
    args <- .merge_SummarizedExperiments(tse_original, tse, join,
                                         assay_name, missing_values)
    return(args)
}

######################## .merge_SummarizedExperiments ##########################
# This function merges the data of two SE objects into one set of arguments that
# can be feed to create a single object.

# Input: Two SEs
# Output: A list of arguments
.merge_SummarizedExperiments <- function(tse_original, tse, join,  
                                         assay_name, missing_values){
    # Merge rowData
    rowdata <- .merge_rowdata(tse_original, tse, join)
    # Merge colData
    coldata <- .merge_coldata(tse_original, tse, join)
    # Merge assay
    assay <- .merge_assay(tse_original, tse, assay_name, join, missing_values, rowdata, coldata)
    assays <- SimpleList(name = assay)
    names(assays) <- assay_name
    # Combine metadata
    metadata <- c( metadata(tse_original), metadata(tse) )
    
    # Create a list of data
    args <- list(assays = assays,
                rowData = rowdata,
                colData = coldata,
                metadata = metadata)
    return(args)
}

################################ .merge_assay ##################################
# This function merges assays.

# Input: Two TreeSEs, the name of the assay, joining method, value to denote
# missing values, merged rowData, and merged colData
# Output: Merged assay
.merge_assay <- function(tse_original, tse, assay_name, join,
                         missing_values, rd, cd){
    # Take assays
    assay1 <- assay(tse_original, assay_name)
    assay2 <- assay(tse, assay_name)
    
    # Merge two assays into one
    assay <- .join_two_tables(assay1, assay2, join)
    
    # Fill missing values
    assay[ is.na(assay) ] <- missing_values
    # Convert into matrix
    assay <- as.matrix(assay)
    
    # Order the assay based on rowData and colData
    assay <- assay[ match(rownames(rd), rownames(assay)), , drop = FALSE ]
    assay <- assay[ , match(rownames(cd), colnames(assay)), drop = FALSE]
    
    return(assay)
}

############################### .merge_rowdata #################################
# This function merges rowDatas,

# Input: Two TreeSEs and joining method
# Output: Merged rowData
.merge_rowdata <- function(tse_original, tse, join){
    # Take rowDatas
    rd1 <- rowData(tse_original)
    rd2 <- rowData(tse)
    
    # Convert column names to lower
    if( length(colnames(rd1)) > 0 ){
        colnames(rd1) <- tolower(colnames(rd1))
    }
    if( length(colnames(rd2)) > 0 ){
        colnames(rd2) <- tolower(colnames(rd2))
    }
    
    # Merge rowdata
    rd <- .join_two_tables(rd1, rd2, join)
    
    # Get column indices that match with taxonomy ranks
    ranks_ind <- match( TAXONOMY_RANKS, colnames(rd) )
    # Remove NAs
    ranks_ind <- ranks_ind[ !is.na(ranks_ind) ]
    # If ranks were found
    if( length(ranks_ind) != 0 ){
        # Get the data in correct order, take only column that have ranks
        rd_rank <- rd[ , ranks_ind, drop = FALSE]
        # Take other columns
        rd_other <- rd[ , !ranks_ind, drop = FALSE]
        # Get rank names
        rank_names <- colnames(rd_rank)
        # Convert names s that they have capital letters
        new_rank_names <- paste(toupper(substr(rank_names, 1, 1)), 
                                substr(rank_names, 2, nchar(rank_names)), sep = "")
        # Add new names to colnames of rd_rank
        colnames(rd_rank) <- new_rank_names
        
        # Combine columns
        rd <- cbind(rd_rank, rd_other)
    }
    return(rd)
}

############################### .merge_coldata #################################
# This function merges colDatas,

# Input: Two TreeSEs and joining method
# Output: Merged colData
.merge_coldata <- function(tse_original, tse, join){
    # Take colDatas
    cd1 <- colData(tse_original)
    cd2 <- colData(tse)
    
    # Merge coldata
    cd <- .join_two_tables(cd1, cd2, join = "full")
    # Convert into DataFrame
    cd <- DataFrame(cd)
    
    return(cd)
}

############################## .join_two_tables ################################
# This general function is used to merge rowDatas, colDatas, and assays.

# Input: Two tables and joining method
# Output: One merged table

#' @importFrom dplyr coalesce
.join_two_tables <- function(df1, df2, join){
    # Get parameter based on join
    all.x <- switch(join,
                    full = TRUE,
                    inner = FALSE,
                    left = TRUE,
                    right = FALSE
    )
    all.y <- switch(join,
                    full = TRUE,
                    inner = FALSE,
                    left = FALSE,
                    right = TRUE
    )
    # Ensure that the data is in correct format
    df1 <- as.data.frame(df1)
    df2 <- as.data.frame(df2)
    # Order tables in alphabetical order
    df1 <- df1[ , order(colnames(df1)) ]
    df2 <- df2[ , order(colnames(df2)) ]
    
    # Get matching variables indices
    matching_variables_ids1 <- match( colnames(df2), colnames(df1) )
    # Get matching variable names
    matching_variables1 <- colnames(df1)[ matching_variables_ids1 ]
    # Remove NAs
    matching_variables1 <- matching_variables1[ !is.na(matching_variables1) ]
    
    # Get matching variables indices
    matching_variables_ids2 <- match( colnames(df1), colnames(df2) )
    # Get matching variable names
    matching_variables2 <- colnames(df2)[ matching_variables_ids2 ]
    # Remove NAs
    matching_variables2 <- matching_variables2[ !is.na(matching_variables2) ]
    
    # Make the matching variables unique
    matching_variables_mod1 <- paste0(matching_variables1, "_X")
    matching_variables_ids1 <- matching_variables_ids1[ !is.na(matching_variables_ids1) ]
    colnames(df1)[ matching_variables_ids1 ] <- matching_variables_mod1
    matching_variables_mod2 <- paste0(matching_variables2, "_Y")
    matching_variables_ids2 <- matching_variables_ids2[ !is.na(matching_variables_ids2) ]
    colnames(df2)[ matching_variables_ids2 ] <- matching_variables_mod2
    
    # Add rownames to one of the columns
    df1$rownames_merge_ID <- rownames(df1)
    df2$rownames_merge_ID <- rownames(df2)
    # Merge data frames into one data frame
    df <- merge(df1, df2, by = "rownames_merge_ID", all.x = all.x, all.y = all.y)
    # Add rownames and remove additional column
    rownames(df) <- df$rownames_merge_ID
    df$rownames_merge_ID <- NULL
    
    # Combine matching variables if found
    if( length(matching_variables1) > 0 ){
        for(i in 1:length(matching_variables1) ){
            # Get columns
            x <- matching_variables_mod1[i]
            y <- matching_variables_mod2[i]
            # Combine information from columns
            x_and_y_combined <- coalesce( df[ , x], df[ , y] )
            # Remove additional columns
            df[ , x ] <- NULL
            df[ , y ] <- NULL
            # Add column that has combined information
            df[ , matching_variables1[i] ] <- x_and_y_combined
        }
    }
    return(df)
}
