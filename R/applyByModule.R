#' Apply a function by feature or sample group
#' 
#' @name applyByModule
#' 
#' @param x A \code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}} object.
#' 
#' @param by \code{Character scalar}. Whether \code{group} information pertains
#'   the \code{"rows"} or \code{"cols"} of \code{x}.
#' 
#' @param group \code{Character vector}. Names of the groups or modules by which
#'   \code{FUN} is applied.
#' 
#' @param FUN \code{Function scalar}. A function which takes \code{x} as the
#'   first argument and \code{...} as additional arguments. To update \code{x}
#'   in place, use \code{get} and not \code{add} functions, as the latter will
#'   populate the altExp slot.
#' 
#' @param ... Additional arguments passed to \code{FUN}.
#' 
#' @param min.group.size \code{Numeric scalar}. The minimum number of features
#'   or samples in a group for \code{FUN} to be applied. (Default: \code{NULL})
#' 
#' @param meta.name \code{Character scalar}. The name of the metadata slot where
#'   results are stored. (Default: \code{"mod.res"})
#' 
#' @returns \code{x} updated with results stored in the appropriate slot
#' depending on the given \code{FUN}.
#' 
#' @examples
#' library(miaViz)
#' 
#' # Import dataset
#' data("Tengeler2020", package = "mia")
#' tse <- Tengeler2020
#' 
#' # Compute alpha diversity indices for each Order
#' tse <- applyByModule(
#'     tse,
#'     by = "rows",
#'     group = "Order",
#'     FUN = getAlpha,
#'     index = c("shannon", "faith"),
#'     min.group.size = 2
#' )
#' 
#' # Compute beta diversity measures for each Order
#' tse <- applyByModule(
#'     tse,
#'     by = "rows",
#'     group = "Order",
#'     FUN = getMDS,
#'     method = "unifrac",
#'     min.group.size = 2
#' )
#' 
#' # Transform assay for each Order, storing results as altExps
#' tse <- applyByModule(
#'     tse,
#'     by = "rows",
#'     group = "Order",
#'     FUN = transformAssay,
#'     method = "clr",
#'     pseudocount = TRUE
#' )
#'
#' rowData(tse)$Var1 <- rowData(tse)$Family == "Ruminococcaceae"
#' rowData(tse)$Var2 <- rowData(tse)$Order == "Clostridiales"
#' 
#' # Plot abundance for each module, specifying metadata name
#' tse <- applyByModule(
#'     tse,
#'     by = "rows",
#'     group = c("Var1", "Var2"),
#'     FUN = plotAbundance,
#'     meta.name = "abund_plots"
#' )
#' 
NULL

#' @export
#' @rdname applyByModule
#' @importFrom dplyr bind_cols
#' @importFrom BiocParallel bplapply
applyByModule <- function(x, by, group, FUN, ..., min.group.size = NULL,
    meta.name = "mod.res"){
    # Check margin
    if( length(by) == 0L || !by %in% c("rows", "cols") ){
        stop("'by' must be either 'rows' or 'cols'.", call. = FALSE)
    }
    # Select side information based on margin
    data_fun <- if( by == "rows" ) rowData else colData
    
    if( !.is_non_empty_character(group) || !all(group %in% colnames(data_fun(x))) ){
        stop("All elements of 'group' must be variables of 'x'.", call. = FALSE)
    }
    
    if( length(group) == 1L ){
        
        gvar <- data_fun(x)[[group]]
        
        gvar <- as.factor(gvar)
        
        group <- levels(gvar)
        
        gvar <- as.numeric(gvar)
        
        idx <- lapply(seq_along(group), function(g) which(gvar == g))
        
        names(idx) <- group
        
    }else{
        
        df <- data_fun(x)[group]
        
        idx <- apply(df, 2L, function(col) which(col != 0))
    }
    
    if( !is.null(min.group.size) ){
        
        to_keep <- which(lengths(idx, use.names = FALSE) >= min.group.size)
        
        idx <- idx[to_keep]
        group <- group[to_keep]
    }
    
    res <- bplapply(idx, function(i) FUN(x[i, ], ...))
    
    fun_name <- FUN |>
        substitute() |>
        deparse()
    
    kwargs <- list(...)
    
    if( by == "rows" && inherits(res[[1L]], "SummarizedExperiment") ){
        
        altExps(x) <- c(altExps(x), res)
        
    }else if( fun_name == "getAlpha" ){
        
        num_index <- ncol(res[[1L]])
        col_labs <- rep(group, each = num_index)
        
        res <- res |>
            lapply(as.data.frame) |>
            bind_cols(.name_repair = "minimal")
        
        colnames(res) <- paste0(colnames(res), ".", col_labs)
        colData(x) <- .update_side_info(colData(x), res)

    }else if( fun_name %in% c("getDivergence", "getDominant") ){
      
      if( "name" %in% names(kwargs) ){
          lab <- kwargs$name
      }else{
          lab <- switch(
              fun_name, getDivergence = "divergence", getDominant = "dominant"
          )
      }
      
      res <- res |>
          lapply(as.data.frame) |>
          bind_cols(.name_repair = "minimal")
      
      colnames(res) <- paste0(lab, ".", group)
      colData(x) <- .update_side_info(colData(x), res)
        
    }else if( fun_name %in% c("getMDS", "getNMDS", "getCCA", "getRDA", "getLDA",
        "getNMF", "getDPCoA", "calculatePCA") ){
        
        reducedDims(x) <- c(reducedDims(x), res)
        
    }else if( fun_name == "getCluster" ){
        
        res <- res |>
            lapply(as.data.frame) |>
            bind_cols(.name_repair = "minimal")
        
        clust_by <- if( "MARGIN" %in% names(kwargs) ) kwargs$MARGIN else "rows"
        clust_col <- if( "clust.col" %in% names(kwargs) ) kwargs$clust.col else "cluster"
        
        colnames(res) <- paste0(clust_col, ".", group)
        
        if( clust_by == "rows" ){
            rowData(x) <- .update_side_info(rowData(x), res)
        }else{
            colData(x) <- .update_side_info(colData(x), res)
        }
        
    }else{
        metadata(x)[[meta.name]] <- res
    }
    
    return(x)
}


# Define function to add new cols to side information
.update_side_info <- function(old, new){
    # Find any duplicate columns
    to_remove <- which(colnames(old) %in% colnames(new))
    # Replace duplicate columns
    if( length(to_remove) != 0L ){
        warning("Some columns were replaced.", call. = FALSE)
        old[to_remove] <- NULL
    }
    # Combine old and new columns
    updated <- cbind(old, new)
    return(updated)
}

