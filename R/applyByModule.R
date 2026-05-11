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
#' @param meta.name \code{Character scalar}. The name of the metadata slot where
#'   results are stored. (Default: \code{"mod.res"})
#' 
#' @returns \code{x} updated with results stored in the appropriate slot
#' depending on the given \code{FUN}.
#' 
#' @examples
#' library(ariadne)
#' library(miaViz)
#' 
#' # Import dataset
#' data("Tengeler2020", package = "mia")
#' tse <- Tengeler2020
#' 
#' # Create taxon labels
#' tax.labs <- getTaxonomyLabels(tse, make.unique = FALSE)
#' tax.labs <- sub("^.+:", "", tax.labs)
#' rowData(tse)$taxname <- tax.labs
#' 
#' # Import ariadne graph
#' graph <- ariadne()
#' 
#' # Link taxa to BugSigDB modules
#' tax2bugsig <- weavePath(graph, taxname ~ bugsig, k = 3, init = tax.labs)
#' 
#' # Add modules to TreeSE
#' tse <- addModules(tse, tax2bugsig, key = "taxname")
#' 
#' # Select a subset of modules to analyse
#' mod.names <- levels(tax2bugsig$bugsig)
#' 
#' # Agglomerate experiment by module
#' altExp(tse, "modules") <- agglomerateByModule(
#'     tse,
#'     by = "rows",
#'     group = mod.names
#' )
#' 
#' # Get names of top modules
#' top.mods <- getTop(altExp(tse, "modules"))
#' 
#' # Compute alpha diversity indices for each module
#' tse <- applyByModule(
#'     tse,
#'     by = "rows",
#'     group = top.mods,
#'     FUN = getAlpha,
#'     index = c("shannon", "faith")
#' )
#' 
#' # Compute beta diversity measures for each module
#' tse <- applyByModule(
#'     tse,
#'     by = "rows",
#'     group = top.mods,
#'     FUN = getMDS,
#'     method = "unifrac"
#' )
#' 
#' # Transform assay for each module, storing results as altExps
#' tse <- applyByModule(
#'     tse,
#'     by = "rows",
#'     group = top.mods,
#'     FUN = transformAssay,
#'     method = "clr",
#'     pseudocount = TRUE
#' )
#' 
#' # Plot abundance for each module, specifying metadata name
#' tse <- applyByModule(
#'     tse,
#'     by = "rows",
#'     group = top.mods,
#'     FUN = plotAbundance,
#'     meta.name = "abund_plots"
#' )
#' 
#' # Plot phylogeny for each module
#' tse <- applyByModule(
#'     tse,
#'     by = "rows",
#'     group = top.mods,
#'     FUN = plotRowTree,
#'     edge.colour.by = "Genus",
#'     add.legend = TRUE,
#'     meta.name = "tree_plots"
#' )
NULL

#' @export
#' @rdname applyByModule
#' @importFrom dplyr bind_cols
#' @importFrom BiocParallel bplapply
applyByModule <- function(x, by, group, FUN, ..., meta.name = "mod.res"){
    # Check margin
    if( length(by) == 0L || !by %in% c("rows", "cols") ){
        stop("'by' must be either 'rows' or 'cols'.", call. = FALSE)
    }
    # Select side information based on margin
    data_fun <- if( by == "rows" ) rowData else colData
    
    if( !all(group %in% colnames(data_fun(x))) ){
        stop("All elements of 'group' must be variables of 'x'.", call. = FALSE)
    }
    
    df <- data_fun(x)[group]
    
    idx <- apply(df, 2L, function(col) which(col != 0))
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
            dplyr::bind_cols(.name_repair = "minimal")
        
        colnames(res) <- paste0(colnames(res), ".", col_labs)
        colData(x) <- cbind(colData(x), res)

    }else if( fun_name %in%  c("getDivergence", "getDominant") ){
      
      if( "name" %in% names(kwargs) ){
          lab <- kwargs$name
      }else{
          lab <- switch(
              fun_name, getDivergence = "divergence", getDominant = "dominant"
          )
      }
      
      res <- res |>
          lapply(as.data.frame) |>
          dplyr::bind_cols(.name_repair = "minimal")
      
      colnames(res) <- paste0(lab, ".", group)
      colData(x) <- cbind(colData(x), res)
        
    }else if( fun_name %in% c("getMDS", "getNMDS", "getCCA", "getRDA", "getLDA",
        "getNMF", "getDPCoA", "calculatePCA") ){
        
        reducedDims(x) <- c(reducedDims(x), res)
        
    }else if( fun_name == "getCluster" ){
        
        res <- res |>
            lapply(as.data.frame) |>
            dplyr::bind_cols(.name_repair = "minimal")
        
        clust_by <- if( "MARGIN" %in% names(kwargs) ) kwargs$MARGIN else "rows"
        clust_col <- if( "clust.col" %in% names(kwargs) ) kwargs$clust.col else "cluster"
        
        colnames(res) <- paste0(clust_col, ".", group)
        
        if( clust_by == "rows" ){
            rowData(x) <- cbind(rowData(x), res)
        }else{
            colData(x) <- cbind(colData(x), res)
        }
        
    }else{
        metadata(x)[[meta.name]] <- res
    }
    
    return(x)
}
