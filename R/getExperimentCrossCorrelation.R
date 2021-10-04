#' Calculate cross-correlation
#' 
#' @param x A
#'   \code{\link[MultiAssayExperiment:MultiAssayExperiment-class]{MultiAssayExperiment}} or
#'   \code{\link[TreeSummarizedExperiment:TreeSummarizedExperiment-class]{TreeSummarizedExperiment}}
#'   object.
#' @param experiment1 ....
#' @param experiment2 .........
#' @param abund_values1 A single character value for selecting the
#'   \code{\link[SummarizedExperiment:SummarizedExperiment-class]{assay}} to be
#'   transformed.
#' @param abund_values2 A single character value for selecting the
#'   \code{\link[SummarizedExperiment:SummarizedExperiment-class]{assay}} to be
#'   transformed.
#' @param method ..............
#' @param mode ................
#' @param p_adj_method .............
#' @param p_adj_threshold ......................
#' @param n_signif .......................
#' @param cth ................
#' @param order ........................
#' @param filter_self_correlations...................
#' 
#' @details
#' Calculate cross-correlation.
#'
#' @references
#' Add references here.
#'
#' @return 
#' Matrix or table.
#'
#' @name getExperimentCrossCorrelation
#' @export
#'
#' @author Leo Lahti and Tuomas Borman. Contact: \url{microbiome.github.io}
#'
#' @examples
#' mae <- microbiomeDataSets::peerj32()
#' 
#' # Subset so that less observations / quicker to run, just for example
#' mae[[1]] <- mae[[1]][1:20, 1:10]
#' mae[[2]] <- mae[[2]][1:20, 1:10]
#' # Calculate cross-correlations
#' result <- getExperimentCrossCorrelation(mae, method = "pearson")
#' # Show first 5 entries
#' head(result, 5)
#' 
#' # Same can be done with TreeSummarizedExperiment and altExp
#' # Create TreeSE with altExp
#' tse <- mae[[1]]
#' altExp(tse, "experiment2") <- mae[[2]]
#' result <- getExperimentCrossCorrelation(tse, experiment2 = "experiment2", method = "pearson")
#' # Show first 5 entries
#' head(result, 5)
#' 
NULL

#' @rdname getExperimentCrossCorrelation
#' @export
setGeneric("getExperimentCrossCorrelation", signature = c("x"),
           function(x,
                    experiment1 = 1,
                    experiment2 = 2,
                    abund_values1 = "counts",
                    abund_values2 = "counts",
                    method = "spearman",
                    mode = "table",
                    p_adj_method = "fdr",
                    p_adj_threshold = 0.05,
                    n_signif = 0,
                    cth = NULL,
                    order = FALSE,
                    filter_self_correlations = FALSE)
             standardGeneric("getExperimentCrossCorrelation"))

#' @rdname getExperimentCrossCorrelation
#' @export
setMethod("getExperimentCrossCorrelation", signature = c(x = "MultiAssayExperiment"),
    function(x,
             experiment1 = 1,
             experiment2 = 2,
             abund_values1 = "counts",
             abund_values2 = "counts",
             method = c("categorical", "pearson","spearman"),
             mode = "table",
             p_adj_method = "fdr",
             p_adj_threshold = 0.05,
             n_signif = 0,
             cth = NULL,
             order = FALSE,
             filter_self_correlations = FALSE){
        ############################# INPUT CHECK ##############################
        # Check experiment1 and experiment2
        if( is.character(experiment1) && !experiment1 %in% names(experiments(x)) || 
            is.numeric(experiment1) && !experiment1 <= length(experiments(x)) ){
          stop("'experiment1' must be numeric or character value specifying", 
               " experiment in experiment(x).", call. = FALSE)
        }
        if( is.character(experiment2) && !experiment2 %in% names(experiments(x)) || 
            is.numeric(experiment2) && !experiment2 <= length(experiments(x)) ){
          stop("'experiment2' must be numeric or character value specifying", 
               " experiment in experiment(x).", call. = FALSE)
        }
        # Fetch tse objects
        tse1 <- mae[[experiment1]]
        tse2 <- mae[[experiment2]]
        # # If object2 is not specified, then correlate object1 with object1
        # if(is.null(experiment2)){
        #   tse2 <- mae[[experiment1]]
        # } else{
        #   tse2 <- mae[[experiment2]]
        # }
        # 
        # Check abund_values1 and abund_values2
        .check_assay_present(abund_values1, tse1)
        .check_assay_present(abund_values2, tse2)
        # Check method
        # If method is not single string, user has not specified method,
        # or has given e.g. a vector
        if(!.is_non_empty_string(method)){
          stop("'method' must be a non-empty single character value.",
               call. = FALSE)
        }
        method <- match.arg(method)
        # Check mode
        if( !.is_non_empty_string(mode) && !mode %in% c("matrix", "table") ){
          stop("'mode' must be 'matrix' or 'table'.", call. = FALSE)
        }
        # p_adj_method is checked in p.adjust
        # Check p_adj_threshold
        if( !is.numeric(p_adj_threshold) && p_adj_threshold<0 || p_adj_threshold>1 ){
          stop("'p_adj_threshold' must be a numeric value >= 0.", call. = FALSE)
        }
        # Check n_signif
        if (!is.numeric(n_signif) && n_signif > 0 ){
          stop("'p_adj_threshold' must be a numeric value greater than or equal to 0.", 
               call. = FALSE)
        }
        # Check cth
        if( ifelse( is.numeric(cth), !(cth>=0 && cth<=1), !is.null(cth)) ){
          stop("'cth' must be a numeric value greater than or equal to 0.", 
               call. = FALSE)
        }
        # Check order
        if( !(order == TRUE || order == FALSE) ){
          stop("'order' must be a boolean value.", 
               call. = FALSE)
        }
        # Check filter_self_correlations
        if( !(filter_self_correlations == TRUE || filter_self_correlations == FALSE) ){
          stop("'filter_self_correlations' must be a boolean value.", 
               call. = FALSE)
        }
        ############################ INPUT CHECK END ###########################
        # Fetch assays to correlate
        assay1 <- assay(tse1, abund_values1)
        assay2 <- assay(tse2, abund_values1)
        # Transposes tables to right format
        assay1 <- t(assay1)
        assay2 <- t(assay2)
        # Converts tables to data frame
        assay1 <- as.data.frame(assay1)
        assay2 <- assay2
        
        # Get assay in right format
        assay1 <- .cor_test_data_type(assay1, method)
        # Calculate correlations
        results <- .calculate_correlation(assay1, assay2, method, p_adj_method)
        correlations <- results$correlations
        p_values <- results$p_values
        p_values_adjusted <- results$p_values_adjusted
        # Do filtering
        result <- .correlation_filter(correlations, 
                                      p_values, 
                                      p_values_adjusted, 
                                      p_adj_threshold,
                                      n_signif, 
                                      cth,
                                      order,
                                      assay1, assay2, filter_self_correlations)
        # Matrix or table?
        if (mode == "table") {
          result <- .correlation_matrix_to_table(result)
        }
        return(result)
    }
)

#' @rdname getExperimentCrossCorrelation
#' @export
setMethod("getExperimentCrossCorrelation", signature = c(x = "TreeSummarizedExperiment"),
    function(x,
             experiment2 = NULL,
             abund_values1 = "counts",
             abund_values2 = "counts",
             method = "spearman",
             mode = "table",
             p_adj_method = "fdr",
             p_adj_threshold = 0.05,
             n_signif = 0,
             cth = NULL,
             order = FALSE,
             filter_self_correlations = FALSE){
        ############################## INPUT CHECK #############################
        if( is.character(experiment2) && !experiment2 %in% names(experiments(x)) || 
            is.numeric(experiment2) && !experiment2 <= length(experiments(x)) ){
          stop("'experiment2' must be numeric or character value specifying", 
               " experiment in altExps(x).", call. = FALSE)
        }
        ############################ INPUT CHECK END ###########################
        # Fetch data sets and create a MAE object
        exp1 <- x
        exp2 <- altExp(x, experiment2)
        x <- MultiAssayExperiment::MultiAssayExperiment(experiments = ExperimentList(exp1 = exp1, exp2 = exp2))
        # Call method with MAE object as an input
        getExperimentCrossCorrelation(x,
                                      experiment1 = 1,
                                      experiment2 = 2,
                                      abund_values1,
                                      abund_values2,
                                      method,
                                      mode,
                                      p_adj_method,
                                      p_adj_threshold,
                                      n_signif,
                                      cth,
                                      order,
                                      filter_self_correlations)
    }
)

################################ HELP FUNCTIONS ################################
.cor_test_data_type <- function(assay, method){
  # Different available methods
  numeric_methods <- c("spearman", "pearson")
  categorical_methods <- c("categorical")
  # Check if method match with values, otherwise give warning.
  # For numeric methods, get only numeric values. For categorical methods, get only factors.
  if (method %in% numeric_methods) {
    inds <- vapply(assay, is.numeric, TRUE)
    # If there are no numeric values, give an error
    if( all(!inds) ){
      stop("Assay, specified by 'abund_values', of 'experiment1' does not include",
           " numeric values. Choose categorical method for 'method'.",
           call. = FALSE)
    }
    else if (any(!inds)) {
      warning("Considering only numeric annotations for \n       
                    pearson/spearman")
    }
  } else if (method %in% categorical_methods) {
    inds <- vapply(assay, is.factor, TRUE)
    # If there are no factor values, give an error
    if( all(!inds) ){
      stop("Assay, specified by 'abund_values', of 'experiment1' does not include",
           " factor values. Choose numeric method for 'method'.",
           call. = FALSE)
    }
    else if (any(!inds)) {
      warning("Considering only categorical annotations for factors")
    }
  }
  # Names of features that are numeric or factor (based on previous step)
  names <- names(which(inds))
  # Subset assay to get only numeric or factor values
  if (!is.vector(assay)) {
    assay <- suppressWarnings(as.matrix(assay[, inds], ncol=length(inds)))
  } else {
    assay <- as.matrix(assay[inds], ncol=length(inds))
  }
  # Add colnames to assay
  colnames(assay) <- names
  return(assay)
}

.calculate_correlation <- function(assay1, assay2, method, p_adj_method){
  # Create empty matrices
  correlations <- matrix(NA, ncol(assay1), ncol(assay2))
  rownames(correlations) <- colnames(assay1)
  colnames(correlations) <- colnames(assay2)
  p_values <- correlations
  
  # Calculate correlations, different methods for numeric and categorical data
  if (method %in% c("pearson", "spearman")) {
    # Specify minimum number of observation needed to do correlation analysis
    minimum_number_of_observations <- 8
    # Loop over every feature in assay2
    for (j in seq_len(ncol(assay2))) {
      # Loop over every feature in assay1
      jc <- apply(assay1, 2, function(xi) {
        # If there are enough observations/samples, do correlation test, and store the result
        # to temporary object
        if (sum(!is.na(xi)) >= minimum_number_of_observations) {
          res <- suppressWarnings(
            cor.test(xi, unlist(assay2[, j], use.names=FALSE), 
                     method=method, use="pairwise.complete.obs"))
          res <- c(res$estimate, res$p.value)
        } 
        # If there are not enough observations/samples, give a warning, and store
        # just NAs
        else {
          warning(paste("Not enough observations (",
                        minimum_number_of_observations, "required); \n   
                        (", 
                        sum(!is.na(xi)), ") \n \n 
                        - skipping correlation estimation"))
          res <- c(NA, NA)
        }
        # Return res which will be stored in jc object
        res
      })
      
      # 'jc' object includes all the correlation between individual feature from assay2
      # and all the features from assay1
      # Finally, store correlations and continue loop with next feature of assay2
      correlations[, j] <- jc[1, ]
      p_values[, j] <- jc[2, ]
    }
  } 
  # If method is categorical
  else if (method == "categorical") {
    # Loop through all the features from assay1
    for (varname in colnames(assay1)) {
      # Loop through all the features from assay2
      for (lev in colnames(assay2)) {
        # Get all the values of individual features
        xvec <- assay1[, varname]
        yvec <- assay2[, lev]
        
        # Keep only those samples that have values in both features
        keep <- rowSums(is.na(cbind(xvec, yvec))) == 0
        xvec <- xvec[keep]
        yvec <- yvec[keep]
        
        # Number of data-annotation samples for
        # calculating the correlations
        n <- sum(keep)
        # Calculate cross-correlation using Goorma and Kruskal tau
        correlations[varname, lev] <- .calculate_gktau(xvec, yvec) 
      }
    }
  }
  
  # If there are p_values that are not NA, adjust them
  if (!all(is.na(p_values))) {
    # Corrected p-values
    p_values_adjusted <- array(NA, dim=dim(p_values))
    p_values_adjusted <- matrix(p.adjust(p_values, method=p_adj_method), nrow=nrow(p_values))
    dimnames(p_values_adjusted) <- dimnames(p_values)
  }
  
  return(list(correlations = correlations, 
              p_values = p_values, 
              p_values_adjusted = p_values_adjusted))
}

.correlation_filter <- function(correlations,
                                p_values,
                                p_values_adjusted,
                                p_adj_threshold,
                                n_signif,
                                cth, 
                                order,
                                assay1, assay2, filter_self_correlations){
  # Filter
  if (!is.null(p_adj_threshold) || !is.null(cth)) {
    
    # Replace NAs with extreme values for filtering purposes
    p_values_adjusted[is.na(p_values_adjusted)] <- 1
    p_values[is.na(p_values_adjusted)] <- 1
    correlations[is.na(correlations)] <- 0
    
    # Filter by adjusted pvalues and correlations
    inds1.q <- inds2.q <- inds1.c <- inds2.c <- NULL
    
    # Which features have more adjusted p-values under the threshold than the 
    # 'n_signif' specifies
    if (!is.null(p_adj_threshold)) {
      
      inds1.q <- apply(p_values_adjusted, 1, function(x) {
        sum(x < p_adj_threshold) >= n_signif
      })
      
      inds2.q <- apply(p_values_adjusted, 2, function(x) {
        sum(x < p_adj_threshold) >= n_signif
      })
    }
    
    # Which features have correlation over correlation threshold?
    if (!is.null(cth)) {
      inds1.c <- apply(abs(correlations), 1, function(x) {
        sum(x > cth | x < cth) >= n_signif
      })
      inds2.c <- apply(abs(correlations), 2, function(x) {
        sum(x > cth | x < cth) >= n_signif
      })
    }
    
    # Combine results from previous step
    if (!is.null(p_adj_threshold) && !is.null(cth)) {
      inds1 <- inds1.q & inds1.c
      inds2 <- inds2.q & inds2.c
    } else if (is.null(p_adj_threshold) && !is.null(cth)) {
      inds1 <- inds1.c
      inds2 <- inds2.c
    } else if (!is.null(p_adj_threshold) && is.null(cth)) {
      inds1 <- inds1.q
      inds2 <- inds2.q
    }
    
    # TODO: add also correlation filter, not only significance
    # Require each has at least n_signif. correlations
    
    # If both features have more TRUEs than n_signif specifies /
    # if there are significant correlations
    if (sum(inds1) >= n_signif && sum(inds2) >= n_signif) {
      # Get names of those features that were TRUE
      rnams <- rownames(correlations)[inds1]
      cnams <- colnames(correlations)[inds2]
      # Subset correlations, p_values, and adjusted p_values table
      correlations <- matrix(correlations[inds1, inds2, drop=FALSE], nrow=sum(inds1))
      p_values <- matrix(p_values[inds1, inds2, drop=FALSE], nrow=sum(inds1))
      p_values_adjusted <- matrix(p_values_adjusted[inds1, inds2, drop=FALSE], nrow=sum(inds1))
      # Add row and column names
      rownames(p_values_adjusted) <- rownames(p_values) <- rownames(correlations) <- rnams
      colnames(p_values_adjusted) <- colnames(p_values) <- colnames(correlations) <- cnams
      
      # If order was specified and there is more than 1 feature left in both feature sets
      if (order && sum(inds1) >= 2 && sum(inds2) >= 2) {
        
        # Order in visually appealing order
        tmp <- correlations
        rownames(tmp) <- NULL
        colnames(tmp) <- NULL
        # Do hierarchical clustering
        rind <- hclust(as.dist(1 - cor(t(tmp),
                                       use="pairwise.complete.obs")))$order
        cind <- hclust(as.dist(1 - cor(tmp,
                                       use="pairwise.complete.obs")))$order
        # Get the order of features from hiearchical clustering
        rnams <- rownames(correlations)[rind]
        cnams <- colnames(correlations)[cind]
        
        # Order the tables based on order of hiearchical clustering
        correlations <- correlations[rind, cind]
        p_values <- p_values[rind, cind]
        p_values_adjusted <- p_values_adjusted[rind, cind]
        # Add column and rownames
        rownames(p_values_adjusted) <- rownames(p_values) <- rownames(correlations) <- rnams
        colnames(p_values_adjusted) <- colnames(p_values) <- colnames(correlations) <- cnams
      }
    } 
    # If there were no significant correlations, give a message
    else {
      message("No significant correlations with the given criteria\n")
      correlations <- p_values <- p_values_adjusted <- NULL
    }
  }
  res <- list(cor=correlations, pval=p_values, p.adj=p_values_adjusted)
  # Filter self correlations if it's specified
  if (nrow(assay1) == nrow(assay2) && ncol(assay1) == ncol(assay2) && filter_self_correlations) {
    diag(res$cor) <- diag(res$pval) <- diag(res$p.adj) <- NA
  }
  return(res)
}

.correlation_matrix_to_table <- function(res) {
  # Melt correlation table
  ctab <- ID <- NULL
  if (!is.null(res$cor)) {
    ctab <- as.data.frame(res$cor)
    ctab$ID <- rownames(res$cor)
    ctab <- reshape2::melt(ctab, "ID")
    colnames(ctab) <- c("X1", "X2", "Correlation")
    ctab$Correlation <- as.numeric(as.character(ctab$Correlation))
  
    # Melt p-values and add them to the melted correlation table
    # If there are adjusted p_values
    if (!is.null(res$p.adj)) {
      ctab2 <- as.data.frame(res$p.adj)
      ctab2$ID <- rownames(res$p.adj)
      ctab2 <- reshape2::melt(ctab2, "ID")
      colnames(ctab2) <- c("X1", "X2", "p.adj")
      ctab2$p.adj <- as.numeric(as.character(ctab2$p.adj))
      
      ctab <- cbind(ctab, ctab2$p.adj)
      colnames(ctab) <- c("X1", "X2", "Correlation", "p.adj")
      ctab <- ctab[order(ctab$p.adj), ]
      colnames(ctab) <- c("X1", "X2", "Correlation", "p.adj")
      
    } 
    # If there are only p-values that are not adjusted
    else {
        ctab2 <- as.data.frame(res$pval)
        ctab2$ID <- rownames(res$pval)
        ctab2 <- reshape2::melt(ctab2, "ID")
        colnames(ctab2) <- c("X1", "X2", "value")
        ctab2$value <- as.numeric(as.character(ctab2$value))
        
        ctab <- cbind(ctab, ctab2$value)
        ctab <- ctab[order(-abs(ctab$Correlation)), ]
        colnames(ctab) <- c("X1", "X2", "Correlation", "pvalue")
    }
  
    # Convert feature names into characters
    ctab$X1 <- as.character(ctab$X1)
    ctab$X2 <- as.character(ctab$X2)
    # Keep the original order of factor levels
    ctab$X1 <- factor(as.character(ctab$X1), levels=rownames(res$cor))
    ctab$X2 <- factor(as.character(ctab$X2), levels=colnames(res$cor))
    # Remove NAs
    ctab <- ctab[!is.na(ctab$Correlation), ]
    
    # Order the table by p-value
    if ("p.adj" %in% colnames(ctab)) {
      ctab <- ctab[order(ctab$p.adj), ]
    } else if ("pvalue" %in% colnames(ctab)) {
      ctab <- ctab[order(ctab$pvalue), ]
    }
  }
  ctab
  
}

.calculate_gktau <- function(x, y){
  # First, compute the IxJ contingency table between x and y
  Nij <- table(x, y, useNA="ifany")
  # Next, convert this table into a joint probability estimate
  PIij <- Nij/sum(Nij)
  # Compute the marginal probability estimates
  PIiPlus <- apply(PIij, MARGIN=1, sum)
  PIPlusj <- apply(PIij, MARGIN=2, sum)
  # Compute the marginal variation of y
  Vy <- 1 - sum(PIPlusj^2)
  # Compute the expected conditional variation of y given x
  InnerSum <- apply(PIij^2, MARGIN=1, sum)
  VyBarx <- 1 - sum(InnerSum/PIiPlus)
  # Compute and return Goodman and Kruskal's tau measure
  tau <- (Vy - VyBarx)/Vy
  tau
}
