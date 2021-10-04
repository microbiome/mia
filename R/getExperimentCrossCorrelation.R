#' Calculate cross-correlation
#' 
#' @param x a
#'   \code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#'   object
#'
#' @details
#'
#' @references
#'
#' @return 
#'
#' @name getExperimentCrossCorrelation
#' @export
#'
#' @author Leo Lahti and Tuomas Borman. Contact: \url{microbiome.github.io}
#'
#' @examples
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
                    verbose = TRUE,
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
        # INPT CHECK
      
        # INPUT CHECK END
        
        # Fetch tse objects
        tse1 <- mae[[experiment1]]
        # If object2 is not specified, then correlate object1 with object1
        if(is.null(experiment2)){
          tse2 <- mae[[experiment1]]
        } else{
          tse2 <- mae[[experiment2]]
        }
        # Fetch assays to correlate
        assay1 <- assay(tse1, abund_values1)
        assay2 <- assay(tse2, abund_values1)
        # Transposes tables to right format
        assay1 <- t(assay1)
        assay2 <- t(assay2)
        # Converts tables to data frame
        assay1 <- as.data.frame(assay1)
        assay2 <- assay2
        # names of features
        feature_names1 <- colnames(assay1)
        feature_names2 <- colnames(assay2)
        
        assay1 <- .cor_test_data_type(assay1, method)
        
        results <- .calculate_correlation(assay1, assay2, method, feature_names1, feature_names2, p_adj_method)
        correlations <- results$correlations
        p_values <- results$p_values
        p_values_adjusted <- results$p_values_adjusted
        
        result <- .correlation_filter(correlations, 
                                      p_values, 
                                      p_values_adjusted, 
                                      p_adj_threshold,
                                      n_signif, 
                                      cth,
                                      order)
        
        # Filter self correlations if it's specified
        if (nrow(assay1) == nrow(assay2) && ncol(assay1) == ncol(assay2) && filter_self_correlations) {
          diag(res$cor) <- diag(res$pval) <- diag(res$p.adj) <- NA
        }
        
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
             verbose = TRUE,
             mode = "table",
             p_adj_method = "fdr",
             p_adj_threshold = 0.05,
             n_signif = 0,
             cth = NULL,
             order = FALSE,
             filter_self_correlations = FALSE){
        exp1 <- x
        exp2 <- altExp(x, experiment2)
        x <- MultiAssayExperiment::MultiAssayExperiment(experiments = ExperimentList(exp1 = exp1, exp2 = exp2))
        getExperimentCrossCorrelation(x,
                                      experiment1 = 1,
                                      experiment2 = 2,
                                      abund_values1,
                                      abund_values2,
                                      method,
                                      verbose,
                                      mode,
                                      p_adj_method,
                                      p_adj_threshold,
                                      n_signif,
                                      cth,
                                      order,
                                      filter_self_correlations)
    }
)

######################## HELP FUNCTIONS #############################
.cor_test_data_type <- function(assay, method){
  # Different available methods
  numeric_methods <- c("spearman", "pearson")
  categorical_methods <- c("categorical")
  # Check if method match with values, othewise give warning.
  # For numeric methods, get only numeric values. For categorical methods, get only factors.
  if (method %in% numeric_methods) {
    inds <- vapply(assay, is.numeric, TRUE)
    if (any(!inds)) {
      warning("Considering only numeric annotations for \n       
                    pearson/spearman")
    }
    inds <- names(which(inds))
  } else if (method %in% categorical_methods) {
    inds <- vapply(assay, is.factor, TRUE)
    if (any(!inds)) {
      warning("Considering only categorical annotations for factors")
    }
    inds <- names(which(inds))
  }
  # Names of features
  names <- inds
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

.calculate_correlation <- function(assay1, assay2, method, feature_names1, feature_names2, p_adj_method){
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
      ###############################################################
      
      # js object includes all the correlation between individual feature from assay2
      # and all the featurs from assay1
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
        
        keep <- rowSums(is.na(cbind(xvec, yvec))) == 0
        xvec <- xvec[keep]
        yvec <- yvec[keep]
        
        # Number of data-annotation samples for
        # calculating the correlations
        n <- sum(keep)
        correlations[varname, lev] <- gktau(xvec, yvec) 
        
      }
    }
  }
  
  if (!all(is.na(p_values))) {
    
    rownames(p_values) <- feature_names1
    colnames(p_values) <- feature_names2
    
    rownames(correlations) <- feature_names1
    colnames(p_values) <- feature_names2
    
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
                                cth, order){
  
  # Filter
  if (!is.null(p_adj_threshold) || !is.null(cth)) {
    
    # Replace NAs with extreme values for filtering purposes
    p_values_adjusted[is.na(p_values_adjusted)] <- 1
    p_values[is.na(p_values_adjusted)] <- 1
    correlations[is.na(correlations)] <- 0
    
    # Filter by adjusted pvalues and correlations
    inds1.q <- inds2.q <- inds1.c <- inds2.c <- NULL
    
    if (!is.null(p_adj_threshold)) {
      
      inds1.q <- apply(p_values_adjusted, 1, function(x) {
        sum(x < p_adj_threshold) >= n_signif
      })
      
      inds2.q <- apply(p_values_adjusted, 2, function(x) {
        sum(x < p_adj_threshold) >= n_signif
      })
    }
    
    if (!is.null(cth)) {
      inds1.c <- apply(abs(correlations), 1, function(x) {
        sum(x > cth) >= n_signif
      })
      inds2.c <- apply(abs(correlations), 2, function(x) {
        sum(x > cth) >= n_signif
      })
    }
    
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
    
    Cmat <- as.matrix(0)
    
    # TODO: add also correlation filter, not only significance
    # Require each has at least n_signif. correlations
    
    if (sum(inds1) >= n_signif && sum(inds2) >= n_signif) {
      
      rnams <- rownames(correlations)[inds1]
      cnams <- colnames(correlations)[inds2]
      
      correlations <- matrix(correlations[inds1, inds2, drop=FALSE], nrow=sum(inds1))
      p_values <- matrix(p_values[inds1, inds2, drop=FALSE], nrow=sum(inds1))
      p_values_adjusted <- matrix(p_values_adjusted[inds1, inds2, drop=FALSE], nrow=sum(inds1))
      
      rownames(p_values_adjusted) <- rownames(p_values) <- rownames(correlations) <- rnams
      colnames(p_values_adjusted) <- colnames(p_values) <- colnames(correlations) <- cnams
      
      if (order && sum(inds1) >= 2 && sum(inds2) >= 2) {
        
        # Order in visually appealing order
        tmp <- correlations
        rownames(tmp) <- NULL
        colnames(tmp) <- NULL
        
        rind <- hclust(as.dist(1 - cor(t(tmp),
                                       use="pairwise.complete.obs")))$order
        cind <- hclust(as.dist(1 - cor(tmp,
                                       use="pairwise.complete.obs")))$order
        
        rnams <- rownames(correlations)[rind]
        cnams <- colnames(correlations)[cind]
        
        correlations <- correlations[rind, cind]
        p_values <- p_values[rind, cind]
        p_values_adjusted <- p_values_adjusted[rind, cind]
        
        rownames(p_values_adjusted) <- rownames(p_values) <- rownames(correlations) <- rnams
        colnames(p_values_adjusted) <- colnames(p_values) <- colnames(correlations) <- cnams
        
      }
      
    } else {
      message("No significant correlations with the given criteria\n")
      correlations <- p_values <- p_values_adjusted <- NULL
    }
  }
  ################################################################################
  
  res <- list(cor=correlations, pval=p_values, p.adj=p_values_adjusted)
  
  return(res)
  
}





.correlation_matrix_to_table <- function(res, verbose=FALSE) {
  
  ctab <- ID <- NULL
  
  if (!is.null(res$cor)) {
    ctab <- as.data.frame(res$cor)
    ctab$ID <- rownames(res$cor)
    ctab <- reshape2::melt(ctab, "ID")
    
    colnames(ctab) <- c("X1", "X2", "Correlation")
    ctab$Correlation <- as.numeric(as.character(ctab$Correlation))
  }
  
  correlation <- NULL  # circumwent warning on globabl vars
  
  if (!is.null(res$p.adj)) {
    
    if (verbose) {
      message("Arranging the table")
    }
    
    ctab2 <- as.data.frame(res$p.adj)
    ctab2$ID <- rownames(res$p.adj)
    ctab2 <- reshape2::melt(ctab2, "ID")
    colnames(ctab2) <- c("X1", "X2", "p.adj")
    ctab2$p.adj <- as.numeric(as.character(ctab2$p.adj))
    
    ctab <- cbind(ctab, ctab2$p.adj)
    colnames(ctab) <- c("X1", "X2", "Correlation", "p.adj")
    ctab <- ctab[order(ctab$p.adj), ]
    colnames(ctab) <- c("X1", "X2", "Correlation", "p.adj")
    
  } else {
    message("No significant adjusted p-values")
    if (!is.null(ctab)) {
      
      ctab2 <- as.data.frame(res$pval)
      ctab2$ID <- rownames(res$pval)
      ctab2 <- reshape2::melt(ctab2, "ID")
      colnames(ctab2) <- c("X1", "X2", "value")
      ctab2$value <- as.numeric(as.character(ctab2$value))
      
      ctab <- cbind(ctab, ctab2$value)
      ctab <- ctab[order(-abs(ctab$Correlation)), ]
      colnames(ctab) <- c("X1", "X2", "Correlation", "pvalue")
    }
  }
  
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
  
  ctab
  
}

