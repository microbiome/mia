#' Perform mediation analysis
#'
#' \code{getMediation} provides a wrapper of \code{\link[mediation:mediate]{mediate}}
#' for \code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}.
#'
#' @param se a \code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}.
#' 
#' @param outcome The name of the colData variable used as outcome in the model.
#' 
#' @param treatment The name of the colData variable used as treatment in the model.
#'
#' @param mediator The name of the colData variable used as mediator in the model.
#'   (default: \code{mediator = NULL})
#'
#' @param assay.type The name of the assay used for feature-wise mediation analysis.
#'   (default: \code{assay.type = NULL})
#' 
#' @param dim.type The name of the reducedDim used for component-wise mediation analysis.
#'   (default: \code{dim.type = NULL})
#'
#' @param family A specification for the outcome model link function.
#'    (default: \code{family = gaussian("identity")})
#' 
#' @param covariates Name or list of names of colData variables used as covariates
#'   in the model. (default: \code{covariates = NULL})
#' 
#' @param p.adj.method A single character value for selecting adjustment method
#'   of p-values. Passed to `p.adjust` function. (default: \code{p.adj.method = "BH"})
#' 
#' @param add.metadata TRUE or FALSE, should the model metadata be returned.
#'   (default: \code{add.metadata = FALSE})
#' 
#' @param message TRUE or FALSE, should execution messages be printed.
#'   (default: \code{message = TRUE})
#' 
#' @param ... additional parameters that can be passed to \code{\link[mediation:mediate]{mediate}}.
#' 
#' @details
#' This wrapper of \code{\link[mediation:mediate]{mediate}} for
#' \code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#' provides a simple method to analyse the effect of a treatment variable on an
#' outcome variable found in \code{colData(se)} through the mediation of either
#' another variable in colData (argument \code{mediator}) or an assay
#' (argument \code{assay.type}) or a reducedDim (argument \code{dim.type}). Notably,
#' those three arguments are mutually exclusive.
#' 
#' @return 
#' A \code{data.frame} of p-values, effect sizes and confidence intervals for
#' the ACMEs and ADEs of every mediator included in the analysis. Columns include:
#' 
#' \describe{
#'   \item{Treatment}{the treatment variable}
#'   \item{Mediator}{the mediator variable}
#'   \item{Outcome}{the outcome variable}
#'   \item{ACME_estimate}{the ACME estimate}
#'   \item{ADE_estimate}{the ADE estimate}
#'   \item{ACME_pval}{the p-value for the ACME estimate}
#'   \item{ADE_pval}{the p-value for the ADE estimate}
#'   \item{ACME_adjpval}{the adjusted p-value for the ACME estimate}
#'   \item{ADE_adjpval}{the adjusted p-value for the ADE estimate}
#'   \item{ACME_CI_lower}{the lower limit of the confidence interval for the ACME estimate}
#'   \item{ACME_CI_upper}{the upper limit of the confidence interval for the ACME estimate}
#'   \item{ADE_CI_lower}{the lower limit of the confidence interval for the ADE estimate}
#'   \item{ADE_CI_upper}{the upper limit of the confidence interval for the ADE estimate}
#' }
#'
#' @name getMediation
#'
#' @examples
#' # Import libraries
#' library(mia)
#' library(microbiomeDataSets)
#' 
#' # Load dataset
#' tse <- LahtiWAData()
#'  
#' # Convert BMI variable to numeric
#' tse$bmi_group <- as.numeric(tse$bmi_group)
#' 
#' # Analyse mediated effect of nationality on BMI through alpha diversity 
#' med_df <- getMediation(tse,
#'                        outcome = "bmi_group",
#'                        treatment = "nationality",
#'                        mediator = "diversity",
#'                        covariates = c("sex", "age"),
#'                        treat.value = "Scandinavia",
#'                        control.value = "CentralEurope",
#'                        boot = TRUE, sims = 1000)
#' 
#' # Apply clr transformation to counts assay
#' tse <- transformAssay(tse,
#'                       method = "clr",
#'                       pseudocount = 1)
#'
#' # Analyse mediated effect of nationality on BMI through clr-transformed features     
#' med_df <- getMediation(tse,
#'                        outcome = "bmi_group",
#'                        treatment = "nationality",
#'                        assay.type = "clr",
#'                        covariates = c("sex", "age"),
#'                        treat.value = "Scandinavia",
#'                        control.value = "CentralEurope",
#'                        boot = TRUE, sims = 300)
#'  
#' # Perform ordination
#' tse <- runNMDS(tse,
#'                method = "euclidean",
#'                assay.type = "clr",
#'                name = "NMDS")
#' 
#' # Analyse mediated effect of nationality on BMI through NMDS components        
#' med_df <- getMediation(tse,
#'                        outcome = "bmi_group",
#'                        treatment = "nationality",
#'                        dim.type = "NMDS",
#'                        covariates = c("sex", "age"),
#'                        treat.value = "Scandinavia",
#'                        control.value = "CentralEurope",
#'                        boot = TRUE, sims = 500)
#' 
NULL

#' @rdname getMediation
#' @export
setGeneric("getMediation", signature = c("se"),
           function(se, ...) standardGeneric("getMediation"))

#' @rdname getMediation
#' @export
setMethod("getMediation", signature = c(se = "SummarizedExperiment"),
      function(se, outcome, treatment,
               mediator = NULL, assay.type = NULL, dim.type = NULL,
               family = gaussian(), covariates = NULL, p.adj.method = "BH",
               add.metadata = FALSE, message = TRUE, ...) {
            
        results <- list(Treatment = c(), Mediator = c(), Outcome = c(),
                        ACME_estimate = c(), ADE_estimate = c(), ACME_pval = c(),
                        ADE_pval = c(), ACME_ci = c(), ADE_ci = c(), Model = list())
            
        med_opts <- sapply(
          list(mediator = mediator, assay.type = assay.type, dim.type = dim.type),
          function(x) !is.null(x)
        )
            
        if (sum(med_opts) == 1) {
          if (med_opts[[1]]) {
                
            if (!mediator %in% names(colData(se))) {
              stop(mediator, " not found in colData(se).", call. = FALSE)
            }
                
            mediators <- mediator
            mat <- NULL
                
          } else if (med_opts[[2]]) {
                
            if (!assay.type %in% assayNames(se)) {
              stop(assay.type, " not found in assays(se).", call. = FALSE)
            }
                
            mat <- assay(se, assay.type)
            mediators <- rownames(mat)
                
          } else if (med_opts[[3]]) {
                
            if (!dim.type %in% reducedDimNames(se)) {
              stop(dim.type, " not found in reducedDims(se).", call. = FALSE)
            }
                
            mat <- t(reducedDim(se, dim.type))
            rownames(mat) <- paste0(dim.type, seq(1, nrow(mat)))
            mediators <- rownames(mat)
                
          }
              
        } else {
          stop("The arguments mediator, assay.type and dim.type are mutually exclusive",
               ", but ", sum(med_opts), " were provided.", call. = FALSE)
        }
            
        i <- 0
            
        for (mediator in mediators) {
              
          i <- i + 1
              
          if (message) {
            print(paste0("Mediator ", i, " out of ", length(mediators), ": ", mediator)) 
          }
              
          med_out <- wrap.mediate(se, outcome, treatment, mediator,
                                  family = family, mat = mat,
                                  covariates = covariates, ...)
              
          results <- update.results(results, med_out, treatment, mediator, outcome)
              
        }
            
        med_df <- make.output(results, p.adj.method, add.metadata)
        return(med_df)
    }
)



check.args <- function(df, ...) {
  
  kwargs <- list(...)
  
  if (any(is.na(df))) {
    
    total <- nrow(df)
    df <- na.omit(df)
    message(paste(total - nrow(df), "samples removed because of missing data."))
    
  }
  
  if (!is.numeric(df[["Treatment"]]) & length(unique((df[["Treatment"]]))) > 2) {
    
    multilevel_message <- paste(
      "Too many treatment levels. Consider specifing a treat.value and a control.value\n"
    )
    
    if (!is.null(kwargs[["boot"]])) {
      
      if (any(sapply(kwargs[c("control.value", "treat.value")], is.null))) {
        stop(multilevel_message, call. = FALSE)
      } else {
        
        if (!all(kwargs[c("control.value", "treat.value")] %in% unique(df[["Treatment"]]))) {
          stop(multilevel_message, call. = FALSE)
        }
        
        keep <- df[["Treatment"]] %in% kwargs[c("control.value", "treat.value")]
        message(paste(nrow(df) - sum(keep), "samples removed because different",
                      "from control and treatment."))
        
        df <- df[keep, ]
        
      }
    }
  }
  return(df)
}



# utility function to run mediation
run.mediation <- function(df, family,
                          covariates = NULL,
                          relation_m, relation_dv,
                          ...) {
  
  df <- check.args(df, ...)
  
  fit_m <- do.call(lm, list(formula = formula(relation_m),
                            data = df))
  
  fit_dv <- do.call(glm, list(formula = formula(relation_dv),
                              family = family,
                              data = df))
  
  med_out <- mediate(fit_m, fit_dv,
                     treat = "Treatment",
                     mediator = "Mediator",
                     covariates = covariates,
                     ...)
  
  return(med_out)
  
}



wrap.mediate <- function(se, outcome, treatment, mediator,
                         family = gaussian(), mat = NULL,
                         covariates = NULL, ...) {
  
  df <- data.frame(Outcome = eval(parse(text = paste0("se$", outcome))),
                   Treatment = eval(parse(text = paste0("se$", treatment))))
  
  if (is.null(mat)) {
    df[["Mediator"]] <- eval(parse(text = paste0("se$", mediator)))
  } else {
    df[["Mediator"]] <- mat[mediator, ]
  }
  
  relation_m <- "Mediator ~ Treatment"
  relation_dv <- "Outcome ~ Treatment + Mediator"
  
  if (!is.null(covariates)) {
    for (covariate in covariates) {
      
      df[[covariate]] <- eval(parse(text = paste0("se$", covariate)))
      
      relation_m <- paste(relation_m, "+", covariate)
      relation_dv <- paste(relation_dv, "+", covariate)
      
    }
  }
  
  med_out <- run.mediation(df, family,
                           relation_m, relation_dv,
                           covariates = covariates, ...)
  
  return(med_out)
}



update.results <- function(results, med_out,
                           treatment, mediator, outcome) {
  
  # variables
  results[["Treatment"]] <- c(results[["Treatment"]], treatment)
  results[["Mediator"]] <- c(results[["Mediator"]], mediator)
  results[["Outcome"]] <- c(results[["Outcome"]], outcome)
  
  # ACME (average causal mediation effect)
  results[["ACME_estimate"]] <- c(results[["ACME_estimate"]], med_out$d.avg)
  results[["ACME_pval"]] <- c(results[["ACME_pval"]], med_out$d.avg.p)
  results[["ACME_ci"]] <- c(results[["ACME_ci"]], med_out$d.avg.ci)
  
  # ADE (average direct effect)
  results[["ADE_estimate"]] <- c(results[["ADE_estimate"]], med_out$z.avg)
  results[["ADE_pval"]] <- c(results[["ADE_pval"]], med_out$z.avg.p)
  results[["ADE_ci"]] <- c(results[["ADE_ci"]], med_out$z.avg.ci)
  
  # Raw model as metadata
  results[["Model"]][[length(results[["Model"]]) + 1]] <- med_out
  
  return(results)
}


make.output <- function(results, p.adj.method, add.metadata) {
  
  med_df <- do.call(data.frame, results[1:(length(results) - 3)])
  
  med_df[["ACME_adjpval"]] <- p.adjust(med_df[["ACME_pval"]], method = p.adj.method)
  med_df[["ADE_adjpval"]] <- p.adjust(med_df[["ADE_pval"]], method = p.adj.method)
  
  med_df[["ACME_CI_lower"]] <- results[["ACME_ci"]][names(results[["ACME_ci"]]) == "2.5%"]
  med_df[["ACME_CI_upper"]] <- results[["ACME_ci"]][names(results[["ACME_ci"]]) == "97.5%"]
  med_df[["ADE_CI_lower"]] <- results[["ADE_ci"]][names(results[["ADE_ci"]]) == "2.5%"]
  med_df[["ADE_CI_upper"]] <- results[["ADE_ci"]][names(results[["ADE_ci"]]) == "97.5%"]
  
  if (add.metadata) {
    attr(med_df, "metadata") <- results[["Model"]]
  }
  
  med_df <- med_df[order(med_df[["ACME_pval"]]), ]
  
  return(med_df)
}