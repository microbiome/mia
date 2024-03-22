#' Perform mediation analysis
#'
#' \code{getMediation} provides a wrapper of \code{\link[mediation:mediate]{mediate}}
#' for \code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}.
#'
#' @param x a \code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}.
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
#' @param verbose TRUE or FALSE, should execution messages be printed.
#'   (default: \code{verbose = TRUE})
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
#'   \item{ACME_estimate}{the Average Causal Mediation Effect (ACME) estimate}
#'   \item{ADE_estimate}{the Average Direct Effect (ADE) estimate}
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
#' library(scater)
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
#' tse <- runMDS(tse, name = "MDS",
#'               method = "euclidean",
#'               assay.type = "clr",
#'               ncomponents = 3)
#' 
#' # Analyse mediated effect of nationality on BMI through NMDS components        
#' med_df <- getMediation(tse,
#'                        outcome = "bmi_group",
#'                        treatment = "nationality",
#'                        dim.type = "MDS",
#'                        covariates = c("sex", "age"),
#'                        treat.value = "Scandinavia",
#'                        control.value = "CentralEurope",
#'                        boot = TRUE, sims = 500)
#' 
NULL

#' @rdname getMediation
#' @export
setGeneric("getMediation", signature = c("x"),
           function(x, ...) standardGeneric("getMediation"))

#' @rdname getMediation
#' @export
setMethod("getMediation", signature = c(x = "SummarizedExperiment"),
        function(x, outcome, treatment,
                 mediator = NULL, assay.type = NULL, dim.type = NULL,
                 family = gaussian(), covariates = NULL, p.adj.method = "BH",
                 add.metadata = FALSE, verbose = TRUE, ...) {

        ###################### Input check ########################
        if( !outcome %in% names(colData(x)) ){
            stop(outcome, " not found in colData(x).", call. = FALSE)
        }
        if( !treatment %in% names(colData(x)) ){
            stop(treatment, " not found in colData(x).", call. = FALSE)
        }
        if( !is.null(covariates) && !all(covariates %in% names(colData(x))) ){
            stop("covariates not found in colData(x).", call. = FALSE)
        }
        if( !.is_a_bool(add.metadata) ){
            stop("add.metadata must be TRUE or FALSE.", call. = FALSE)
        }
        if( !.is_a_bool(verbose) ){
            stop("verbose must be TRUE or FALSE.", call. = FALSE)
        }
        
        # Check that arguments can be passed to mediate and remove unused samples
        x <- .check.mediate.args(x, outcome, treatment, mediator, covariates, verbose, ...)
        
        # Check which mediator was provided (colData, assay or reducedDim)
        med_opts <- sapply(
            list(mediator = mediator, assay.type = assay.type, dim.type = dim.type),
            function(x) !is.null(x)
        )
        
        if( sum(med_opts) != 1 ){
            # Throw error if none or multiple mediation options are specified
            stop(
                "The arguments mediator, assay.type and dim.type are mutually exclusive",
                ", but ", sum(med_opts), " were provided.", call. = FALSE
            )
        
        } else {
            if ( med_opts[[1]] ){
                # Check that mediator is in colData
                if( !mediator %in% names(colData(x)) ) {
                    stop(mediator, " not found in colData(x).", call. = FALSE)
                }
                # Use mediator for analysis  
                mediators <- mediator
                mat <- NULL
                
            } else if( med_opts[[2]] ){
                # Check that assay is in assays
                if( !assay.type %in% assayNames(x) ){
                    stop(assay.type, " not found in assays(x).", call. = FALSE)
                }
                # Define matrix for analysis
                mat <- assay(x, assay.type)
                # Use assay for analysis
                mediators <- rownames(mat)
                
            } else if( med_opts[[3]] ){
                # Check that reducedDim is in reducedDims
                if(!dim.type %in% reducedDimNames(x)){
                    stop(dim.type, " not found in reducedDims(x).", call. = FALSE)
                }
                # Define matrix for analysis
                mat <- t(reducedDim(x, dim.type))
                # Give component names to matrix rows
                rownames(mat) <- paste0(dim.type, seq(1, nrow(mat)))
                # Use reducedDim for analysis
                mediators <- rownames(mat)
            }
        }
        
        # Create template list of results
        results <- list(
            Treatment = c(), Mediator = c(), Outcome = c(),
            ACME_estimate = c(), ADE_estimate = c(), ACME_pval = c(),
            ADE_pval = c(), ACME_ci = c(), ADE_ci = c(), Model = list()
        )
        
        # Set initial index  
        i <- 0
            
        for( mediator in mediators ){
             
            # Update index 
            i <- i + 1
          
            if( verbose ){
                message("\rMediator ", i, " out of ", length(mediators), ": ", mediator, appendLF = FALSE) 
            }
          
            # Run mediation analysis for current mediator
            med_out <- .run.mediate(
                x, outcome, treatment, mediator,
                family = family, mat = mat,
                covariates = covariates, ...
            )
            
            # Update list of results
            results <- .update.results(results, med_out, treatment, mediator, outcome)
        }
        
        # Combine results into dataframe
        med_df <- .make.output(results, p.adj.method, add.metadata)
        return(med_df)
    }
)


# Check that arguments can be passed to mediate and remove unused samples
.check.mediate.args <- function(x, outcome, treatment, mediator, covariates, verbose = TRUE, ...) {
  
    # Create dataframe from selected columns of colData
    df <- as.data.frame(colData(x)[ , names(colData(x)) %in% c(outcome, treatment, mediator, covariates)])
    # Store kwargs into variable
    kwargs <- list(...)
  
    # Remove missing data from df
    df <- na.omit(df)
    diff <- ncol(x) - nrow(df)
  
    if( diff != 0 ){
        # Remove missing data from se
        x <- x[ , rownames(df)]
    
        if( verbose ){
            message(diff, " samples removed because of missing data.")
        }
    }
  
    # If the treatment variable has three or more levels
    if( !is.numeric(df[[treatment]]) && length(unique((df[[treatment]]))) > 2 ) {
        # and boot is TRUE
        if( !is.null(kwargs[["boot"]]) ){
            # and control and treat value are not specified
            if( any(sapply(kwargs[c("control.value", "treat.value")], is.null)) ){
                stop("Too many treatment levels. Consider specifing a treat.value and a control.value", call. = FALSE)
            # but if they are specified
            } else {
                # and they also appear in the treatment variable
                if( !all(kwargs[c("control.value", "treat.value")] %in% unique(df[[treatment]])) ){
                    stop("treat.value and/or control.value not found in the levels of the treatment variable.", call. = FALSE)
                }
        
                # Find indices of samples that belong to either control or treatment
                keep <- df[[treatment]] %in% kwargs[c("control.value", "treat.value")]
        
                # Remove samples different from control and treatment from df
                df <- df[keep, ]
                diff <- ncol(x) - nrow(df)
        
                # Remove samples different from control and treatment from se
                x <- x[ , rownames(df)]
        
                if( verbose ){
                    message(diff, " samples removed because different from control and treatment.")
                }
            }
        }
    }
    return(x)
}

# Run mediation analysis
#' @importFrom mediation mediate
.run.mediate <- function(x, outcome, treatment, mediator = NULL, mat = NULL,
                        family = gaussian(), covariates = NULL, ...) {
  
    # Create initial dataframe with outcome and treatment variables
    df <- data.frame(Outcome = colData(x)[[outcome]], Treatment = colData(x)[[treatment]])
  
    if( is.null(mat) ){
        # If matrix not given, fetch mediator from colData
        df[["Mediator"]] <- colData(x)[[mediator]]
    } else {
        # If matrix given, use it as mediators
        df[["Mediator"]] <- mat[mediator, ]
    }
  
    # Define basic formula mediation model
    relation_m <- "Mediator ~ Treatment"
    # Define basic formula outcome model
    relation_dv <- "Outcome ~ Treatment + Mediator"
  
    if( !is.null(covariates) ){
        for( covariate in covariates ){
        # Fetch covariate from colData and store it in dataframe
        df[[covariate]] <- colData(x)[[covariate]]
      
        # Add covariate to formula of mediation model
        relation_m <- paste(relation_m, "+", covariate)
        # Add covariate to formula of outcome model
        relation_dv <- paste(relation_dv, "+", covariate)
        }
    }

    # Fit mediation model
    fit_m <- do.call(lm, list(formula = formula(relation_m), data = df))
    # Fit outcome model
    fit_dv <- do.call(glm, list(formula = formula(relation_dv), family = family, data = df))
    # Run mediation analysis
    med_out <- mediate(fit_m, fit_dv, treat = "Treatment", mediator = "Mediator", covariates = covariates, ...)
  
    return(med_out)
}


# Update list of results
.update.results <- function(results, med_out, treatment, mediator, outcome) {
  
    # Update model variables
    results[["Treatment"]] <- c(results[["Treatment"]], treatment)
    results[["Mediator"]] <- c(results[["Mediator"]], mediator)
    results[["Outcome"]] <- c(results[["Outcome"]], outcome)
  
    # Update stats of ACME (average causal mediation effect)
    results[["ACME_estimate"]] <- c(results[["ACME_estimate"]], med_out$d.avg)
    results[["ACME_pval"]] <- c(results[["ACME_pval"]], med_out$d.avg.p)
    results[["ACME_ci"]] <- c(results[["ACME_ci"]], med_out$d.avg.ci)
  
    # Update stats of ADE (average direct effect)
    results[["ADE_estimate"]] <- c(results[["ADE_estimate"]], med_out$z.avg)
    results[["ADE_pval"]] <- c(results[["ADE_pval"]], med_out$z.avg.p)
    results[["ADE_ci"]] <- c(results[["ADE_ci"]], med_out$z.avg.ci)
  
    # Add current model to metadata
    results[["Model"]][[length(results[["Model"]]) + 1]] <- med_out
  
    return(results)
}


# Combine results into output dataframe
.make.output <- function(results, p.adj.method, add.metadata) {
  
    # Create dataframe with model variables, effect sizes and p-values
    med_df <- do.call(data.frame, results[1:(length(results) - 3)])
  
    # Compute adjusted p-values and add them to dataframe
    med_df[["ACME_adjpval"]] <- p.adjust(med_df[["ACME_pval"]], method = p.adj.method)
    med_df[["ADE_adjpval"]] <- p.adjust(med_df[["ADE_pval"]], method = p.adj.method)
  
    # Split CI lists into lower and upper limits and add them to dataframe
    med_df[["ACME_CI_lower"]] <- results[["ACME_ci"]][names(results[["ACME_ci"]]) == "2.5%"]
    med_df[["ACME_CI_upper"]] <- results[["ACME_ci"]][names(results[["ACME_ci"]]) == "97.5%"]
    med_df[["ADE_CI_lower"]] <- results[["ADE_ci"]][names(results[["ADE_ci"]]) == "2.5%"]
    med_df[["ADE_CI_upper"]] <- results[["ADE_ci"]][names(results[["ADE_ci"]]) == "97.5%"]
  
    if( add.metadata ){
        # If desired, the models for every mediator are saved into the metadata attribute
        attr(med_df, "metadata") <- results[["Model"]]
    }
  
    # Order output dataframe by ACME p-values
    med_df <- med_df[order(med_df[["ACME_pval"]]), ]
  
    return(med_df)
}