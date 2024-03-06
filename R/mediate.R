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


check.args <- function(df, ...) {
  
  kwargs <- list(...)
  
  if (any(is.na(df))) {
    
    df <- na.omit(df)
    message(paste(ncol(tse) - nrow(df), "samples removed because of missing data."))
    
  }
  
  if (!is.numeric(df$Treatment) & n_distinct(df$Treatment) > 2) {
    
    multilevel_message <- paste(
      "Too many treatment levels. Consider specifing a treat.value and a control.value\n"
    )
    
    if (!is.null(kwargs[["boot"]])) {
      
      if (any(sapply(kwargs[c("control.value", "treat.value")], is.null))) {
        stop(multilevel_message, call. = FALSE)
      } else {
        
        if (!all(kwargs[c("control.value", "treat.value")] %in% unique(df$Treatment))) {
          stop(multilevel_message, call. = FALSE)
        }
        
        keep <- df$Treatment %in% kwargs[c("control.value", "treat.value")]
        message(paste(nrow(df) - sum(keep), "samples removed because different",
                      "from control and treatment."))
        
        df <- df[keep, ]
        
      }
    }
  }
  return(df)
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
  
  return(results)
}


make.output <- function(results, p.adj.method) {
  
  med_df <- do.call(data.frame, results[1:(length(results) - 2)])
  
  med_df[["ACME_adjpval"]] <- p.adjust(med_df[["ACME_pval"]], method = p.adj.method)
  med_df[["ADE_adjpval"]] <- p.adjust(med_df[["ADE_pval"]], method = p.adj.method)
  
  med_df[["ACME_CI_lower"]] <- results[["ACME_ci"]][names(results[["ACME_ci"]]) == "2.5%"]
  med_df[["ACME_CI_lower"]] <- results[["ACME_ci"]][names(results[["ACME_ci"]]) == "97.5%"]
  med_df[["ADE_CI_lower"]] <- results[["ADE_ci"]][names(results[["ADE_ci"]]) == "2.5%"]
  med_df[["ADE_CI_lower"]] <- results[["ADE_ci"]][names(results[["ADE_ci"]]) == "97.5%"]
  
  med_df <- arrange(med_df, ACME_pval, ACME_estimate)
  
  return(med_df)
}


# main function to mediate coldata
mediateColData <- function(tse, outcome, treatment, mediator,
                           family = gaussian(), mat = NULL,
                           covariates = NULL, ...) {
  
  df <- data.frame(Outcome = eval(parse(text = paste0("tse$", outcome))),
                   Treatment = eval(parse(text = paste0("tse$", treatment))))
  
  if (is.null(mat)) {
    df[["Mediator"]] <- eval(parse(text = paste0("tse$", mediator)))
  } else {
    df[["Mediator"]] <- mat[mediator, ]
  }
  
  relation_m <- "Mediator ~ Treatment"
  relation_dv <- "Outcome ~ Treatment + Mediator"
  
  if (!is.null(covariates)) {
    for (covariate in covariates) {
      
      df[[covariate]] <- eval(parse(text = paste0("tse$", covariate)))
      
      relation_m <- paste(relation_m, "+", covariate)
      relation_dv <- paste(relation_dv, "+", covariate)
      
    }
  }
  
  med_out <- run.mediation(df, family,
                           relation_m, relation_dv,
                           covariates = covariates, ...)
  
  return(med_out)
}


# main function to mediate assay or reduced dimension
mediateAssay <- function(tse, outcome, treatment,
                         assay.type = NULL, dim.type = NULL,
                         family = gaussian(), covariates = NULL,
                         p.adj.method = "BH", ...) {
  
  results <- list(Treatment = c(), Mediator = c(), Outcome = c(),
                  ACME_estimate = c(), ADE_estimate = c(), ACME_pval = c(),
                  ADE_pval = c(), ACME_ci = c(), ADE_ci = c())
  
  if (!is.null(assay.type)) {
    
    mat <- assay(tse, assay.type)
    
  } else if (!is.null(dim.type)) {
    
    mat <- t(reducedDim(tse, dim.type))
    rownames(mat) <- paste0(dim.type, seq(1, nrow(mat)))
    
  }
  
  mediators <- rownames(mat)
  i <- 0
  
  for (mediator in mediators) {
    
    print(paste(length(mediators) - i, "left"))
    print(paste("Current mediator:", mediator))
    i <- i + 1
    
    med_out <- mediateColData(tse, outcome, treatment, mediator,
                              family = family, mat = mat,
                              covariates = covariates, ...)
    
    results <- update.results(results, med_out, treatment, mediator, outcome)
    
  }
  
  med_df <- make.output(results, p.adj.method)
  return(med_df)
}
