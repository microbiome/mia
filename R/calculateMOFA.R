#' calculateMOFA
#'
#' \code{calculateMOFA} produces a prepared MOFA2 model from a MAE.
#'
#' @param mae a \code{\link[MultiAssayExperiment:MultiAssayExperiment]{MultiAssayExperiment}}
#'
#' @name calculateMOFA
NULL

#' @rdname calculateMOFA
#' @export
setGeneric("calculateMOFA", signature = c("mae"),
           function(mae, ...) standardGeneric("calculateMOFA")
)

#' @rdname calculateMOFA
#' @export
setMethod("calculateMOFA", signature = c(mae = "MultiAssayExperiment"),
    function(mae, assay.types = rep("counts", length(mae)),
             groups = NULL, extract_metadata = FALSE, ...){
        # Select assays of each experiment for MOFA
        mae <- .select_assays(mae, assay.types)
        # Create MOFA from selected experiments
        mdl <- create_mofa_from_MultiAssayExperiment(
            mae,
            groups = groups,
            extract_metadata = extract_metadata
        )
        # Make a list of arguments for MOFA
        mofa_args <- list(
            object = mdl,
            data_options = .set_opts(get_default_data_options(mdl), ...),
            model_options = .set_opts(get_default_model_options(mdl), ...),
            training_options = .set_opts(get_default_training_options(mdl), ...),
            mefisto_options = .set_opts(get_default_mefisto_options(mdl), ...)
        )
        # Add stochastic options if stochastic is turned on
        if ( mofa_args[["training_options"]][["stochastic"]] ){
            mofa_args[["stochastic_options"]] <- .set_opts(get_default_stochastic_options(mdl), ...)
        }
        # Prepare MOFA
        prep_mdl <- do.call("prepare_mofa", mofa_args)
        return(prep_mdl)
    }
)

########################## HELP FUNCTIONS ##########################

# Select assays to be included in MOFA
#' @importFrom MultiAssayExperiment experiments
#' @importFrom SummarizedExperiment assay assays assayNames
.select_assays <- function(mae, assay.types) {
  # Give corresponding experiment names to assay.types
  names(assay.types) <- names(experiments(mae))
  # For every experiment in MAE
  for ( exp in names(experiments(mae)) ){
    # Keep only selected assay.type from a given experiment
    assays(mae[[exp]]) <- list(assay(mae[[exp]], assay.types[[exp]]))
    # Update assay names
    assayNames(mae[[exp]]) <- assay.types[[exp]]
  }
  return(mae)
}

# Combine custom options found in ... with default options
.set_opts <- function(default, ...) {
  # For every option in a set (data, model, train, ...)
  for ( opt in names(default) ){
    # If that option is found among arguments
    if ( opt %in% names(list(...)) ){
      # Replace default with value specified in arguments
      default[[opt]] <- list(...)[[opt]]
    }
  }
  return(default)
}
