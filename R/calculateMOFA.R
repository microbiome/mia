# Transforming microbiome data with rclr
mae[[1]] <- transformAssay(mae[[1]], method = "relabundance")
mae[[1]] <- transformAssay(mae[[1]], assay.type = "relabundance", method = "rclr")

# Transforming metabolomic data with log10
mae[[2]] <- transformAssay(mae[[2]], assay.type = "nmr",
                           MARGIN = "samples",
                           method = "log10")

# Transforming biomarker data with z-transform
mae[[3]] <- transformAssay(mae[[3]], assay.type = "signals",
                           MARGIN = "features",
                           method = "z", pseudocount = 1)

# Building our mofa model

calculateMOFA <- function(mae, assay.types = rep("counts", length(mae)),
                          groups = NULL, extract_metadata = FALSE, ...) {
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

.set_opts <- function(default, ...) {
    for ( opt in names(default) ){
        if ( opt %in% names(list(...)) ){
            default[[opt]] <- list(...)[[opt]]
        }
    }
    return(default)
}

.select_assays <- function(mae, assay.types) {
    names(assay.types) <- names(experiments(mae))
    for ( exp in names(experiments(mae)) ){
        assays(mae[[exp]]) <- list(assay(mae[[exp]], assay.types[[exp]]))
        assayNames(mae[[exp]]) <- assay.types[[exp]]
    }
    return(mae)
}

prep_model <- calculateMOFA(mae, groups = "Diet",
                            assay.types = c("rclr", "log10", "z"),
                            extract_metadata = TRUE,
                            num_factors = 5, seed = 50,
                            spikeslab_factors = TRUE, verbose = FALSE,
                            stochastic = TRUE)

# Some systems may require the specification `use_basilisk = TRUE`
# so it has been added to the following code
model.trained <- run_mofa(model.prepared, use_basilisk = TRUE)

mae <- .select_assays(mae, c("rclr", "log10", "z"))
