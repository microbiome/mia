#' Universal Joint RPCA Wrapper
#'
#' @param data Input object: MultiAssayExperiment, TreeSummarizedExperiment, SummarizedExperiment, list of matrices, or single matrix.
#' @param assays Character vector of assay names to extract (only used for MultiAssayExperiment). If NULL, use all assays.
#' @param ... Additional arguments passed to jointRPCA()
#'
#' @return Output from jointRPCA()
#' @export
jointRPCAuniversal <- function(data, assays = NULL, ...) {
    if (inherits(data, "MultiAssayExperiment")) {
        if (is.null(assays)) assays <- names(experiments(data))
        tables <- lapply(assays, function(a) assay(data, a))
        names(tables) <- assays  
    } else if (inherits(data, "TreeSummarizedExperiment") || inherits(data, "SummarizedExperiment")) {
        tables <- list(assay(data))
        nm <- tryCatch(assayNames(data)[1], error = function(e) NULL)
        names(tables) <- if (is.null(nm) || is.na(nm)) "assay1" else nm  
    } else if (is.list(data) && all(sapply(data, is.matrix))) {
        tables <- data
        if (is.null(names(tables))) names(tables) <- paste0("view", seq_along(tables))  
    } else if (is.matrix(data)) {
        tables <- list(data); names(tables) <- "assay1"  
    } else {
        stop("Unsupported input type ...")
    }
  
    jointRPCA(tables = tables, ...)
}