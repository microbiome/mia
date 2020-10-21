#' @title meltAssay
#'
#' Converts a \code{\link[=MicrobiomeExperiment-class]{MicrobiomeExperiment}} object into a
#' long data.frame which can be used for \code{\link[tidyverse]{tidyverse}}-tools.
#'
#' @param x a numeric matrix or a \code{\link[=MicrobiomeExperiment-class]{MicrobiomeExperiment}}
#'          object containing a tree.
#' @param abund_values Default is counts, must be one of "counts", "relabundance", "both".
#' @param colData_choice Choice of colData columns to add. By default it will consider `"all"`
#' columns. Specific columns can be provided as `c("col1", "col2", ...)`. If `colData_choice=NA`
#' only the specified `assayName` will be converted to long data.frame
#' @param ... optional arguments not used.
#' @return A data.frame.
#' @importFrom SummarizedExperiment assayNames rowData rowData<-
#' @importFrom tibble rownames_to_column
#' @importFrom tidyr pivot_longer
#' @importFrom dplyr left_join mutate_if %>%
#' @export
setGeneric("meltAssay",
  signature = "x",
  function(x, colData_choice = "all", abund_values = "counts", ...) {
    standardGeneric("meltAssay")
  }
)


#' @rdname meltAssay
#' @aliases meltAssay
#'
#' @importFrom SummarizedExperiment rowData rowData<- colData colData<-
#'
#' @export
setMethod("meltAssay",
  signature = c(x = "MicrobiomeExperiment"),
  function(x, colData_choice = "all", abund_values = "counts") {
    # input check
    if (isFALSE(any(abund_values == c("counts", "relabundance", "both")))) {
      stop(paste0(
        "'abund_values' must be one:  ",
        list(assayNames(x))
      ), " or", " 'both' ")
    }

    uTaxaID <- uSamId <- Abundance <- NULL

    if (is.null(rowData(x))) {
      warning("'rowData(x)' is NULL!")
    } else {
      # Get taxonomic data stored in rowData
      me_tax <- rowData(x) %>%
        data.frame() %>%
        # uTaxaID is unique TaxaID
        rownames_to_column("uTaxaID")
      # may need to enforce in case of user supplied duplicate
    }

    if (is.null(assay(x))) {
      stop("'rowData(x)' is NULL!")
    }

    if (abund_values == "counts") {
      dfSE <- assays(x)$counts %>%
        data.frame() %>%
        rownames_to_column("uTaxaID") %>%
        # uSamId is unique sample id
        pivot_longer(!uTaxaID,
          values_to = "Abundance",
          names_to = "uSamId"
        )
    } else if (abund_values == "relabundance") {
      dfSE <- assays(x)$relabundance %>%
        data.frame() %>%
        rownames_to_column("uTaxaID") %>%
        pivot_longer(!uTaxaID,
          values_to = "Abundance",
          names_to = "uSamId"
        )
    } else if (abund_values == "both") {
      me_count_ab <- assays(x)$counts %>%
        data.frame() %>%
        rownames_to_column("uTaxaID") %>%
        pivot_longer(!uTaxaID,
          values_to = "CountAbundance",
          names_to = "uSamId"
        )

      me_rel_ab <- assays(x)$relabundance %>%
        data.frame() %>%
        rownames_to_column("uTaxaID") %>%
        pivot_longer(!uTaxaID,
          values_to = "RelAbundance",
          names_to = "uSamId"
        )

      dfSE <- me_count_ab %>%
        left_join(me_rel_ab, by = c("uTaxaID", "uSamId"))
    }

    if (is.null(rowData(x)) & is.null(colData_choice) | is.na(colData_choice)) {
      return(dfSE)
    } else if (!is.null(rowData(x)) & is.null(colData_choice) | is.na(colData_choice)) {
      dfSE <- dfSE %>%
        left_join(me_tax, by = "uTaxaID") %>%
        mutate_if(is.factor, as.character)

      return(dfSE)
    }

    # check coldata
    if (colData_choice == "all") {
      me_sam <- colData(x) %>%
        data.frame() %>%
        rownames_to_column("uSamId")

      dfSE <- dfSE %>%
        left_join(me_sam, by = "uSamId") %>%
        mutate_if(is.factor, as.character)

      return(dfSE)
    } else if (isTRUE(any(colData_choice == colnames(colData(x))))) {
      colData_choice.n <- intersect(colData_choice, colnames(colData(x)))
      sub.colData <- colData(x)[, colData_choice.n]

      me_sam <- sub.colData %>%
        data.frame() %>%
        rownames_to_column("uSamId")

      dfSE <- dfSE %>%
        left_join(me_sam, by = "uSamId") %>%
        mutate_if(is.factor, as.character)

      return(dfSE)
    }
  }
)
