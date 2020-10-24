#' @title meltAssay
#'
#' @description Converts a \code{\link[=TreeSummarizedExperiment-class]{TreeSummarizedExperiment}} object into a
#' long data.frame which can be used for \code{\link[tidyverse]{tidyverse}}-tools.
#'
#' @param x a numeric matrix or a \code{\link[=TreeSummarizedExperiment-class]{TreeSummarizedExperiment}}
#'          object containing a tree.
#' @param abund_values Must be one of the values returned by assayNames.
#' @param add_col_data Choice of colData columns to add. Default="none".
#' @param ... optional arguments not used.
#' @return A data.frame tbl_df.
#' @examples
#' \dontrun{
#' data(GlobalPatterns)
#' data(GlobalPatterns, package = "MicrobiomeExperiment")
#' se <- GlobalPatterns
#' se.rel <- relAbundanceCounts(se)
#' molten_se <- meltAssay(se.rel,
#'   add_col_data = colnames(colData(se.rel)),
#'   abund_values = "counts"
#' )
#' head(molten_se)
#' }
#' @importFrom SummarizedExperiment assayNames rowData rowData<-
#' @importFrom tibble rownames_to_column
#' @importFrom tidyr pivot_longer
#' @importFrom dplyr left_join mutate_if %>%
#' @keywords Utilities
#' @export
setGeneric("meltAssay",
  signature = "x",
  function(x, add_col_data = "none", abund_values = "counts", ...) {
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
  signature = c(x = "TreeSummarizedExperiment"),
  function(x, add_col_data = "none", abund_values = "counts") {

    # input check
    .check_abund_values(abund_values, x)

    uTaxaID <- uSamId <- Abundance <- NULL

    dfSE <- .melt_assay(x, abund_values)
    dfSE2 <- .add_taxonomic_data_to_molten_assay(x, dfSE)
    dfSE3 <- .add_col_data_to_molten_assay(dfSE2, x, add_col_data)
    return(dfSE3)
  }
)



#' @description Melts assay in \code{\link[=TreeSummarizedExperiment-class]{TreeSummarizedExperiment}}
#'        object to long data format.
#' @param x a numeric matrix or a \code{\link[=TreeSummarizedExperiment-class]{TreeSummarizedExperiment}}
#'          object containing a tree.
#' @param abund_values Must be one of the values returned by assayNames.
#'
.melt_assay <- function(x, abund_values) {
  uTaxaID <- NULL

  # input check
  .check_abund_values(abund_values, x)

  molten_assay <- assays(x)[[abund_values]] %>%
    data.frame() %>%
    rownames_to_column("uTaxaID") %>%
    # uSamId is unique sample id
    pivot_longer(!uTaxaID,
      values_to = "Abundance",
      names_to = "uSamId"
    )
  return(molten_assay)
}


#' @description Combines molten assay with rowData i.e. taxonomy table.
#' @param dfSE Molten SE.
#' @param x a \code{\link[=TreeSummarizedExperiment-class]{TreeSummarizedExperiment}}
#'          object.
#'
.add_taxonomic_data_to_molten_assay <- function(x, dfSE) {
  dfSE2 <- uTaxaID <- NULL
  if (is.null(rowData(x))) {
    message("Missing taxonomic information")
    # message("Missing taxonomic information \n returning only molten assay")

    dfSE2 <- dfSE
    return(dfSE2)
  } else {

    # Get taxonomic data stored in rowData
    me_tax <- rowData(x) %>%
      data.frame() %>%
      # uTaxaID is unique TaxaID
      rownames_to_column("uTaxaID")
    # may need to enforce in case of user supplied duplicate

    dfSE2 <- dfSE %>%
      left_join(me_tax, by = "uTaxaID") %>%
      mutate_if(is.factor, as.character)

    return(dfSE2)
  }
}

#' @description Combines molten assay and rowData i.e. taxonomy table with
#'              ColData i.e. sample information.
#' @param dfSE2 Molten SE.
#' @param x a \code{\link[=TreeSummarizedExperiment-class]{TreeSummarizedExperiment}}
#'          object.
#' @param add_col_data Choice of colData columns to add. Default="none".
#'
.add_col_data_to_molten_assay <- function(dfSE2, x, add_col_data) {
  dfSE3 <- uSamId <- NULL

  #if(add_col_data[1] == "none") {
   # add_col_data <- NULL
  #}

  # since it can be a vector consisting of more than one colnames(colData(x))
  # values, choose only 1
  if (is.null(add_col_data[1]) | is.na(add_col_data[1]) | add_col_data[1] == "none") {
    message("Returning molten assay without colData")

    dfSE3 <- dfSE2

  } else if (isTRUE(all(add_col_data %in% colnames(colData(x))))) {
    sub.colData <- colData(x)[, add_col_data]

    me_sam <- sub.colData %>%
      data.frame() %>%
      rownames_to_column("uSamId")

    # dfSE3 <- .add_taxonomic_data_to_molten_assay(x,dfSE2)

    dfSE3 <- dfSE2 %>%
      left_join(me_sam, by = "uSamId") %>%
      mutate_if(is.factor, as.character)
  }

  return(dfSE3)
}
