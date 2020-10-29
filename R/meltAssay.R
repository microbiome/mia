#' Converts a \code{\link[=SummarizedExperiment-class]{SummarizedExperiment}}
#'     object into a long data.frame
#'
#' Converts a \code{\link[=SummarizedExperiment-class]{SummarizedExperiment}}
#' object into a long data.frame which can be used for
#' \code{\link[tidyverse]{tidyverse}}-tools.
#'
#' @param x A numeric matrix or \code{\link[=SummarizedExperiment-class]{SummarizedExperiment}}
#'     object containing a tree.
#'
#' @param abund_values Must be one of the values returned by assayNames.
#'
#' @param add_col_data Choice of colData columns to add, can be a single
#'    column name in colData or multiple c("col1", "col2",...). Default=NULL.
#'
#' @param add_row_data Choice of rowData i.e. taxonomic levels to add, can be a single
#'    column name in colData or multiple c("col1", "col2",...). Default=NULL.
#'
#' @param ... optional arguments not used.
#'
#' @return A data.frame
#'
#' @name meltAssay
#'
#' @export
#'
#' @importFrom SummarizedExperiment assayNames rowData rowData<-
#'
#' @importFrom tibble rownames_to_column
#'
#' @importFrom tidyr pivot_longer
#'
#' @importFrom dplyr left_join rename mutate_if %>%
#'
#' @importFrom MicrobiomeExperiment taxonomyRanks
#'
#' @author
#' Sudarshan A. Shetty
#'
#' @examples
#' data(GlobalPatterns)
#' data(GlobalPatterns, package = "MicrobiomeExperiment")
#' se <- GlobalPatterns
#' se.rel <- relAbundanceCounts(se)
#' molten_se <- meltAssay(se.rel,
#'   add_row_data= TRUE,
#'   add_col_data = colnames(colData(se.rel)),
#'   abund_values = "counts"
#' )
#' head(molten_se)

#'
NULL
setGeneric("meltAssay",
           signature = "x",
           function(x,
                    add_row_data = NULL,
                    add_col_data = NULL,
                    abund_values = "counts", ...) {
             standardGeneric("meltAssay")
           }
)


#' @rdname meltAssay
#'
#' @importFrom SummarizedExperiment rowData rowData<- colData colData<-
#'
#' @export
#'
setMethod("meltAssay",
          signature = c(x = "SummarizedExperiment"),
          function(x,
                   add_row_data = NULL,
                   add_col_data = NULL,
                   abund_values = "counts", ...) {

            # input check add_col_dat
            MicrobiomeExperiment:::.check_abund_values(abund_values, x)

            # check and enforce reserved names
            #x <- .check_enforce_names(x)

            #FeatureID <- SampleID <- Abundance <- NULL

            dfSE <- .melt_assay(x, abund_values)

            if(!is.null(add_row_data)){

              if(is.logical(add_row_data) && length(add_row_data) == 1L && add_row_data){

                add_row_data <- taxonomyRanks(x)

              } else if (isFALSE(all(add_row_data %in% taxonomyRanks(x)))) {

                stop("Please provide valid column names matching
             those in 'taxonomyRanks(se)'")

              }
              dfSE <- .add_row_data_to_molten_assay(dfSE, x, add_row_data)
            }

            #dfSE <- .add_row_data_to_molten_assay(x, dfSE)

            if(!is.null(add_col_data)){

              if(is.logical(add_col_data) && length(add_col_data) == 1L && add_col_data){

                add_col_data <- colnames(colData(x))

              } else if (isFALSE(all(add_col_data %in% colnames(colData(x))))) {

                stop("Please provide valid column names matching
             those in 'colData(x)'")

              }
              dfSE <- .add_col_data_to_molten_assay(dfSE, x, add_col_data)
            }

            #dfSE <- .add_col_data_to_molten_assay(dfSE, x, add_col_data)
            return(dfSE)
          }
)



#' @description Melts assay in \code{\link[=SummarizedExperiment-class]{SummarizedExperiment}}
#'        object to long data format.
#'
#' @param x A numeric matrix or a \code{\link[=SummarizedExperiment-class]{SummarizedExperiment}}
#'          object containing a tree.
#'
#' @param abund_values Must be one of the values returned by assayNames.
#'
#'
.melt_assay <- function(x, abund_values) {

  # input check
  #.check_abund_values(abund_values, x)

  molten_assay <- assays(x)[[abund_values]] %>%
    data.frame() %>%
    rownames_to_column("FeatureID") %>%
    # SampleID is unique sample id
    pivot_longer(!FeatureID,
                 values_to = "Abundance",
                 names_to = "SampleID"
    )
  return(molten_assay)
}


#' @description Combines molten assay with rowData i.e. taxonomy table.
#' @param dfSE Molten SE.
#' @param x \code{\link[=SummarizedExperiment-class]{SummarizedExperiment}}
#'          object.
#'
.add_row_data_to_molten_assay <- function(dfSE, x, add_row_data) {


  # since it can be a vector consisting of more than one colnames(colData(x))
  # values, choose only 1
  if (is.null(add_row_data[1]) | is.na(add_row_data[1])) {
    message("Returning molten assay without colData")

    dfSE <- dfSE

  } else if (isTRUE(all(add_row_data %in% colnames(rowData(x))))) {

    me_tax <- .get_row_data_frame(x=x)

    me_tax <- me_tax[,c("FeatureID", add_row_data)]

    dfSE <- dfSE %>%
      left_join(me_tax, by = "FeatureID") %>%
      mutate_if(is.factor, as.character)
  }

  return(dfSE)
}


#' @description Converts rowData to data.frame to avoid issues with
#'      subsetting and rownames.
#'
#' @param x \code{\link[=SummarizedExperiment-class]{SummarizedExperiment}}
#'          object.
#'
#' @return data.frame
#'
.get_row_data_frame <- function(x){

  if(isTRUE(any(taxonomyRanks(x)=="FeatureID"))){
    message("To avoid issues with SE plotting
            existing `FeatureID` column has been
            renamed as `OldFeatureID` ")

    me_tax <- SummarizedExperiment::rowData(x) %>%
      data.frame() %>%
      rename(OldFeatureID=FeatureID) %>%
      rownames_to_column("FeatureID")

    return(me_tax)

  } else {

    me_tax <- SummarizedExperiment::rowData(x) %>%
      data.frame() %>%
      rownames_to_column("FeatureID")

    return(me_tax)

  }
}


#' @description Combines molten assay and rowData i.e. taxonomy table with
#'              ColData i.e. sample information.
#' @param dfSE2 Molten SE.
#' @param x \code{\link[=SummarizedExperiment-class]{SummarizedExperiment}}
#'          object.
#' @param add_col_data Choice of colData columns to add. Default="none".
#'
.add_col_data_to_molten_assay <- function(dfSE, x, add_col_data) {

  #dfSE3 <- SampleID <- NULL

  #if(add_col_data[1] == "none") {
  # add_col_data <- NULL
  #}

  # since it can be a vector consisting of more than one colnames(colData(x))
  # values, choose only 1
  if (is.null(add_col_data[1]) | is.na(add_col_data[1])) {
    message("Returning molten assay without colData")

    dfSE <- dfSE

  } else if (isTRUE(all(add_col_data %in% colnames(colData(x))))) {

    me_sam <- .get_col_data_frame(x)
    me_sam <- me_sam[,c("SampleID", add_col_data)]

    dfSE <- dfSE %>%
      left_join(me_sam, by = "SampleID") %>%
      mutate_if(is.factor, as.character)
  }

  return(dfSE)
}


#' @description Converts colData to data.frame to avoid issues with
#'      subsetting and rownames.
#'
#' @param x \code{\link[=SummarizedExperiment-class]{SummarizedExperiment}}
#'          object.
#'
#' @return data.frame
#'
.get_col_data_frame <- function(x){

  if(isTRUE(any(colnames(colData(x))=="SampleID"))){

    message("To avoid issues with SE plotting
            existing `SampleID` column has been
            renamed as `OldSampleID` ")

    me_sam <- colData(x) %>%
      data.frame() %>%
      rename(OldSampleID=SampleID) %>%
      rownames_to_column("SampleID")
    return(me_sam)

  } else {

    me_sam <- colData(x) %>%
      data.frame() %>%
      rownames_to_column("SampleID")

    return(me_sam)

  }
}







