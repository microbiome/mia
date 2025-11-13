#' IBDMDB 2-omic demo (MGX + MTX)
#'
#' Compact example objects prepared from IBDMDB/HMP2 for quick examples/vignettes:
#' \itemize{
#'   \item \code{se_mgx}: MGX \link[SummarizedExperiment]{SummarizedExperiment}
#'   \item \code{se_mtx}: MTX \link[SummarizedExperiment]{SummarizedExperiment}
#'   \item \code{mae2}:   \link[MultiAssayExperiment]{MultiAssayExperiment} with MGX and MTX assays
#' }
#'
#' Load with \code{data(ibdmdb_2omic_demo)}; this will place the objects
#' \code{se_mgx}, \code{se_mtx}, and \code{mae2} into your workspace.
#'
#' @name ibdmdb_2omic_demo
#' @docType data
#' @usage data(ibdmdb_2omic_demo)
#' @keywords datasets
#' @format An .rda file containing \code{se_mgx}, \code{se_mtx}, \code{mae2}.
#' @source Prepared by \code{inst/scripts/prepare_ibdmdb_demo.R}.
NULL

#' IBDMDB metadata (demo subset)
#'
#' Optional IBDMDB metadata used with the demo objects.
#' Load with \code{data(ibdmdb_meta_demo)}; this will place the object
#' \code{ibdmdb_meta_demo} (a \code{data.frame}).
#'
#' @name ibdmdb_meta_demo
#' @docType data
#' @usage data(ibdmdb_meta_demo)
#' @keywords datasets
#' @format A \code{data.frame}.
#' @source Prepared by \code{inst/scripts/prepare_ibdmdb_demo.R}.
NULL

#' IBDMDB 3-omic demo (16S + MGX + MTX)
#'
#' Compact example objects prepared from IBDMDB/HMP2:
#' \itemize{
#'   \item \code{se_16s}, \code{se_mgx2}, \code{se_mtx2}: \link[SummarizedExperiment]{SummarizedExperiment}s
#'   \item \code{mae3}: \link[MultiAssayExperiment]{MultiAssayExperiment} with all three assays
#' }
#'
#' Load with \code{data(ibdmdb_3omic_demo)}; this will place the objects
#' \code{se_16s}, \code{se_mgx2}, \code{se_mtx2}, and \code{mae3}.
#'
#' @name ibdmdb_3omic_demo
#' @docType data
#' @usage data(ibdmdb_3omic_demo)
#' @keywords datasets
#' @format An .rda file containing \code{se_16s}, \code{se_mgx2}, \code{se_mtx2}, \code{mae3}.
#' @source Prepared by \code{inst/scripts/prepare_ibdmdb_demo.R}.
NULL