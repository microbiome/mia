#' IBDMDB 2-omic demo dataset (MGX + MTX)
#'
#' A compact example derived from the Integrative Human Microbiome Project (iHMP)
#' Inflammatory Bowel Disease (IBD) cohort.  
#'
#' This demo contains two experiments:
#' \itemize{
#'   \item \code{se_mgx}: metagenomic taxonomic profiles (MGX)
#'   \item \code{se_mtx}: metatranscriptomic taxonomic profiles (MTX)
#'   \item \code{mae2}: a
#'     \link[MultiAssayExperiment]{MultiAssayExperiment}
#'     containing both experiments
#' }
#'
#' These compact objects are intended for quick examples and vignettes.
#' Load with:
#' \code{data(ibdmdb_2omic_demo)}.
#'
#' @name ibdmdb_2omic_demo
#' @docType data
#' @usage data(ibdmdb_2omic_demo)
#' @keywords datasets
#'
#' @format A `.rda` file containing:
#' \describe{
#'   \item{se_mgx}{A \code{SummarizedExperiment} with MGX abundance data.}
#'   \item{se_mtx}{A \code{SummarizedExperiment} with MTX abundance data.}
#'   \item{mae2}{A \code{MultiAssayExperiment} with two experiments: MGX and MTX.}
#' }
#'
#' @source Prepared by \code{inst/scripts/prepare_ibdmdb_demo.R} from
#'   publicly available IBDMDB/HMP2 taxonomic profiles.
#'
#' @references
#' Lloyd-Price J, Arze C, Ananthakrishnan AN, et al.  
#' \emph{Multi-omics of the gut microbial ecosystem in inflammatory bowel diseases.}  
#' Nature 569, 655–662 (2019).  
#' \doi{10.1038/s41586-019-1237-9}
NULL

#' IBDMDB demo metadata subset
#'
#' Sample-level metadata corresponding to the compact IBDMDB 2-omic demo
#' (\code{ibdmdb_2omic_demo}). This is essentially
#' \code{as.data.frame(colData(mae2))} and is provided as a convenience
#' for examples and vignettes.
#'
#' Loaded with:
#' \code{data(ibdmdb_meta_demo)}.
#'
#' @name ibdmdb_meta_demo
#' @docType data
#' @usage data(ibdmdb_meta_demo)
#' @keywords datasets
#'
#' @format A \code{data.frame} with selected sample metadata fields.
#'
#' @source Prepared by \code{inst/scripts/prepare_ibdmdb_demo.R} from
#'   publicly available iHMP/IBDMDB metadata.
#'
#' @references
#' Lloyd-Price J, Arze C, Ananthakrishnan AN, et al.  
#' \emph{Multi-omics of the gut microbial ecosystem in inflammatory bowel diseases.}  
#' Nature 569, 655–662 (2019).  
#' \doi{10.1038/s41586-019-1237-9}
NULL