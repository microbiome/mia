#' Convert data from various formats to a TreeSummarizedExperiment object or 
#' Phyloseq object
#'
#' @param x Input data to be converted. See details for more information on 
#' supported formats.
#' 
#' @param ... See \code{mergePairs} function for more details 
#' (x : dada2 object)
#' 
#' @param removeTaxaPrefixes \code{TRUE} or \code{FALSE}: Should
#' taxonomic prefixes be removed? The prefixes is removed only from detected
#' taxa columns meaning that \code{rankFromPrefix} should be enabled in the most cases.
#' (default \code{removeTaxaPrefixes = FALSE}) 
#' (x : biom object)
#' 
#' @param rankFromPrefix \code{TRUE} or \code{FALSE}: If file does not have
#' taxonomic ranks on feature table, should they be scraped from prefixes?
#' (default \code{rankFromPrefix = FALSE}) 
#' (x : biom object)
#' 
#' @param remove.artifacts \code{TRUE} or \code{FALSE}: If file have
#' some taxonomic character naming artifacts, should they be removed.
#' (default \code{remove.artifacts = FALSE}) 
#' (x : biom object)
#' 
#' @param ... additional arguments 
#'   \itemize{
#'        \item{\code{patter}}{\code{character} value specifying artifacts
#'        to be removed. If \code{patterns = "auto"}, special characters
#'        are removed. (default: \code{pattern = "auto"}) 
#'        (x : biom object)
#'        } 
#'    }
#'
#' @param assay.type A single character value for selecting the
#'   \code{\link[SummarizedExperiment:SummarizedExperiment-class]{assay}} to be
#'   included in the phyloseq object that is created. 
#'   (By default: \code{assay.type = "counts"})
#'   (x : phyloseq object)
#'   
#' @param assay_name a single \code{character} value for specifying which
#'   assay to use for calculation.
#'   (Please use \code{assay.type} instead. At some point \code{assay_name}
#'   will be disabled.)
#'   (x : phyloseq object)
#'   
#' @param tree_name a single \code{character} value for specifying which
#'   tree will be included in the phyloseq object that is created, 
#'   (By default: \code{tree_name = "phylo"})
#'   (x : phyloseq object)
#'
#' @param ... Additional arguments
#'
#' @details
#' The `convert` function supports the conversion of data from the following
#' formats:
#' 
#' - Coerce a \code{phyloseq} object to a \code{TreeSummarizedExperiment}
#' - Coerce \sQuote{DADA2} results to a \code{TreeSummarizedExperiment}
#' - Coerce a \sQuote{biom} object to a \code{TreeSummarizedExperiment}
#' - Coerce a \code{TreeSummarizedExperiment} object to a \code{phyloseq} object
#' 
#' The appropriate conversion is automatically selected based on the class of 
#' the input data.
#'
#' @return An object of class TreeSummarizedExperiment or Phyloseq.
#'
#' @importFrom Biostrings DNAStringSet
#' @importFrom S4Vectors SimpleList DataFrame make_zero_col_DFrame
#' @importFrom SummarizedExperiment colData colData<-
#'
#' @name convert
#' 
#' @seealso 
#' \code{\link[=importQIIME2]{importQIIME2}}
#' \code{\link[=importMothur]{importMothur}}
#'
#' @examples
#' ### Make a TreeSummarizedExperiment from a Phyloseq object
#' 
#' if( requireNamespace("phyloseq") ) {
#'     data(GlobalPatterns, package="phyloseq")
#'     convert(GlobalPatterns)
#'     data(enterotype, package="phyloseq")
#'     convert(enterotype)
#'     data(esophagus, package="phyloseq")
#'     convert(esophagus)
#' }
#' 
#' ### Make a TreesummarizedExperiment from dada2 results
#' 
#' if( requireNamespace("dada2") ) {
#'     fnF <- system.file("extdata", "sam1F.fastq.gz", package="dada2")
#'     fnR = system.file("extdata", "sam1R.fastq.gz", package="dada2")
#'     dadaF <- dada2::dada(fnF, selfConsist=TRUE)
#'     dadaR <- dada2::dada(fnR, selfConsist=TRUE)
#'
#'     tse <- convert(dadaF, fnF, dadaR, fnR)
#'     tse
#' }
#' 
#' ### Make a TreeSummarizedExperiment from a biom file
#' 
#' # Load biom file
#' library(biomformat)
#' biom_file <- system.file("extdata", "rich_dense_otu_table.biom",
#'                          package = "biomformat")
#' 
#' # Make a TreeSummarizedExperiment from a biom object
#' biom_object <- biomformat::read_biom(biom_file)
#' tse <- convert(biom_object)
#' 
#' ### Make a Phyloseq object from a TreeSummarizedExperiment
#' 
#' # Get tse object
#' data(GlobalPatterns)
#' tse <- GlobalPatterns
#'
#' # Create a phyloseq object from it
#' phy <- convert(tse)
#' phy
#'
#' # By default the chosen table is counts, but if there are other tables,
#' # they can be chosen with assay.type.
#'
#' # Counts relative abundances table
#' tse <- transformAssay(tse, method = "relabundance")
#' phy2 <- convert(tse, assay.type = "relabundance")
#' phy2
#' @export
NULL

#' @rdname convert
#' @export
setGeneric("convert", signature = c("x"),
           function(x,...)
               standardGeneric("convert"))

#' @export
setMethod("convert", signature = c(x = "SummarizedExperiment"),
          function(x,...){
              .make_phyloseq_from_SE(x,...)
          }
)

#' @export
setMethod("convert", signature = c(x = "TreeSummarizedExperiment"),
          function(x,...){
              .make_phyloseq_from_TreeSE(x,...)
          }
)

#' @export
setMethod("convert", signature = c(x = "dada"),
          function(x,...){
              .make_TreeSE_from_DADA2(x,...)
          }
)

#' @export
setMethod("convert", signature = c(x = "phyloseq"),
          function(x,...){
              .make_TreeSE_from_phyloseq(x,...)
          }
)

#' @export
setMethod("convert", signature = c(x = "biom"),
          function(x,...){
              .make_TreeSE_from_biom(x,...)
          }
)
