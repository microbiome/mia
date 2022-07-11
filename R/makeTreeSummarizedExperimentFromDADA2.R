#' Coerce \sQuote{DADA2} results to \code{TreeSummarizedExperiment}
#'
#' \code{makeTreeSEFromDADA2} is a wrapper for the
#' \code{mergePairs} function from the \code{dada2} package.
#'
#' @param ... See \code{mergePairs} function for
#'   more details.
#'
#' @details
#' A count matrix is contructed via \code{makeSequenceTable(mergePairs(...))}
#' and rownames are dynamically created as \code{ASV(N)} with \code{N} from
#' 1 to \code{nrow} of the count tables. The colnames and rownames from the
#' output of \code{makeSequenceTable} are stored as \code{colnames} and in the
#' \code{referenceSeq} slot of the \code{TreeSummarizedExperiment},
#' respectively.
#'
#' @return An object of class \code{TreeSummarizedExperiment}
#'
#' @importFrom S4Vectors SimpleList
#' @importFrom Biostrings DNAStringSet
#'
#' @name makeTreeSEFromDADA2
#' @seealso
#' \code{\link[=makeTreeSEFromPhyloseq]{makeTreeSEFromPhyloseq}}
#' \code{\link[=makeTreeSEFromBiom]{makeTreeSEFromBiom}}
#' \code{\link[=loadFromQIIME2]{loadFromQIIME2}}
#' \code{\link[=loadFromMothur]{loadFromMothur}}
#'
#' @export
#'
#' @examples
#' if(requireNamespace("dada2")) {
#'   fnF <- system.file("extdata", "sam1F.fastq.gz", package="dada2")
#'   fnR = system.file("extdata", "sam1R.fastq.gz", package="dada2")
#'   dadaF <- dada2::dada(fnF, selfConsist=TRUE)
#'   dadaR <- dada2::dada(fnR, selfConsist=TRUE)
#'
#'   tse <- makeTreeSEFromDADA2(dadaF, fnF, dadaR, fnR)
#'   tse
#' }
makeTreeSEFromDADA2 <- function(...) {
    # input checks
    .require_package("dada2")
    .require_package("stringr")
    #
    mergers <- dada2::mergePairs(...)
    seqtab <- dada2::makeSequenceTable(mergers)
    seqtab <- t(seqtab)
    # generate row and col names
    rName <- paste0("ASV",
                    stringr::str_pad(seq.int(1L,nrow(seqtab)),
                                     nchar(nrow(seqtab)) + 1L,
                                     pad="0"))
    cName <- colnames(seqtab)
    # retrieve count data and reference sequence
    assays <- S4Vectors::SimpleList(counts = unname(seqtab))
    refseq <- Biostrings::DNAStringSet(rownames(seqtab))
    # construct ME an name rows and cols
    output <- TreeSummarizedExperiment(assays = assays,
                                       referenceSeq = refseq)
    colnames(output) <- cName
    rownames(output) <- rName
    output
}
