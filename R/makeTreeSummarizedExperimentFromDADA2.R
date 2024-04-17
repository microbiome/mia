.make_TreeSE_from_DADA2 <- function(...) {
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
