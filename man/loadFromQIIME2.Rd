% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/loadFromQIIME2.R
\name{loadFromQIIME2}
\alias{loadFromQIIME2}
\title{Import QIIME2 results to \code{TreeSummarizedExperiment}}
\usage{
loadFromQIIME2(
  featureTableFile,
  taxonomyTableFile = NULL,
  sampleMetaFile = NULL,
  featureNamesAsRefSeq = TRUE,
  refSeqFile = NULL,
  phyTreeFile = NULL,
  ...
)
}
\arguments{
\item{featureTableFile}{a single \code{character} value defining the file
path of the feature table to be imported.}

\item{taxonomyTableFile}{a single \code{character} value defining the file
path of the taxonomy table to be imported. (default:
\code{taxonomyTableFile = NULL}).}

\item{sampleMetaFile}{a single \code{character} value defining the file path
of the sample metadata to be imported. The file has to be in tsv format.
(default: \code{sampleMetaFile = NULL}).}

\item{featureNamesAsRefSeq}{\code{TRUE} or \code{FALSE}: Should the feature
names of the feature table be regarded as reference sequences? This setting
will be disregarded, if \code{refSeqFile} is not \code{NULL}. If the
feature names do not contain valid DNA characters only, the reference
sequences will not be set.}

\item{refSeqFile}{a single \code{character} value defining the file path of
the reference sequences for each feature. (default: \code{refSeqFile =
  NULL}).}

\item{phyTreeFile}{a single \code{character} value defining the file path of
the phylogenetic tree. (default: \code{phyTreeFile = NULL}).}

\item{...}{additional arguments:
\itemize{
\item{\code{temp}:} {the temporary directory used for decompressing the
data. (default: \code{tempdir()})}
\item{\code{removeTaxaPrefixes}:} {\code{TRUE} or \code{FALSE}: Should
taxonomic prefixes be removed? (default:
\code{removeTaxaPrefixes = FALSE})}
}}
}
\value{
A
\code{\link[TreeSummarizedExperiment:TreeSummarizedExperiment-class]{TreeSummarizedExperiment}}
object
}
\description{
Results exported from QIMME2 can be imported as a
\code{TreeSummarizedExperiment} using \code{loadFromQIIME2}. Except for the
\code{featureTableFile}, the other data types, \code{taxonomyTableFile},
\code{refSeqFile} and \code{phyTreeFile}, are optional, but are highly
encouraged to be provided.
}
\details{
Both arguments \code{featureNamesAsRefSeq} and \code{refSeqFile} can be used
to define reference sequences of features. \code{featureNamesAsRefSeq} is
only taken into account, if \code{refSeqFile} is \code{NULL}. No reference
sequences are tried to be created, if \code{featureNameAsRefSeq} is
\code{FALSE} and \code{refSeqFile} is \code{NULL}.
}
\examples{
featureTableFile <- system.file("extdata", "table.qza", package = "mia")
taxonomyTableFile <- system.file("extdata", "taxonomy.qza", package = "mia")
sampleMetaFile <- system.file("extdata", "sample-metadata.tsv", package = "mia")
phyTreeFile <- system.file("extdata", "tree.qza", package = "mia")
refSeqFile <- system.file("extdata", "refseq.qza", package = "mia")
tse <- loadFromQIIME2(
  featureTableFile = featureTableFile,
  taxonomyTableFile = taxonomyTableFile,
  sampleMetaFile = sampleMetaFile,
  refSeqFile = refSeqFile,
  phyTreeFile = phyTreeFile
)

tse
}
\references{
Bolyen E et al. 2019: Reproducible, interactive, scalable and extensible
microbiome data science using QIIME 2. Nature Biotechnology 37: 852–857.
\url{https://doi.org/10.1038/s41587-019-0209-9}

\url{https://qiime2.org}
}
\seealso{
\code{\link[=makeTreeSummarizedExperimentFromphyloseq]{makeTreeSummarizedExperimentFromphyloseq}}
\code{\link[=makeSummarizedExperimentFromBiom]{makeSummarizedExperimentFromBiom}}
\code{\link[=makeTreeSummarizedExperimentFromDADA2]{makeTreeSummarizedExperimentFromDADA2}}
\code{\link[=loadFromMothur]{loadFromMothur}}
}
\author{
Yang Cao
}
