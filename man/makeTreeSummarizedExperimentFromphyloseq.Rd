% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/makeTreeSummarizedExperimentFromphyloseq.R
\name{makeTreeSummarizedExperimentFromphyloseq}
\alias{makeTreeSummarizedExperimentFromphyloseq}
\title{Coerce a \code{phyloseq} object to a \code{TreeSummarizedExperiment}}
\usage{
makeTreeSummarizedExperimentFromphyloseq(obj)
}
\arguments{
\item{obj}{a \code{phyloseq} object}
}
\value{
An object of class \code{TreeSummarizedExperiment}
}
\description{
\code{makeTreeSummarizedExperimentFromphyloseq} converts \code{phyloseq}
objects into \code{TreeSummarizedExperiment} objects.
}
\details{
All data stored in a \code{phyloseq} object is transfered.
}
\examples{
if (requireNamespace("phyloseq")) {
    data(GlobalPatterns, package="phyloseq")
    makeTreeSummarizedExperimentFromphyloseq(GlobalPatterns)
    data(enterotype, package="phyloseq")
    makeTreeSummarizedExperimentFromphyloseq(enterotype)
    data(esophagus, package="phyloseq")
    makeTreeSummarizedExperimentFromphyloseq(esophagus)
}
}
\seealso{
\code{\link[=makeSummarizedExperimentFromBiom]{makeSummarizedExperimentFromBiom}}
\code{\link[=makeTreeSummarizedExperimentFromDADA2]{makeTreeSummarizedExperimentFromDADA2}}
\code{\link[=loadFromQIIME2]{loadFromQIIME2}}
\code{\link[=loadFromMothur]{loadFromMothur}}
}
