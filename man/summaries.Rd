% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/summaries.R
\name{summaries}
\alias{summaries}
\alias{getTopTaxa}
\alias{getTopTaxa,SummarizedExperiment-method}
\alias{getUniqueTaxa}
\alias{getUniqueTaxa,SummarizedExperiment-method}
\alias{countDominantTaxa}
\alias{countDominantTaxa,SummarizedExperiment-method}
\alias{summary,SummarizedExperiment-method}
\title{Summarizing microbiome data}
\usage{
getTopTaxa(
  x,
  top = 5L,
  method = c("mean", "sum", "median"),
  abund_values = "counts"
)

\S4method{getTopTaxa}{SummarizedExperiment}(
  x,
  top = 5L,
  method = c("mean", "sum", "median", "prevalence"),
  abund_values = "counts"
)

getUniqueTaxa(x, ...)

\S4method{getUniqueTaxa}{SummarizedExperiment}(x, rank = NULL)

countDominantTaxa(x, group = NULL, ...)

\S4method{countDominantTaxa}{SummarizedExperiment}(x, group = NULL, ...)

\S4method{summary}{SummarizedExperiment}(object, abund_values = "counts")
}
\arguments{
\item{x}{A
\code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}} object.}

\item{top}{Numeric value, how many top taxa to return. Default return top
five taxa.}

\item{method}{Specify the method to determine top taxa. Either sum, mean,
median or prevalence. Default is 'mean'.}

\item{abund_values}{a \code{character} value to select an
\code{\link[SummarizedExperiment:SummarizedExperiment-class]{assayNames}}
By default it expects count data.}

\item{...}{Additional arguments passed on to \code{agglomerateByRank()} when
\code{rank} is specified for \code{countDominantTaxa}.}

\item{rank}{A single character defining a taxonomic rank. Must be a value of
the output of \code{taxonomicRanks()}.}

\item{group}{With group, it is possible to group the observations in an
overview. Must be one of the column names of \code{colData}.}

\item{object}{A
\code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}} object.}
}
\value{
The \code{getTopTaxa} returns a vector of the most \code{top} abundant
\dQuote{FeatureID}s

The \code{getUniqueTaxa} returns a vector of unique taxa present at a
particular rank

The \code{countDominantTaxa} returns an overview in a tibble. It contains dominant taxa
in a column named \code{*name*} and its abundance in the data set.

The \code{summary} returns a list with two \code{tibble}s
}
\description{
To query a \code{SummarizedExperiment} for interesting features, several
functions are available.
}
\details{
The \code{getTopTaxa} extracts the most \code{top} abundant \dQuote{FeatureID}s
in a \code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
object.

The \code{getUniqueTaxa} is a basic function to access different taxa at a
particular taxonomic rank.

\code{countDominantTaxa} returns information about most dominant
taxa in a tibble. Information includes their absolute and relative
abundances in whole data set.

The \code{summary} will return a summary of counts for all samples and
features in
\code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
object.
}
\examples{
data(GlobalPatterns)
top_taxa <- getTopTaxa(GlobalPatterns,
                       method = "mean",
                       top = 5,
                       abund_values = "counts")
top_taxa

# Gets the overview of dominant taxa
dominant_taxa <- countDominantTaxa(GlobalPatterns,
                                   rank = "Genus")
dominant_taxa

# With group, it is possible to group observations based on specified groups
# Gets the overview of dominant taxa
dominant_taxa <- countDominantTaxa(GlobalPatterns,
                                   rank = "Genus",
                                   group = "SampleType",
                                   na.rm= TRUE)

dominant_taxa

# Get an overview of sample and taxa counts
summary(GlobalPatterns)

# Get unique taxa at a particular taxonomic rank
getUniqueTaxa(GlobalPatterns, "Phylum")

}
\seealso{
\code{\link[=getPrevalence]{getPrevalentTaxa}}

\code{\link[scuttle:perCellQCMetrics]{perCellQCMetrics}},
\code{\link[scuttle:perFeatureQCMetrics]{perFeatureQCMetrics}},
\code{\link[scuttle:addPerCellQC]{addPerCellQC}},
\code{\link[scuttle:addPerFeatureQC]{addPerFeatureQC}},
\code{\link[scuttle:quickPerCellQC]{quickPerCellQC}}
}
\author{
Leo Lahti, Tuomas Borman and Sudarshan A. Shetty
}
