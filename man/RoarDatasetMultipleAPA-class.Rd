\name{RoarDatasetMultipleAPA-class}
\Rdversion{1.1}
\docType{class}
\alias{RoarDatasetMultipleAPA-class}
\alias{computePvals,RoarDatasetMultipleAPA-method}
\alias{computePairedPvals,RoarDatasetMultipleAPA,numeric,numeric-method}
\alias{computeRoars,RoarDatasetMultipleAPA-method}
\alias{cores,RoarDatasetMultipleAPA-method}
\alias{countPrePost,RoarDatasetMultipleAPA,logical-method}
\alias{countPrePost,RoarDatasetMultipleAPA-method}
\alias{fpkmResults,RoarDatasetMultipleAPA-method}
\alias{countResults,RoarDatasetMultipleAPA-method}
\alias{totalResults,RoarDatasetMultipleAPA-method}
\alias{standardFilter,RoarDatasetMultipleAPA,numeric-method}
\alias{pvalueFilter,RoarDatasetMultipleAPA,numeric,numeric-method}
\alias{pvalueCorrectFilter,RoarDatasetMultipleAPA,numeric,numeric,character-method}

\title{Class \code{"RoarDatasetMultipleAPA"}}
\description{
RoarDataset - a class to perform 3'UTR shortening analyses
}
\section{Objects from the Class}{
Objects of thiss class should be created using the functions
\code{\link{RoarDatasetMultipleAPA}} or \code{\link{RoarDatasetMultipleAPAFromFiles}}, ideally
the raw \code{\link{new}} method should never be invoked by end users. Then
to perform the analysis the user should call, in order: countPrePost, computeRoars, computePvals and
one of the methods to format results. This class is used to allow
efficient analyses that allow to study more than one APA site for each gene: internally
it uses a \code{\link{RoarDataset}} object that stores PRE/POST counts for all possible alternative
APA choices for each gene.
}
\section{Slots}{
  \describe{
    \item{\code{treatmentBams}:}{Object of class \code{"list"} - a list of GappedAlignment objects for the first condition (by convention it is considered the ``treated'' condition) in analysis. }
    \item{\code{controlBams}:}{Object of class \code{"list"} - a list of GappedAlignment objects for the second condition (by convention it is considered the ``control'' condition) in analysis. }
    \item{\code{geneCoords}:}{Object of class \code{"GRangesList"} - private object that represents the exon structures of genes in study. }
    \item{\code{apaCoords}:}{Object of class \code{"GRangesList"} - private object that represents the APA fallin on genes in study.  }
    \item{\code{fragments}:}{Object of class \code{"GRangesList"} - private object used to efficiently count reads falling on short and long isoforms. }
    \item{\code{prePostDef}:}{Object of class \code{"list"} - private object representing all possible short and long isoforms. }
    \item{\code{roars}:}{Object of class \code{"list"} - private object with a list of \code{\link{RoarDataset}} objects, each one representing all possible PRE/POST choices for a single gene. }
    \item{\code{corrTreatment}:}{\code{"numeric"} slot - private, integer representing the mean length of reads for the treatment samples. }
    \item{\code{corrControl}:}{\code{"numeric"} slot - private, integer representing the mean length of reads for the control samples. }
    \item{\code{paired}:}{\code{"logical"} slot - private. }
    \item{\code{step}:}{\code{"numeric"} slot - private. }
    \item{\code{cores}:}{\code{"numeric"} slot - private. }
  }
}

\section{Methods}{
  \describe{
    \item{\code{\link{countPrePost}}}{\code{signature(rds = "RoarDatasetMultipleAPA", stranded = "logical")}: Counts reads falling over all the possible PRE/POST portions of the given transcripts. WARNING: stranded = TRUE is still unsupported and could give unpredictable results. }
    \item{\code{\link{computeRoars}}}{\code{signature(rds = "RoarDatasetMultipleAPA")}: Computes m/M and roar values for this \code{\linkS4class{RoarDatasetMultipleAPA}} object. }
    \item{\code{\link{computePvals}}}{\code{signature(rds = "RoarDatasetMultipleAPA")}: Computes pvalues (Fisher test) for this \code{\linkS4class{RoarDatasetMultipleAPA}} object. }
    \item{totalResults}{\code{signature(rds = "RoarDatasetMultipleAPA")}: Returns a dataframe with results of the analysis for a \code{\linkS4class{RoarDatasetMultipleAPA}} object. }
    \item{\code{\link{fpkmResults}}}{\code{signature(rds = "RoarDatasetMultipleAPA")}: 
    The last step of a classical Roar analyses: it returns a dataframe containing m/M values, roar 
    values, pvalues and estimates of expression (a measure recalling FPKM). }
    \item{\code{\link{countResults}}}{\code{signature(rds = "RoarDatasetMultipleAPA")}: 
    The last step of a classical Roar analyses: it returns a dataframe containing m/M values, roar 
    values, pvalues and estimates of expression (counts of reads falling over a gene). }
    \item{\code{\link{standardFilter}}}{\code{signature(rds = "RoarDatasetMultipleAPA", fpkmCutoff = "double")}: Returns a dataframe with results of the analysis for a \code{\linkS4class{RoarDatasetMultipleAPA}} object. }
    \item{\code{\link{pvalueFilter}}}{\code{signature(rds = "RoarDatasetMultipleAPA", fpkmCutoff = "double", pvalCutoff = "double")}: ... }
    \item{\code{\link{cores}}}{\code{signature(rds = "RoarDatasetMultipleAPA")}: returns the number of cores used for computation, right now always 1. }
    }
}

\author{
Elena Grassi, PhD student in  Biomedical Sciences and Oncology - Dept. of Molecular Biotechnologies and Health Sciences, Molecular Biotechnology Center, Torino
}

\examples{
showClass("RoarDatasetMultipleAPA")
}
\keyword{classes}
