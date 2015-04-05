\name{RoarDataset-class}
\Rdversion{1.1}
\docType{class}
\alias{RoarDataset-class}
\alias{computePvals,RoarDataset-method}
\alias{computePairedPvals,RoarDataset,numeric,numeric-method}
\alias{computeRoars,RoarDataset-method}
\alias{cores,RoarDataset-method}
\alias{countPrePost,RoarDataset,logical-method}
\alias{countPrePost,RoarDataset-method}
\alias{fpkmResults,RoarDataset-method}
\alias{countResults,RoarDataset-method}
\alias{totalResults,RoarDataset-method}
\alias{standardFilter,RoarDataset,numeric-method}
\alias{pvalueFilter,RoarDataset,numeric,numeric-method}
\alias{pvalueCorrectFilter,RoarDataset,numeric,numeric,character-method}

\title{Class \code{"RoarDataset"}}
\description{
RoarDataset - a class to perform 3'UTR shortening analyses
}
\section{Objects from the Class}{
Objects of thiss class should be created using the functions
\code{\link{RoarDataset}} or \code{\link{RoarDatasetFromFiles}}, ideally
the raw \code{\link{new}} method should never be invoked by end users. Then
to perform the analysis the user should call, in order: countPrePost, computeRoars, computePvals and
one of the methods to format results.
}
\section{Slots}{
  \describe{
    \item{\code{treatmentBams}:}{Object of class \code{"list"} - a list of GappedAlignment objects for the first condition (by convention it is considered the ``treated'' condition) in analysis. }
    \item{\code{controlBams}:}{Object of class \code{"list"} - a list of GappedAlignment objects for the second condition (by convention it is considered the ``control'' condition) in analysis. }
    \item{\code{prePostCoords}:}{Object of class \code{"GRanges"} - represents the APA sites coords, defining "PRE" (last exon coords up until the alternative APA, defining the shorter isoform) and "POST" (from the alternative APA to the ``standard'' one) regions of the genes. }
    \item{\code{postCoords}:}{Object of class \code{"GRanges"} - private object. }
    \item{\code{countsTreatment}:}{Object of class \code{"SummarizedExperiment"} - private object. }
    \item{\code{countsControl}:}{Object of class \code{"SummarizedExperiment"} - private object. }
    \item{\code{pVals}:}{Object of class \code{"SummarizedExperiment"} - private object. }
    \item{\code{paired}:}{\code{"logical"} slot - private. }
    \item{\code{step}:}{\code{"numeric"} slot - private. }
    \item{\code{cores}:}{\code{"numeric"} slot - private. }
    \item{\code{exptData}:}{Object of class \code{"SimpleList"} - private object. }
    \item{\code{rowData}:}{Object of class \code{"GRangesORGRangesList"} - private object. }
    \item{\code{colData}:}{Object of class \code{"DataFrame"} - private object. }
    \item{\code{assays}:}{Object of class \code{"Assays"} - private object. }
  }
}
\section{Extends}{
Class \code{"\linkS4class{SummarizedExperiment}"}, directly.
}
\section{Methods}{
  \describe{
    \item{\code{\link{countPrePost}}}{\code{signature(rds = "RoarDataset", stranded = "logical")}: Counts reads falling over PRE/POST portions of the given transcripts. }
    \item{\code{\link{computeRoars}}}{\code{signature(rds = "RoarDataset")}: Computes m/M and roar values for this \code{\linkS4class{RoarDataset}} object. }
    \item{\code{\link{computePvals}}}{\code{signature(rds = "RoarDataset")}: Computes pvalues (Fisher test) for this \code{\linkS4class{RoarDataset}} object. }
    \item{totalResults}{\code{signature(rds = "RoarDataset")}: Returns a dataframe with results of the analysis for a \code{\linkS4class{RoarDataset}} object. }
    \item{\code{\link{fpkmResults}}}{\code{signature(rds = "RoarDataset")}: 
    The last step of a classical Roar analyses: it returns a dataframe containing m/M values, roar 
    values, pvalues and estimates of expression (a measure recalling FPKM). }
    \item{\code{\link{countResults}}}{\code{signature(rds = "RoarDataset")}: 
    The last step of a classical Roar analyses: it returns a dataframe containing m/M values, roar 
    values, pvalues and estimates of expression (counts over PRE portions). }
    \item{\code{\link{standardFilter}}}{\code{signature(rds = "RoarDataset", fpkmCutoff = "double")}: Returns a dataframe with results of the analysis for a \code{\linkS4class{RoarDataset}} object. }
    \item{\code{\link{pvalueFilter}}}{\code{signature(rds = "RoarDataset", fpkmCutoff = "double", pvalCutoff = "double")}: ... }
    \item{\code{\link{cores}}}{\code{signature(rds = "RoarDataset")}: returns the number of cores used for computation, right now always 1. }
    }
}
%\references{
% website for now? XXX
%}
\author{
Elena Grassi, PhD student in  Biomedical Sciences and Oncology - Dept. of Molecular Biotechnologies and Health Sciences, Molecular Biotechnology Center, Torino
}

\examples{
showClass("RoarDataset")
}
\keyword{classes}