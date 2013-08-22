\name{RoarDataset-class}
\Rdversion{1.1}
\docType{class}
\alias{RoarDataset-class}
\alias{computePvals,RoarDataset-method}
\alias{computeRoars,RoarDataset-method}
\alias{cores<-,RoarDataset,numeric-method}
\alias{cores,RoarDataset-method}
\alias{countPrePost,RoarDataset,logical-method}
\alias{filteringInfoResults,RoarDataset-method}
\alias{totalResults,RoarDataset-method}

\title{Class \code{"RoarDataset"}}
\description{
RoarDataset - a class to perform 3'UTR shortening analyses
}
\section{Objects from the Class}{
Objects of thiss class should be created using the functions
\code{\link{RoarDataset}} or \code{\link{RoarDatasetFromFiles}}, ideally
the raw \code{\link{new}} method should never be invoked by end users.
%%  ~~ describe objects here ~~ 
}
\section{Slots}{
  \describe{
    \item{\code{rightBams}:}{Object of class \code{"list"} ~~ }
    \item{\code{leftBams}:}{Object of class \code{"list"} ~~ }
    \item{\code{prePostCoords}:}{Object of class \code{"GRanges"} ~~ }
    \item{\code{postCoords}:}{Object of class \code{"GRanges"} ~~ }
    \item{\code{countsRight}:}{Object of class \code{"SummarizedExperiment"} ~~ }
    \item{\code{countsLeft}:}{Object of class \code{"SummarizedExperiment"} ~~ }
    \item{\code{pVals}:}{Object of class \code{"SummarizedExperiment"} ~~ }
    \item{\code{step}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{cores}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{exptData}:}{Object of class \code{"SimpleList"} ~~ }
    \item{\code{rowData}:}{Object of class \code{"GenomicRangesORGRangesList"} ~~ }
    \item{\code{colData}:}{Object of class \code{"DataFrame"} ~~ }
    \item{\code{assays}:}{Object of class \code{"Assays"} ~~ }
  }
}
\section{Extends}{
Class \code{"\linkS4class{SummarizedExperiment}"}, directly.
}
\section{Methods}{
  \describe{
    \item{computePvals}{\code{signature(rds = "RoarDataset")}: ... }
    \item{computeRoars}{\code{signature(rds = "RoarDataset")}: ... }
    \item{cores<-}{\code{signature(rds = "RoarDataset", value = "numeric")}: ... }
    \item{cores}{\code{signature(rds = "RoarDataset")}: ... }
    \item{countPrePost}{\code{signature(rds = "RoarDataset", stranded = "logical")}: ... }
    \item{filteringInfoResults}{\code{signature(rds = "RoarDataset")}: ... }
    \item{totalResults}{\code{signature(rds = "RoarDataset")}: ... }
    }
}
\references{
%%  ~~put references to the literature/web site here~~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%%  ~~objects to See Also as \code{\link{~~fun~~}}, ~~~
%%  ~~or \code{\linkS4class{CLASSNAME}} for links to other classes ~~~
}
\examples{
showClass("RoarDataset")
}
\keyword{classes}