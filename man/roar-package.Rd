\name{roar-package}
\alias{roar-package}
\alias{roar}
\docType{package}
\title{
Identify differential APA usage from RNA-seq alignments
}
\description{
Identify preferential usage of APA sites, 
comparing two biological conditions, starting from known alternative sites and
alignments obtained from standard RNA-seq experiments.
}
\details{
The code{\linkS4class{RoarDataset}} class exposes methods to perform the whole analysis, in order to 
identify genes with preferential expression of long/short isoforms in a condition with respect to 
another one. The needed input data are alignments deriving from RNA-seq experiments of the
two conditions and a set of coordinates of APA sites for genes with an alternative APA site
proximal to the one used ``normally''.
}
\author{
Elena Grassi <grassi.e@gmail.com>
}
\keyword{ package }

