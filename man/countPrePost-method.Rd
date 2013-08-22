\name{countPrePost}
     \alias{countPrePost}
     \title{Counts reads falling over PRE/POST portions of the given transcripts}
     \description{
       This is the first step in the Roar analyses: it counts reads overlapping with the
       PRE/POST portions defined in the given gtf/GRanges annotation. See \code{\link{RoarDataset}} for details
       on how to define these portions.
       Reads of the given bam annotation files that falls over this portion are accounted for with the following
       rules:
       1- reads falling over different gene_id PRE or POST portions are discarded
       2- reads falling on a single gene_id PRE (or POST) portion are considered, even if the fall 
       outside of the portion boundaries
       3- reads falling over a single gene coordinates but both on PRE and POST are assigned to the POST portion,
       as long as the could be derived only from the long isoform.
       If the stranded argument is set to TRUE then strandness is considered when counting reads.
     }
     \usage{
      countPrePost(rds, stranded
     }
     \arguments{
       \item{rds}{
           The \code{\linkS4class{RoarDataset}} which contains the alignments and annotation informations over which
           counts will be performed.}
       \item{stranded}{
           A logical indicating if strandness should be considered when counting reads or not.    
       }
      }
      \value{
         The \code{\linkS4class{RoarDataset}} object given as rds with the counting reads phase of the analysis done.
      }
     \examples{
     #rds <- countPrePost(rds, FALSE)
    }
\keyword{countPrePost}