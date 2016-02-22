\name{countPrePost} 
     \alias{countPrePost}
     \alias{countPrePost,RoarDataset,logical}
     \alias{countPrePost,RoarDatasetMultipleAPA}
     \title{Counts reads falling over PRE/POST portions of the given transcripts}
     \description{
       This is the first step in the Roar analyses: it counts reads overlapping with the
       PRE/POST portions defined in the given gtf/GRanges annotation. See \code{\link{RoarDataset}} for details
       on how to define these portions.
       Reads of the given bam annotation files that falls over this portion are accounted for with the following
       rules:
       
       1- reads that align on only one of the given features are assigned to that feature, even if the
       overlap is not complete 
       2- reads that align on both a PRE and a POST feature of the same gene (spanning reads) are 
       assigned to the POST one, considering that they have clearly been obtained from the longest isoform
       
       If the stranded argument is set to TRUE then strandness is considered when counting reads.
       When rds is a \code{\link{RoarDatasetMultipleAPA}} counts are obtained on more than two portions for
       each transcript in order to be able to efficiently evaluate multiple APA sites. 
       The options stranded=TRUE is still not implemented for \code{\linkS4class{RoarDatasetMultipleAPA}}.
       
     }
     \usage{
      %\S4method{countPrePost}{RoarDataset,logical}(rds, stranded=FALSE)
      %\S4method{countPrePost}{RoarDatasetMultipleAPA,logical}(rds, stranded=FALSE)
      countPrePost(rds, stranded=FALSE)
     }
     \arguments{
       \item{rds}{
           The \code{\linkS4class{RoarDataset}} or the \code{\linkS4class{RoarDatasetMultipleAPA}} which contains the             alignments and annotation informations over which
           counts will be performed.}
       \item{stranded}{
           A logical indicating if strandness should be considered when counting reads or not.
           Default=FALSE. \bold{WARNING}: not implemented (ignored) when using \code{\linkS4class{RoarDatasetMultipleAPA}}.
       }
      }
      \value{
         The \code{\linkS4class{RoarDataset}} object given as rds with the counting reads phase of the analysis done. Counts will be held in the RoarDataset object itself in the case of single samples, while
      }
     \examples{
   library(GenomicAlignments)
   gene_id <- c("A_PRE", "A_POST", "B_PRE", "B_POST")
   features <- GRanges(
      seqnames = Rle(c("chr1", "chr1", "chr2", "chr2")),
      strand = strand(rep("+", length(gene_id))),
      ranges = IRanges(
         start=c(1000, 2000, 3000, 3600),
         width=c(1000, 900, 600, 300)),
      DataFrame(gene_id)
   )
   rd1 <- GAlignments("a", seqnames = Rle("chr1"), pos = as.integer(1000), cigar = "300M", strand = strand("+"))
   rd2 <- GAlignments("a", seqnames = Rle("chr1"), pos = as.integer(2000), cigar = "300M", strand = strand("+"))
   rd3 <- GAlignments("a", seqnames = Rle("chr2"), pos = as.integer(3000), cigar = "300M", strand = strand("+"))
   rds <- RoarDataset(list(c(rd1,rd2)), list(rd3), features)
   rds <- countPrePost(rds)
    }
\keyword{countPrePost}
