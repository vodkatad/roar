\name{computeRoars}
     \alias{computeRoars}
     \alias{computeRoars, RoarDataset, numeric, numeric}
     \alias{computeRoars, RoarDatasetMultipleAPA, numeric, numeric}
     \title{Computes m/M and roar values}
     \description{
       This is the second step in the Roar analyses: it computes the ratio of prevalence
       of the short and long isoforms for every gene in the treatment and control condition (m/M) and their
       ratio, roar, that indicates if there is a relative shortening-lengthening in a condition
       over the other one. A roar > 1 for a given gene means that in the treatment condition
       that gene has an higher ratio of short vs long isoforms with respect to the control condition (and
       the opposite for roar < 1). Negative or NA m/M or roar occurs in not definite situations,
       such as counts equal to zero for PRE or POST portions.
       If for one of the conditions there are more than one samples then calculations are performed
       on average counts.
     }
     \usage{
         %\S4method{computeRoars}{RoarDataset,numeric,numeric}(rds, qwidthTreatment=NA, qwidthControl=NA)
         %\S4method{computeRoars}{RoarDataset,numeric,numeric}(rds, qwidthTreatment, qwidthControl)
         computeRoars(rds, qwidthTreatment=NA, qwidthControl=NA)
         computeRoars(rds, qwidthTreatment, qwidthControl)
     }
     \arguments{
       \item{rds}{
           The \code{\linkS4class{RoarDataset}} or the \code{\linkS4class{RoarDatasetMultipleAPA}} which contains the   
           counts over PRE-POST portions in the two conditions to be compared via roar.}
         \item{qwidthTreatment}{ The mean length of the reads in the treatment bam files
            - used internally for the interaction between
            \code{\linkS4class{RoarDataset}} and \code{\linkS4class{RoarDatasetMultipleAPA}} objects.
            The default (NA) calculates this value from the bam and should not be changed.}
         \item{qwidthControl}{The mean length of the reads in the control bam files 
            - used internally for the interaction between
            \code{\linkS4class{RoarDataset}} and \code{\linkS4class{RoarDatasetMultipleAPA}} objects. 
            The default (NA) calculates this value from the bam and should not be changed.}
      }
      \value{
         The \code{\linkS4class{RoarDataset}} or the \code{\linkS4class{RoarDatasetMultipleAPA}} object given as rds             with the computeRoars phase of the  
         analysis done. m/M and roars will be held in the RoarDataset object itself in the case of single samples,
         while
         in two slots otherwise, but end user normally should not analyze those directly but use
         \code{\link{totalResults}} or \code{\link{fpkmResults}} at the end of the analysis. 
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
   rds <- countPrePost(rds, FALSE)
   rds <- computeRoars(rds)
    }
\keyword{computeRoars}
