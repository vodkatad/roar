\name{computeRoars}
     \alias{computeRoars}
     \title{Compute m/M and roar values for this \code{\linkS4class{RoarDataset}} object}
     \description{
       This is the second step in the Roar analyses: it compute the ratio of prevalence
       of the short and long isoforms for every gene in the right and left condition (m/M) and their
       ratio, roar, that indicates if there is a relative shortening-lengthening in a condition
       over the other one. A roar > 1 for a given gene means that in the right condition
       that gene has an higher ratio of short vs long isoforms with respect to the left condition (and
       the opposite for roar < 1). Negative or NA m/M or roar occurs in not-valuable situations,
       such as counts equal to zero for PRE or POST portions.
       If for one of the conditions there are more than one samples then calculations are performed
       on the average counts.
     }
     \usage{
      computeRoars(rds)
     }
     \arguments{
       \item{rds}{
           The \code{\linkS4class{RoarDataset}} which contains the counts over PRE-POST portions
           in the two conditions to be compared via roar.}
      }
      \value{
         The \code{\linkS4class{RoarDataset}} object given as rds with the compute Roars phase of the           analysis done. m/M and roars will be held in the RoarDataset object itself in the case of single samples, while
         in two slots otherwise, but end user normally should not analyze those directly but use
         \code{\link{totalResults}} or \code{\link{filteringInfoResults}} at the end of the analysis. 
      }
     \examples{
        gene_id <- c("A_PRE", "A_POST", "B_PRE", "B_POST")
   features <- GRanges(
      seqnames = Rle(c("chr1", "chr1", "chr2", "chr2")),
      strand = strand(rep("+", length(gene_id))),
      ranges = IRanges(
         start=c(1000, 2000, 3000, 3600),
         width=c(1000, 900, 600, 300)),
      DataFrame(gene_id)
   )
   rd1 <- GappedAlignments("a", seqnames = Rle("chr1"), pos = as.integer(1000), cigar = "300M", strand = strand("+"))
   rd2 <- GappedAlignments("a", seqnames = Rle("chr1"), pos = as.integer(2000), cigar = "300M", strand = strand("+"))
   rd3 <- GappedAlignments("a", seqnames = Rle("chr2"), pos = as.integer(3000), cigar = "300M", strand = strand("+"))
   rds <- RoarDataset(list(c(rd1,rd2)), list(rd3), features)
   rds <- countPrePost(rds, FALSE)
   rds <- computeRoars(rds)
    }
\keyword{computeRoars}