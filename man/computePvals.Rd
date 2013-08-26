\name{computePvals}
     \alias{computePvals}
     \title{Computes pvalues (Fisher test) for this \code{\linkS4class{RoarDataset}} object}
     \description{
       This is the third step in the Roar analyses: it applies a Fisher test comparing counts falling
       on the PRE and POST portion in the right and left for every gene. If there are multiple samples for
       a condition every combinations of comparisons between the samples lists are considered.
     }
     \usage{
      computePvals(rds)
     }
     \arguments{
       \item{rds}{
           The \code{\linkS4class{RoarDataset}} which contains the counts over PRE-POST portions
           in the two conditions to be compared via pvalues.}
      }
      \value{
         The \code{\linkS4class{RoarDataset}} object given as rds with the compute Roars phase of the           analysis done. Pvalues will be held in the RoarDataset object itself in the case of single samples, while
         in a separate slot otherwise, but end user normally should not analyze those directly but use
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
   rds <- computePvals(rds)
    }
\keyword{computeRoars}