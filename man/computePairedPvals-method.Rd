\name{computePairedPvals}
     \alias{computePairedPvals}
     \title{Computes pvalues (Fisher test) for this \code{\linkS4class{RoarDataset}} object}
     \description{
       This is the third step in the Roar analyses: it applies a Fisher test comparing counts falling
       on the PRE and POST portion in the treatment and control conditions for every gene. 
       The paired method should be used when the experimental setup offers multiple paired samples
       for the two conditions: that is foreach sample of the control condition there is a naturally
       paired one for the treatment (i.e. cells derived from the same plate divided in two groups
       and treated or not). For example in the below code sample treatment sample n.1 (rd1) is paired with control n.2 (rd4)
       and rd2 with rd3. The pvalue resulting from Fisher test applied on the different samples pairings will
       be combined with the Fisher method, therefore the pairs of samples should be independent between each other.
      }
     \usage{
      computePairedPvals(rds, treatmentSamples, controlSamples)
     }
     \arguments{
       \item{rds}{
           The \code{\linkS4class{RoarDataset}} which contains the counts over PRE-POST portions
           in the two conditions to be compared via pvalues.}
      \item{treatmentSamples}{Numbers that represent the indexes of the treatmentBams/GappedAlign parameter given to
         the RoarDataset costructor and the order in which they are paired with control samples.}
      \item{controlSamples}{Numbers that represent the indexes of the controlBams/GappedAlign parameter given to
         the RoarDataset costructor and the order in which they are paired with treatment samples.}
      }      
      \value{
         The \code{\linkS4class{RoarDataset}} object given as rds with the compute pvalues phase of the analysis done. Pvalues will be held in the RoarDataset object itself in the case of single samples, while
         in a separate slot otherwise, but end user normally should not analyze those directly but use
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
   rd4 <- GAlignments("a", seqnames = Rle("chr2"), pos = as.integer(3400), cigar = "300M", strand = strand("+"))
   rds <- RoarDataset(list(rd1,rd2), list(rd3, rd4), features)
   rds <- countPrePost(rds, FALSE)
   rds <- computeRoars(rds)
   rds <- computePairedPvals(rds, c(1,2), c(2,1))
    }
\keyword{computePairedPvals}
