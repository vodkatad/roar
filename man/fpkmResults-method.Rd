\name{fpkmResults}
     \alias{fpkmResults}
     \alias{fpkmResults, RoarDataset}
     \alias{fpkmResults, RoarDatasetMultipleAPA}
     \title{Returns a dataframe with results of the analysis for a \code{\linkS4class{RoarDataset}} object
          or a \code{\linkS4class{RoarDatasetMultipleAPA}} object}
     \description{
         The last step of a classical Roar analyses: it returns a dataframe containing m/M values, roar 
         values, pvalues and estimates of expression (a measure recalling FPKM).
     }
     \usage{
      fpkmResults(rds)
     }
     \arguments{
       \item{rds}{
           The \code{\linkS4class{RoarDataset}} or the \code{\linkS4class{RoarDatasetMultipleAPA}} 
           with all the analysis steps (\code{\link{countPrePost}},
            \code{\link{computeRoars}}, \code{\link{computePvals}}) performed.
            If one or more steps hadn't been performed they will be called automatically.}
      }
      \value{
         The resulting dataframe will be identical to that returned by \code{\link{totalResults}} but
         with two columns added: "treatmentValue" and "controlValue". These columns will contain a number
         that indicates the level of expression of the relative gene in the treatment (or control) condition.
         For \code{\linkS4class{RoarDataset}} this number derives from the counts (averages across samples when 
         applicable) obtained for the PRE portion of the gene and is similar
         to the RPKM standard measure of expression used in RNAseq experiment. Specifically
         we correct the counts on the PRE portions dividing them by portion length and total numer
         of reads aligned on all PRE portions and the multiply the results for 1000000000.
         See the vignette for more details.
         
         For \code{\linkS4class{RoarDatasetMultipleAPA}} the same procedure is applied to all the possible
         PRE choices for genes. Note that summing all the counts for every PRE portion assigned to a gene
         could lead to count some reads multiple times when summing all the PRE portions counts therefore
         this measure is not completely comparable with the one obtained with the single APA analysis.
         The length column added in this case contains the length of the PRE portions (counting only
         exonic bases).
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
   rds <- computePvals(rds)
   dat <- fpkmResults(rds)
    }
\keyword{fpkmResults}
