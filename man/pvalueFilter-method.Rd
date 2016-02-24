\name{pvalueFilter}
     \alias{pvalueFilter}
     \alias{pvalueFilter, RoarDataset}
     \alias{pvalueFilter, RoarDatasetMultipleAPA}
     \title{Returns a dataframe with results of the analysis for a \code{\linkS4class{RoarDataset}} object
          or a \code{\linkS4class{RoarDatasetMultipleAPA}} object}
     \description{
         The last step of a classical Roar analyses: it returns a dataframe containing m/M values, roar 
         values, pvalues and estimates of expression (a measure recalling FPKM). Only the genes
         with an expression estimate bigger than a given cutoff will be considered. Also pvalues
         will be considered for filtering.
     }
     \usage{
      pvalueFilter(rds, fpkmCutoff, pvalCutoff)
     }
     \arguments{
       \item{rds}{
           The \code{\linkS4class{RoarDataset}} or the \code{\linkS4class{RoarDatasetMultipleAPA}} 
           with all the analysis steps (\code{\link{countPrePost}},
            \code{\link{computeRoars}}, \code{\link{computePvals}}) performed.
            If one or more steps hadn't been performed they will be called automatically
            (\bold{WARNING}: only working for \code{\linkS4class{RoarDataset}} right now).}
       \item{fpkmCutoff}{
         The cutoff that will be used to determine if a gene is expressed or not.}
       \item{pvalCutoff}{
         The cutoff that will be used to determine if a pvalue is significative or not.}
      }
      \value{
         For \code{\linkS4class{RoarDataset}}: 
         
         
         The resulting dataframe will be identical to that returned by \code{\link{standardFilter}} but
         after gene expression filtering and pvalue correction another step will be performed:
         for single samples comparisons only genes with a nominal pvalue smaller than the given cutoff
         will be considered, while for multiple samples a column (nUnderCutoff) will be added to the
         dataframe. This column will contain an integer number representing the number of comparisons
         between the samples of the two conditions that results in a nominal pvalue lower than the given    
         cutoff (pvalCutoff). For multiple samples with a paired design (i.e. if computePairedPvals was used)
         the pvalues of the requested pairings will be listed together with the combined pvalued
         obtained with the Fisher method and the filtering will be done on this pvalue.
        % \bold{WARNING}: the given cutoff will be used to determine significativeness of corrected pvalues in one situation (single sample for both conditions) and of nominal one in the other one.
        
         For \code{\linkS4class{RoarDatasetMultipleAPA}}: TODO
      } 
     \examples{
   library("GenomicAlignments")
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
   dat <- pvalueFilter(rds, 1, 0.05)
    }
\keyword{pvalueFilter}