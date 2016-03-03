\name{standardFilter}
     \alias{standardFilter}
     \alias{standardFilter, RoarDataset}
     \alias{standardFilter, RoarDatasetMultipleAPA}
     \title{Returns a dataframe with results of the analysis for a \code{\linkS4class{RoarDataset}} object
          or a \code{\linkS4class{RoarDatasetMultipleAPA}} object}
     \description{
         The last step of a classical Roar analyses: it returns a dataframe containing m/M values, roar 
         values, pvalues and estimates of expression (a measure recalling FPKM). Only the genes
         with an expression estimate bigger than a given cutoff will be considered.
         %and the pvalue
         %will be corrected for multiple comparisons (with the Bonferroni method).
     }
     \usage{
      standardFilter(rds, fpkmCutoff)
     }
     \arguments{
       \item{rds}{
           The \code{\linkS4class{RoarDataset}} or the \code{\linkS4class{RoarDatasetMultipleAPA}} 
           with all the analysis steps (\code{\link{countPrePost}},
            \code{\link{computeRoars}}, \code{\link{computePvals}}) performed.
            If one or more steps hadn't been performed they will be called automatically.}
      \item{fpkmCutoff}{
         The cutoff that will be used to determine if a gene is expressed or not.}
      }
      \value{
         For \code{\linkS4class{RoarDataset}} and \code{\linkS4class{RoarDatasetMultipleAPA}}: 
         
         The resulting dataframe will be identical to that returned by \code{\link{fpkmResults}} but
         it will contains rows relative only with genes with an expression estimate (treatment or controlValue)
         bigger than the given fpkmCutoff in both the conditions and with sensitive m/M and roar values (it
         removes negative or NA
         m/M values/roar - these values arise when there aren't enough information to draw a conclusion
         about the shortening/lengthening of the gene).
         %The pvalues (for multiple samples their multiplication, i.e.
         %the value in the "pval" column) will then be corrected for multiple comparisons (Bonferroni method)
         %only considering those genes, that are considered as "expressed" in both the biological conditions
         %in analysis.
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
   dat <- standardFilter(rds, 1)
    }
\keyword{standardFilter}
