\name{standardFilter}
     \alias{standardFilter}
     \title{Returns a dataframe with results of the analysis for a \code{\linkS4class{RoarDataset}} object}
     \description{
         The last step of a classical Roar analyses: it returns a dataframe containing m/M values, roar 
         values, pvalues and estimates of expression (a measure recalling FPKM). Only the genes
         with an expression estimate bigger than a given cutoff will be considered and the pvalue
         will be corrected for multiple comparisons (with the Bonferroni method).
     }
     \usage{
      totalResults(rds)
     }
     \arguments{
       \item{rds}{
           The \code{\linkS4class{RoarDataset}} with all the analysis steps (\code{\link{countPrePost}},
            \code{\link{computeRoars}}, \code{\link{computePvals}}) performed.
            If one or more steps hadn't been performed they will be called automatically.}
      \item{fpkmCutoff}{
         The cutoff that will be used to determine if a gene is expressed or not.}
      }
      \value{
         The resulting dataframe will be identical to that returned by \code{link{filteringInfoResults}} but
         it will contains rows relative only with genes with an expression estimate (right or leftValue)
         bigger than the given fpkmCutoff in both the conditions.
         The pvalues (for multiple samples their multiplication, i.e.
         the value in the "pval" column) will then be corrected for multiple comparisons (Bonferroni method)
         only considering those genes, that are considered as "expressed" in both the biological conditions
         in analysis.
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
   df <- standardFilter(rds)
    }
\keyword{standardFilter}