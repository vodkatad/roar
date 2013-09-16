\name{totalResults}
     \alias{totalResults}
     \title{Returns a dataframe with results of the analysis for a \code{\linkS4class{RoarDataset}} object}
     \description{
         The last step of a classical Roar analyses: it returns a dataframe containing m/M values, roar 
         values and pvalues.
     }
     \usage{
      totalResults(rds)
     }
     \arguments{
       \item{rds}{
           The \code{\linkS4class{RoarDataset}} with all the analysis steps (\code{\link{countPrePost}},
            \code{\link{computeRoars}}, \code{\link{computePvals}}) performed.}
      }
      \value{
         The \code{\linkS4class{RoarDataset}} object given as rds with all the analysis steps performed.
         If one or more steps hadn't been performed they will be called automatically.
         The resulting dataframe will have the "gene_id" of the initial annotation as row names (without
         the trailing "_PRE"/"_POST") and as columns the m/M ratio for the right and left conditions,
         the roar value and the Fisher test pvalue (respectively: mM_right, mM_left, roar, pval).
         If more than one sample has been given for a condition the "pval" column will contain the
         multiplication of all the comparisons pvalue and there will be other columns containing
         the pvalues resulting from all the pairwise right vs left contrasts, with names "pvalue_X_Y" where
         X represent the position of the sample in the right list of bam files (or GappedAlignment) and
         Y the position for the left list.
         \bold{WARNING}: this method does not filter in any way the results, therefore there will be negative
         m/M values/ROAR and also NA - in these cases there aren't enough information to draw a conclusion
         about the shortening/lengthening of the gene in the given samples and thus the pvalues
         should not be kept in consideration. Furthermore there isn't any filter on the expression level
         of the genes. See \code{\link{fpkmResults}}, \code{\link{standardFilter}} and
         \code{\link{pvalueFilter}} about results filtering possibilities.
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
   df <- totalResults(rds)
    }
\keyword{totalResults}