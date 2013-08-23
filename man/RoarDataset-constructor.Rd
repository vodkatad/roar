\name{RoarDataset}
     \alias{RoarDataset}
     \title{Creates a \code{\link[=RoarDataset-class]{RoarDataset}} object}
     \description{
       This function creates an \code{\link[=RoarDataset-class]{RoarDataset}} object from two lists of 
       of \code{\link{GappedAlignments}} and a \code{\link{GRanges}} containing a suitable annotation of 
       alternative APA sites.
     }
     \usage{
      RoarDataset(rightGappedAlign, leftGappedAlign, gtfGRanges)
     }
     \arguments{
       \item{rightGappedAlign}{
           A list of \code{\link{GappedAlignments}} representing alignment of samples for the first condition
           to be considered.}
        \item{leftGappedAlign}{
           A list of \code{\link{GappedAlignments}} representing alignment of samples for the second condition
           to be considered.}
       \item{gtfGRanges}{
           A \code{\link{GRanges}} object with coordinates for the portions of transcripts
           that has to be considered pertaining to the short (or long) isoform.
           This \code{\link{GRanges}} object must have a character elementMetadata called "gene_id" 
           that ends with "_PRE" or "_POST" to address respectively the short and the long isoform.
           An element in the annotation
           is considered "PRE" (i.e. common to the short and long isoform of the transcript) if its gene_id
           ends with "_PRE". If it ends with "_POST" it is considered the portion present only in the long
           isoform.
           The prefix of gene_id should be an identifier for the gene and each identifier has to be
           associated with only one "_PRE" and one "_POST", leading to two genomic region associated to each   
           gene_id.
        }
      }
      \value{
         A \code{\link[=RoarDataset-class]{RoarDataset}} object ready to be analyzed via the other methods.
         %\code{"\linkS4class{RoarDataset}"}
      }
     \seealso{
       %\code{\link{read.HTSeqCounts}}
       %TODO
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
    }
\keyword{RoarDataset}