\name{RoarDataset}
     \alias{RoarDataset}
     \title{Creates a \code{\link[=RoarDataset-class]{RoarDataset}} object}
     \description{
       This function creates an \code{\link[=RoarDataset-class]{RoarDataset}} object from two lists of 
       of \code{\link{GAlignments}} and a \code{\link{GRanges}} containing a suitable annotation of 
       alternative APA sites.
     }
     \usage{
      RoarDataset(treatmentGappedAlign, controlGappedAlign, gtfGRanges)
     }
     \arguments{
       \item{treatmentGappedAlign}{
           A list of \code{\link{GAlignments}} representing alignment of samples for the treatment 
           condition (by convention it is considered the ``treated'' condition: this simply means
            that the package will compute roar values (ratios of the m/M) using this condition as the numerator)
           to be considered.}
        \item{controlGappedAlign}{
           A list of \code{\link{GAlignments}} representing alignment of samples for the control condition
           to be considered.}
       \item{gtfGRanges}{
           A \code{\link{GRanges}} object with coordinates for the portions of transcripts
           that has to be considered pertaining to the short (or long) isoform.
           This \code{\link{GRanges}} object must have a character metadata column called "gene_id" 
           that ends with "_PRE" or "_POST" to address respectively the short and the long isoform.
           An element in the annotation
           is considered "PRE" (i.e. common to the short and long isoform of the transcript) if its gene_id
           ends with "_PRE". If it ends with "_POST" it is considered the portion present only in the long
           isoform.
           The prefix of gene_id should be a unique identifier for the gene and each identifier has to be
           associated with only one "_PRE" and one "_POST", leading to two genomic region associated to each   
           gene_id.
           The \code{\link{GRanges}} object can also contain a numeric 
           metadata column that represents the lengths of PRE and POST portions on the transcriptome.
           If this is omitted the lengths on the genome are used instead. Note that right now every gtf entry 
           (or none of them) should have it.
        }
      }
      \value{
         A \code{\link[=RoarDataset-class]{RoarDataset}} object ready to be analyzed via the other methods.
         %\code{"\linkS4class{RoarDataset}"}
      }
     \seealso{
       \code{\link{RoarDatasetFromFiles}}
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
    }
\keyword{RoarDataset}
