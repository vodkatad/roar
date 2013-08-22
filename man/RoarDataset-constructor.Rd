\name{RoarDataset}
     \alias{RoarDataset}
     \title{Creates an RoarDataset object}
     \description{
       This function creates an RoarDataset object from two list of filenames
       of alignments (bam format) and a gtf with a suitable annotation of 
       alternative APA sites.
     }
     \usage{
      RoarDataset(rightBams, leftBams, gtf)
     }
     \arguments{
       \item{rightBams}{
           A list of filenames of bam alignments with data for the first condition
           to be considered.}
        \item{leftBams}{
           A list of filenames of bam alignments with data for the second condition
           to be considered.}
       \item{gtf}{
           A filename of a gtf with coordinates for the portions of transcripts
           that has to be considered pertaining to the short (or long) isoform.
           This gtf must have an attribute called "gene_id" that ends with "_PRE"
           or "_POST" to address respectively the short and the long isoform.
           A ready-to-go gtf, with coordinates derived from the PolyADB on
           the human genome (version hg19), is available TODO.
        }
      }
      \value{
         A RoarDataset object ready to be analyzed via the other methods.
      }
     \seealso{
       %\code{\link{read.HTSeqCounts}}
       %TODO
     }
     \examples{
     #rds <- RoarDataset(rightBams, leftBams, gtf)
    }
\keyword{RoarDataset}