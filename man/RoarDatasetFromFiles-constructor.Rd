\name{RoarDatasetFromFiles}
     \alias{RoarDatasetFromFiles}
     \title{Creates a \code{\link[=RoarDataset-class]{RoarDataset}} object}
     \description{
       This function creates an \code{\link[=RoarDataset-class]{RoarDataset}} object from two lists and a gtf with a suitable annotation of alternative APA sites.
     }
     \usage{
      RoarDatasetFromFiles(treatmentBams, controlBams, gtf)
     }
     \arguments{
            \item{treatmentBams}{
           A list of filenames of bam alignments with data for the treatment condition  (by convention it is 
            considered the ``treated'' condition: this simply means
            that the package will compute roar values (ratios of the m/M) using this condition as the numerator)
           to be considered.}
        \item{controlBams}{
           A list of filenames of bam alignments with data for the control condition
           to be considered.}
       \item{gtf}{
           A filename of a gtf with coordinates for the portions of transcripts
           that has to be considered pertaining to the short (or long) isoform.
           This gtf must have an attribute called "gene_id" that ends with "_PRE"
           or "_POST" to address respectively the short and the long isoform.
           A ready-to-go gtf, with coordinates derived from the PolyADB on
           the human genome (version hg19), is available in the "examples" package directory.
           An element in the annotation
           is considered "PRE" (i.e. common to the short and long isoform of the transcript) if its gene_id
           feature in the gtf ends with "_PRE". If
           it ends with "_POST" it is considered the portion present only in the long isoform.
           The prefix of gene_id should be an identifier for the gene and each identifier has to be
           associated with only one "_PRE" and one "_POST", leading to two genomic region associated to each   
           gene_id.
           The gtf can also contain an attribute that represents the lengths of PRE and POST portions on the transcriptome.
           If this is omitted the lengths on the genome are used instead. Note that right now every gtf entry
           (or none of them) should have it.
        }
      }
      \value{
         A \code{\link[=RoarDataset-class]{RoarDataset}} object ready to be analyzed via the other methods.
      }
     \seealso{
       \code{\link{RoarDataset}}
     }
     \examples{
     #rds <- RoarDatasetFromFiles(treatmentBams, controlBams, gtf)
    }
\keyword{RoarDatasetFromFiles}