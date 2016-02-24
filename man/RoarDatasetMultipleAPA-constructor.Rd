\name{RoarDatasetMultipleAPA}
     \alias{RoarDatasetMultipleAPA}
     \title{Creates a \code{\link[=RoarDatasetMultipleAPA-class]{RoarDatasetMultipleAPA}} object}
     \description{
       This function creates an \code{\link[=RoarDatasetMultipleAPA-class]{RoarDatasetMultipleAPA}} object from two lists of 
       of \code{\link{GAlignments}} and a \code{\link{GRanges}} containing a suitable annotation of 
       alternative APA sites and gene exon structure.
      A MultipleAPA analysis computes several roar values and p-values for each gene: one 
      for every possible combination of APA-canonical end of a gene (i.e. the
      end of its last exon). This is more efficient than performing several different ``standard'' roar analyses
      choosing the PRE and POST portions corresponding to different APAs because reads overlaps are computed only once.
     }
     \usage{
      RoarDatasetMultipleAPA(treatmentBamsGenomicAlignments, controlBamsGenomicAlignments, gtfGRanges)
     }
     \arguments{
       \item{treatmentBamsGenomicAlignments}{
           A list of \code{\link{GAlignments}} representing alignment of samples for the treatment 
           condition (by convention it is considered the ``treated'' condition: this simply means
            that the package will compute roar values (ratios of the m/M) using this condition as the numerator)
           to be considered.}
        \item{controlBamsGenomicAlignments}{
           A list of \code{\link{GAlignments}} representing alignment of samples for the control condition
           to be considered.}
       \item{gtfGRanges}{
           A \code{\link{GRanges}} containing a suitable annotation of alternative APA sites and gene exonic structure.
           Minimal requirements are: metadata columns called "gene", "apa" and "type."
           APA should be single bases falling over one of the given genes and need to have the
           metadata column "type" equal to "apa" and the "apa"  column
           composed of unambiguous id and the corresponding gene id pasted together with an underscore. The "gene"
           metadata columns for these entries should not be initialized.
           All the studied gene exons need to be reported, in this case the metadata column "gene" should contain 
           the gene id (the same
           one reported for each gene APAs) while "type" should be set to "gene" and "apa" to NA.
           All apa entries assigned to a gene should have coordinates that falls inside it (otherwise XXX TODO CHECK) and             every gene that
           appears should contain at least one APA.
           
        }
      }
      \value{
         A \code{\link[=RoarDatasetMultipleAPA-class]{RoarDatasetMultipleAPA}} object ready to be analyzed via the other methods.
      }
     \seealso{
       \code{\link{RoarDatasetMultipleAPAFromFiles}}
     }
     \examples{   
      library(GenomicAlignments)
      gene <- c("A", "B", NA, NA)
      type <- c("gene","gene","apa", "apa")
      apa <- c(NA, NA, "apa1_A", "apa2_B")
      features <- GRanges(
         seqnames = Rle(c("chr1", "chr2", "chr1", "chr2")),
         strand = strand(rep("+", length(gene))),
         ranges = IRanges(
            start=c(1000, 2000, 1300, 2050),
            width=c(500, 900, 1, 1)),
         DataFrame(gene, apa, type)
      )
      rd1 <- GAlignments("a", seqnames = Rle("chr1"), pos = as.integer(1000), cigar = "300M", strand = strand("+"))
      rds <- RoarDatasetMultipleAPA(list(c(rd1,rd1)), list(c(rd1,rd1)), features)
    }
\keyword{RoarDataset}
