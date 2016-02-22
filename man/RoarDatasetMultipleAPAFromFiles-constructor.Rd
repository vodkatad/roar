\name{RoarDatasetMultipleAPAFromFiles}
     \alias{RoarDatasetMultipleAPAFromFiles}
     \title{Creates a \code{\link[=RoarDatasetMultipleAPA-class]{RoarDatasetMultipleAPA}} object}
     \description{
       This function creates an \code{\link[=RoarDatasetMultipleAPA-class]{RoarDatasetMultipleAPA}} object from two lists and a gtf with a suitable annotation of 
       alternative APA sites and exonic structures of genes.
                  A MultipleAPA analysis computes several roar values and p-values for each gene: one 
           for every possible combination of APA-canonical end of a gene (i.e. the
           end of its last exon). This is more efficient than performing several different ``standard'' roar analyses
           choosing the PRE and POST portions corresponding to different APAsbecause reads overlaps are computed only once.
     }
     \usage{
      RoarDatasetMultipleAPAFromFiles(treatmentBams, controlBams, gtf)
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
           A filename of a gtf with coordinates for alternative APA sites and gene exonic structure.
           This gtf must have three attributes called "gene", "apa" and "type" to distinguish different features.
           APA should be single bases falling over one of the given genes and need to have the attribute "type" equal 
           to "apa" and the "apa" attribute
           composed of unambiguous id and the corresponding gene id pasted together with an underscore. The "gene"
           attributes for these entries should not be initialized.
           All the studied gene exons need to be reported, in this case the attribute "gene" should contain the gene id (the same
           one reported for each gene APAs) while "type" should be set to "gene" and "apa" to NA.
           All apa entries assigned to a gene should have coordinates that falls inside it (otherwise XXX TODO CHECK) and every gene that
           appears should contain at least one APA.
           A ready-to-go gtf, with coordinates derived from the PolyADB on
           the human genome (version hg19), is available in the "examples" package directory.
        }
      }
      \value{
         A \code{\link[=RoarDatasetMultipleAPA-class]{RoarDatasetMultipleAPA}} object ready to be analyzed via the other methods.
      }
     \seealso{
       \code{\link{RoarDatasetMultipleAPA}}
     }
     \examples{
     #rds <- RoarDatasetMultipleAPAFromFiles(treatmentBams, controlBams, gtf)
    }
\keyword{RoarDatasetFromFiles}