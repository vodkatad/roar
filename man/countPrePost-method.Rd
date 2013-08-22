\name{countPrePost}
     \alias{countPrePost}
     \title{Counts reads falling over PRE/POST portions of the given transcripts}
     \description{
       This is the first step in the Roar analyses: it counts reads overlapping with the
       PRE/POST portions defined in the given gtf/GRanges annotation. An element in the annotation
       is considered "PRE" (i.e. common to the short and long isoform of the transcript) if its gene_id
       feature in the gtf (loaded in the elementMetadata of the GRanges object) ends with "_PRE". If
       it ends with "_POST" it is considered the portion present only in the long isoform.
       The prefix of gene_id should be an identifier for the gene and each identifier has to be
       associated with only one "_PRE" and one "_POST", leading to two genomic region associated to each gene_id.
       Reads of the given bam annotation files that falls over this portion are accounted for with the following
       rules:
       1- reads falling over different gene_id PRE or POST portions are discarded
       2- reads falling on a single gene_id PRE (or POST) portion are considered, even if the fall 
       outside of the portion boundaries
       3- reads falling over a single gene coordinates but both on PRE and POST are assigned to the POST portion,
       as long as the could be derived only from the long isoform.
       If the stranded argument is set to TRUE then strandness is considered when counting reads.
     }
     \usage{
      countPrePost(rds, stranded
     }
     \arguments{
       \item{rds}{
           The RoarDataset which contains the alignments and annotation informations over which
           counts will be performed.}
       \item{stranded}{
           A logical indicating if strandness should be considered when counting reads or not.    
       }
      }
      \value{
         The RoadDataset object given as rds with the counting reads phase of the analysis done.
      }
     \examples{
     #rds <- countPrePost(rds, FALSE)
    }
\keyword{countPrePost}