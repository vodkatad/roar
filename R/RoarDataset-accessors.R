# Methods for the RoarDataset class.
# R/AllGenerics.R contains the signature of the not overriding methods.

RoarDataset <- function(rightBams, leftBams, gtf) {
   genome = "hg19" # FIXME
   # The format will be assumed using the file extension. Will work everytime?
   # Do we force hg19?
   gtfGRanges<- import(gtf, genome=genome, asRangedData=F)
   rightBamsGenomicAlignments <- lapply(rightBams, readGappedAlignments)
   leftBamsGenomicAlignments <- lapply(leftBams, readGappedAlignments)
   #rightBamsGenomicAlignments <- BamFileList(rightBams)
   #leftBamsGenomicAlignments <- BamFileList(leftBams)
   new("RoarDataset", rightBams=rightBamsGenomicAlignments, leftBams=leftBamsGenomicAlignments, prePostCoords=gtfGRanges, cores=1)
}

# Could have used setMethod("initialize", "xx",) but in this way should have had a gtf filename slot.

setMethod("countPrePost", signature(rds="RoarDataset"),
   function(rds){
      #if (!is(rds, "RoarDataset")) {
      #   stop("countPrePost could be applied only to RoarDataset objects")
      #} # Why is this needed? Is it needed?
        # To me it does not seem nedeed.
      summOv <- function(x) {
         summarizeOverlaps(features=rds@prePostCoords, reads=x, mode='Union', ignore.strand=T, mc.cores=rds@coresproa)
      } 
      # Now we need to keep means and totals of counts over PRE/POST for the two lists.
      # In the simpler case with a single alignment for both conditions we just keep the counts.
      if (length(rds@rightBams) == 1 && length(rds@leftBams)) {
         se <- SummarizedExperiment(assays = matrix(nrow=length(rds@prePostCoords), ncol=2),
                              rowData=rds@prePostCoords, 
                              colData=DataFrame(row.names=c("right","left"))
                              )
         assay(se,1)[,"right"] <- assays(summOv(rds@rightBams[[1]]))$counts
         assay(se,1)[,"left"] <- assays(summOv(rds@leftBams[[1]]))$counts
         rowData(rds) <- rowData(se)
         colData(rds) <- colData(se)
         assays(rds) <- assays(se)
         names(assays(rds)) <- "counts"
      } else {
         stop("TODO")
         # Ideally here will wet counts for all right and left bams, compute means and totals
         # and obtain a SE with 2 assays, called means and totals.
         #summarizedRight <- lapply(rds@rightBamsGenomicAlignments, summOv)
         #summarizedLeft <- lapply(rds@leftBamsGenomicAlignments, summOv)
      }

      return(rds)
   }
)

# Simple getters and setters. Arf Arf!
setMethod("cores",  signature(rds="RoarDataset"),
   function(rds) {
      return(rds@cores)
   }
)

setReplaceMethod("cores",  signature(rds="RoarDataset", value="numeric"),
   function(rds, value) {
      if (value < 1)
         stop("You can specify only a positive number of cores")
      rds@cores <- value
      return(rds)
   }
)