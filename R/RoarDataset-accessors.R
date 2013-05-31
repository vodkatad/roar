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
             pre_elems <- grep("_PRE", elementMetadata(rds@prePostCoords)$gene_id)
             post_elems <- grep("_POST", elementMetadata(rds@prePostCoords)$gene_id)
             reducedGRangesPre <- rds@prePostCoords[pre_elems,]
             if (length(rds@rightBams) == 1 && length(rds@leftBams)) {
                se <- SummarizedExperiment(assays = matrix(nrow=length(rds@prePostCoords)/2, ncol=4),
                                           rowData=reducedGRangesPre, 
                                           colData=DataFrame(row.names=c("right_pre","right_post","left_pre", "left_post"))
                )
                rightSE <- summOv(rds@rightBams[[1]])
                leftSE <- summOv(rds@leftBams[[1]])  
                assay(se,1)[,"right_pre"] <- assays(rightSE)$counts[pre_elems,]
                assay(se,1)[,"right_post"] <- assays(rightSE)$counts[post_elems,]
                assay(se,1)[,"left_pre"] <- assays(leftSE)$counts[pre_elems,]
                assay(se,1)[,"left_post"] <- assays(leftSE)$counts[post_elems,]
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

setMethod("computeRoars", signature(rds="RoarDataset"),
          function(rds){
              # roar is the m/M of right condition divided by the m/M of the left one.
              # m/M = ((Lpost*Cpre)/(Lpre*Cpost))-1
              # Negative m/M are discarded.
              # roar = (m/M_right)/(m/M_left)
              # We must obtain the list of lengths.
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