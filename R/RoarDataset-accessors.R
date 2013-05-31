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
      preElems <- grep("_PRE", elementMetadata(rds@prePostCoords)$gene_id)
      postElems <- grep("_POST", elementMetadata(rds@prePostCoords)$gene_id)
      preCoords <- rds@prePostCoords[preElems,]
      rds@postCoords <- rds@prePostCoords[postElems,]
      if (length(rds@rightBams) == 1 && length(rds@leftBams)) {
         se <- SummarizedExperiment(assays = matrix(nrow=length(rds@prePostCoords)/2, ncol=4),
                                    rowData=preCoords, 
                                    colData=DataFrame(row.names=c("right_pre","right_post","left_pre", "left_post"))
         )
         rightSE <- summOv(rds@rightBams[[1]])
         leftSE <- summOv(rds@leftBams[[1]])  
         assay(se,1)[,"right_pre"] <- assays(rightSE)$counts[preElems,]
         assay(se,1)[,"right_post"] <- assays(rightSE)$counts[postElems,]
         assay(se,1)[,"left_pre"] <- assays(leftSE)$counts[preElems,]
         assay(se,1)[,"left_post"] <- assays(leftSE)$counts[postElems,]
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
      #roar = (m/M_right)/(m/M_left)
      # We must obtain the list of lengths from rds@postCoords and rowData(rds) (which is pre).
      preLen <- end(rowData(rds)) - start(rowData(rds))
      postLen <- end(rds@postCoords) - end(rds@postCoords)
      # Then the bam lengths to correct our lengths, ie: postLen+ReadLength-1
      if (length(rds@rightBams) == 1 && length(rds@leftBams)) {
         corrRight <- mean(qwidth(r@rightBams[[1]]))
         # qwidth(x): Returns an integer vector of length length(x) containing the length 
         # of the query *after* hard clipping (i.e. the length of the query sequence 
         # that is stored in the corresponding SAM/BAM record).
         corrLeft <- mean(qwidth(r@leftBams[[1]]))
         postLenRight <- postLen + corrRight - 1
         postLenLeft <- postLen + corrLeft - 1
         mMright <- (assay(rds,1)[,"right_pre"]*preLen)/(assay(rds,1)[,"right_post"]*postLenRight)-1
         mMleft <- (assay(rds,1)[,"left_pre"]*preLen)/(assay(rds,1)[,"left_post"]*postLenLeft)-1
         roar <- mMright / mMleft
         pVal <- rep(1, length(roar))
         assay(rds,2) <- as.matrix(data.frame(right_pre=mMright, right_post=mMleft, left_pre=roar, left_post=pVal))
      } else {      
         stop("TODO")
      }
      names(assays(rds)) <- c("counts", "stats")
      return(rds)
   }
)

setMethod("computePvals", signature(rds="RoarDataset"),
   function(rds){
      if (length(rds@rightBams) == 1 && length(rds@leftBams)) {
         if (cores(rds) == 1) {
            assay(rds,2)[,"left_post"] <- apply(assay(rds,1), 1, get_fisher)
         } else {
            stop("TODO")
         }
      } else {      
         stop("TODO")
      }
      return(rds)
   }
)

setMethod("totalResults", signature(rds="RoarDataset"),
   function(rds){
      return(data.frame(row.names=sub("^\\s+","",sub("_POST","",elementMetadata(rds@postCoords)$gene_id)), 
                        mM_right=assay(rds,2)[,"right_pre"], 
                        mM_left=assay(rds,2)[,"right_post"],
                        roar=assay(rds,2)[,"left_pre"],
                        pval=assay(rds,2)[,"left_post"]))
   }
)

# TODO ADD check on order of function calls!

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