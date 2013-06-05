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
   new("RoarDataset", rightBams=rightBamsGenomicAlignments, leftBams=leftBamsGenomicAlignments, 
       prePostCoords=gtfGRanges, step = 0, cores=1)
}

# Could have used setMethod("initialize", "xx",) but in this way should have had a gtf filename slot.

setMethod("countPrePost", signature(rds="RoarDataset"),
   function(rds){
      #if (!is(rds, "RoarDataset")) {
      #   stop("countPrePost could be applied only to RoarDataset objects")
      #} # Why is this needed? Is it needed?
      # To me it does not seem nedeed.
      goOn <- checkStep(rds, 0)
      if (!goOn[[1]]) {
         return(rds)
      }
      rds <- goOn[[2]]
      summOv <- function(x) {
         summarizeOverlaps(features=rds@prePostCoords, reads=x, ignore.strand=T, mc.cores=rds@cores)
      }
      
      summOvPost <- function(x) {
         summarizeOverlaps(features=rds@postCoords, reads=x, ignore.strand=T, mc.cores=rds@cores)
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
         rightSEpost <- summOvPost(rds@rightBams[[1]])
         leftSEpost <- summOvPost(rds@leftBams[[1]])
         assay(se,1)[,"right_pre"] <- assays(rightSE)$counts[preElems,]
         assay(se,1)[,"right_post"] <- assays(rightSEpost)$counts 
         assay(se,1)[,"left_pre"] <- assays(leftSE)$counts[preElems,]
         assay(se,1)[,"left_post"] <- assays(leftSEpost)$counts # is the order conserved?
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
      rds@step <- 1;
      return(rds)
   }
)

setMethod("computeRoars", signature(rds="RoarDataset"),
   function(rds){
      #if (rds@step < 1) {
      #   rds <- countPrePost(rds)         
      #} else {
      #   # warning
      #}
      goOn <- checkStep(rds, 1)
      if (!goOn[[1]]) {
         return(rds)
      }
      rds <- goOn[[2]]
      # roar is the m/M of right condition divided by the m/M of the left one.
      # m/M = ((Lpost*Cpre)/(Lpre*Cpost))-1
      # Negative m/M are discarded.
      #roar = (m/M_right)/(m/M_left)
      # We must obtain the list of lengths from rds@postCoords and rowData(rds) (which is pre).
      preLen <- end(rowData(rds)) - start(rowData(rds))
      postLen <- end(rds@postCoords) - start(rds@postCoords)
      # Then the bam lengths to correct our lengths, ie: postLen+ReadLength-1
      if (length(rds@rightBams) == 1 && length(rds@leftBams)) {
         corrRight <- mean(qwidth(rds@rightBams[[1]]))
         # qwidth(x): Returns an integer vector of length length(x) containing the length 
         # of the query *after* hard clipping (i.e. the length of the query sequence 
         # that is stored in the corresponding SAM/BAM record).
         corrLeft <- mean(qwidth(rds@leftBams[[1]]))
         postLenRight <- postLen + corrRight - 1
         postLenLeft <- postLen + corrLeft - 1
         mMright <- (assay(rds,1)[,"right_pre"]*postLenRight)/(assay(rds,1)[,"right_post"]*preLen)-1
         mMleft <- (assay(rds,1)[,"left_pre"]*postLenLeft)/(assay(rds,1)[,"left_post"]*preLen)-1
         roar <- mMright / mMleft
         pVal <- rep(1, length(roar))
         assay(rds,2) <- as.matrix(data.frame(right_pre=mMright, right_post=mMleft, left_pre=roar, left_post=pVal))
      } else {      
         stop("TODO")
      }
      names(assays(rds)) <- c("counts", "stats")
      rds@step <- 2
      return(rds)
   }
)

setMethod("computePvals", signature(rds="RoarDataset"),
   function(rds){
#       if (rds@step < 2) {
#          if (rds@step < 1) {
#             rds <- countPrePost(rds)         
#          }
#          rds <- computeRoars(rds)
#       } else {
#          # warning
#       }
      goOn <- checkStep(rds, 2)
      if (!goOn[[1]]) {
         return(rds)
      }
      rds <- goOn[[2]]
      if (length(rds@rightBams) == 1 && length(rds@leftBams)) {
         if (cores(rds) == 1) {
            assay(rds,2)[,"left_post"] <- apply(assay(rds,1), 1, get_fisher)
         } else {
            stop("TODO")
         }
      } else {      
         stop("TODO")
      }
      rds@step <- 3
      return(rds)
   }
)

setMethod("totalResults", signature(rds="RoarDataset"),
   function(rds){
#       if (rds@step < 3) {
#          if (rds@step < 1) {
#             rds <- countPrePost(rds)         
#          }
#          if (rds@step < 2) {
#             rds <- computeRoars(rds)         
#          }
#          rds <- computePvals(rds)
#       } else {
#          # warning
#       }
      goOn <- checkStep(rds, 3)
      rds <- goOn[[2]]
      rds@step <- 4
      return(data.frame(row.names=sub("^\\s+","",sub("_POST","",elementMetadata(rds@postCoords)$gene_id)), 
                        mM_right=assay(rds,2)[,"right_pre"], 
                        mM_left=assay(rds,2)[,"right_post"],
                        roar=assay(rds,2)[,"left_pre"],
                        pval=assay(rds,2)[,"left_post"]))
   }
)

setMethod("filteredResults", signature(rds="RoarDataset"),
   function(rds){
      df <- totalResults(rds)
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