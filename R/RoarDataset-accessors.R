# Methods for the RoarDataset class.
# R/AllGenerics.R contains the signature of the not overriding methods.

RoarDataset <- function(rightBams, leftBams, gtf) {
   # The format will be assumed using the file extension. Will work everytime?
   # Do we need to force a genome(eg. hg19)? It doesn't seem so.
   gtfGRanges<- import(gtf, asRangedData=F)
   rightBamsGenomicAlignments <- lapply(rightBams, readGappedAlignments)
   leftBamsGenomicAlignments <- lapply(leftBams, readGappedAlignments)
   #rightBamsGenomicAlignments <- BamFileList(rightBams)
   #leftBamsGenomicAlignments <- BamFileList(leftBams)
   new("RoarDataset", rightBams=rightBamsGenomicAlignments, leftBams=leftBamsGenomicAlignments, 
       prePostCoords=gtfGRanges, step = 0, cores=1)
}

# Could have used setMethod("initialize", "xx",) but in this way should have had a gtf filename slot.

setMethod("countPrePost", signature(rds="RoarDataset", stranded="logical"),
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
         summarizeOverlaps(features=rds@prePostCoords, reads=x, ignore.strand=stranded, mc.cores=rds@cores)
      }
      summOvPost <- function(x) {
         summarizeOverlaps(features=rds@postCoords, reads=x, ignore.strand=stranded, mc.cores=rds@cores)
      } 
      
      # Now we need to keep means and totals of counts over PRE/POST for the two lists.
      # In the simpler case with a single alignment for both conditions we just keep the counts.
      preElems <- grep("_PRE", elementMetadata(rds@prePostCoords)$gene_id)
      postElems <- grep("_POST", elementMetadata(rds@prePostCoords)$gene_id)
      preCoords <- rds@prePostCoords[preElems,]
      rds@postCoords <- rds@prePostCoords[postElems,]
      se <- SummarizedExperiment(assays = matrix(nrow=length(rds@prePostCoords)/2, ncol=4),
                                 rowData=preCoords, 
                                 colData=DataFrame(row.names=c("right_pre","right_post","left_pre", "left_post"))
      )
      if (length(rds@rightBams) == 1 && length(rds@leftBams) == 1) {
         # We obtain counts for both conditions on PRE and POST coords.
         rightSE <- summOv(rds@rightBams[[1]])
         leftSE <- summOv(rds@leftBams[[1]])    
         # We obtain counts for both conditions only on POST coords (notwithstanding overlap with PRE) -
         # this is needed to implement the prefer-POST policy where a read overlapping PRE and POST
         # is counted on the POST portion (it has to be derived from a long isoform).
         # Note that in the first case we need to use the whole gtf with PRE and POST to avoid
         # assigning to PRE the overlapping reads, even if the counts on POST gotten afterwards
         # will be discarded.
         rightSEpost <- summOvPost(rds@rightBams[[1]])
         leftSEpost <- summOvPost(rds@leftBams[[1]])
         assay(se,1)[,"right_pre"] <- assays(rightSE)$counts[preElems,]
         assay(se,1)[,"right_post"] <- assays(rightSEpost)$counts 
         assay(se,1)[,"left_pre"] <- assays(leftSE)$counts[preElems,]
         assay(se,1)[,"left_post"] <- assays(leftSEpost)$counts # Is the order conserved?
         rowData(rds) <- rowData(se)
         colData(rds) <- colData(se)
         assays(rds) <- assays(se)
         names(assays(rds)) <- "counts"
      } else {
         # Ideally here will wet counts for all right and left bams, compute means and totals
         # and obtain a SE with 2 assays, called means and totals.
         # As long as we need all the raw counts for the Fisher tests it is better to
         # keep them here as separate assays and compute (and keep) means in the computeRoar function.
         # The structure here is different from the single sample case for each condition, we will
         # see if it will be enough general to be used even there. Right now I ought to fallback
         # to that structure with the average counts.
         len <- length(rds@rightBams) + length(rds@leftBams)
         # Do we need preallocation for the list of matrixes?
         # x <- vector(mode = "list", length = 10)
         # ma <- matrix(nrow=5, ncol=2)
         # test <- list(ma, ma)
         # testse <- SummarizedExperiment(assays=test, rowData=gtfGRanges[c(1,2,3,4,5)], colData=DataFrame(row.names=c("a","b")))
         counts <- SummarizedExperiment(assays = matrix(nrow=length(rds@prePostCoords)/2, ncol=2),
                                    rowData=preCoords, 
                                    colData=DataFrame(row.names=c("pre","post"))
         )
         for (i in 1:length(rds@rightBams)) {
            rightSE <- summOv(rds@rightBams[[i]])
            rightSEpost <- summOvPost(rds@rightBams[[i]])
            assay(counts,i)[,"pre"] <- assays(rightSE)$counts[preElems,]
            assay(counts,i)[,"post"] <- assays(rightSEpost)$counts 
         }
         for (i in 1:length(rds@leftBams)) {
            leftSE <- summOv(rds@leftBams[[i]])
            leftSEpost <- summOvPost(rds@leftBams[[i]])
            assay(counts,i)[,"pre"] <- assays(leftSE)$counts[preElems,]
            assay(counts,i)[,"post"] <- assays(leftSEpost)$counts 
         }
         rds@counts <- counts
         names(assays(rds@counts)) <- "multiple_counts"
      }
      rowData(rds) <- rowData(se)
      colData(rds) <- colData(se)
      assays(rds) <- assays(se)
      names(assays(rds)) <- "counts"
      # We keep as "our" rds objects this SE, while those with counts are just slots.
      # In the case of multiple samples the primary object will be empty right now and will be filled
      # during the computeRoars step.
      rds@step <- 1;
      return(rds)
   }
)

setMethod("computeRoars", signature(rds="RoarDataset"),
   function(rds){
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
      if (length(rds@rightBams) == 1 && length(rds@leftBams) == 1) {
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
             goOn <- checkStep(rds, 2)
             if (!goOn[[1]]) {
                return(rds)
             }
             rds <- goOn[[2]]
             if (length(rds@rightBams) == 1 && length(rds@leftBams) == 1) {
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
      goOn <- checkStep(rds, 3)
      rds <- goOn[[2]]
      return(data.frame(row.names=sub("^\\s+","",sub("_POST","",elementMetadata(rds@postCoords)$gene_id)), 
                        mM_right=assay(rds,2)[,"right_pre"], 
                        mM_left=assay(rds,2)[,"right_post"],
                        roar=assay(rds,2)[,"left_pre"],
                        pval=assay(rds,2)[,"left_post"]))
   }
)

# This function will add to the totalResults dataframe RPKM gotten on the pre portions, then
# the user will be able to apply its preferred filtering/pvalue correction strategy.
# RPKM are gotten from counts over the PRE portions working on means across replicates.
# As total number of mapped reads we use the total number of reads mapped over all PRE portions.
setMethod("filteringInfoResults", signature(rds="RoarDataset"),
   function(rds){
      df <- totalResults(rds)
      preLen <- end(rowData(rds)) - start(rowData(rds))
      sumPreRight <- sum(assay(rds, 1)[,"right_pre"])
      sumPreLeft <- sum(assay(rds, 1)[,"left_pre"])
      df$rightValue <- (assay(rds, 1)[,"right_pre"]*1000000000)/(preLen*sumPreRight)
      df$leftValue <- (assay(rds, 1)[,"left_pre"]*1000000000)/(preLen*sumPreLeft)
      return(df)
   }
)

# Add a simple function to filter and compute corrected pvalues TODO

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