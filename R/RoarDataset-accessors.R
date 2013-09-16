# Methods for the RoarDataset class.
# R/AllGenerics.R contains the signature of the not overriding methods.

RoarDatasetFromFiles <- function(rightBams, leftBams, gtf) {
   # The format will be assumed using the file extension. Will work everytime?
   # Do we need to force a genome(eg. hg19)? It doesn't seem so.
   gtfGRanges<- import(gtf, asRangedData=FALSE)
   # here check _pre/ / _post
   ordered <- order(elementMetadata(gtfGRanges)$gene_id)
   gtfGRanges <- gtfGRanges[ordered]
   rightBamsGenomicAlignments <- lapply(rightBams, readGappedAlignments)
   leftBamsGenomicAlignments <- lapply(leftBams, readGappedAlignments)
   new("RoarDataset", rightBams=rightBamsGenomicAlignments, leftBams=leftBamsGenomicAlignments, 
       prePostCoords=gtfGRanges, step = 0, cores=1)
}

RoarDataset <- function(rightGappedAlign, leftGappedAlign, gtfGRanges) {
   if (length(rightGappedAlign) == 0 || length(leftGappedAlign) == 0) {
      stop("Lists of GappedAlignments could not be empty")
   }
   # here check _pre/ / _post
   ordered <- order(elementMetadata(gtfGRanges)$gene_id)
   gtfGRanges <- gtfGRanges[ordered]
   new("RoarDataset", rightBams=rightGappedAlign, leftBams=leftGappedAlign, 
       prePostCoords=gtfGRanges, step = 0, cores=1)
}

# Could have used setMethod("initialize", "xx",) but in this way should have had a gtf filename slot.

setMethod("countPrePost", signature(rds="RoarDataset", stranded="logical"),
   function(rds, stranded) {
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
         summarizeOverlaps(features=rds@prePostCoords, reads=x, ignore.strand=!stranded, mc.cores=rds@cores)
      }
      summOvPost <- function(x) {
         summarizeOverlaps(features=rds@postCoords, reads=x, ignore.strand=!stranded, mc.cores=rds@cores)
      } 
      
      # Now we need to keep means and totals of counts over PRE/POST for the two lists.
      # In the simpler case with a single alignment for both conditions we just keep the counts.
      preElems <- grep("_PRE$", elementMetadata(rds@prePostCoords)$gene_id)
      postElems <- grep("_POST$", elementMetadata(rds@prePostCoords)$gene_id)
      if (length(preElems)+length(postElems) != length(elementMetadata(rds@prePostCoords)$gene_id)) {
         stop("The prePostCoords given for this RoarDataset are wrong, some of the gene_id
              does not end in _PRE/_POST.")
      }
      # The number of elements has to be checked because some horrible special case with recycling
      # of vector elements in the comparison are possible.
      if (!all(sub("_PRE", "", elementMetadata(rds@prePostCoords)$gene_id[preElems]) ==
                     sub("_POST", "", elementMetadata(rds@prePostCoords)$gene_id[postElems]))) {
         stop("The prePostCoords given for this RoarDataset are wrong, not all prefixes of PRE-POST
              correspond.")
      }
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
         # This works but if a PRE portion of a gene overlaps with a POST portion of another one 
         # reads falling there are counted two times. FIXME
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
         countsLeft <- SummarizedExperiment(assays = matrix(nrow=length(rds@prePostCoords)/2, ncol=2),
                                            rowData=preCoords, 
                                            colData=DataFrame(row.names=c("pre","post"))
                                            )
         countsRight <- SummarizedExperiment(assays = matrix(nrow=length(rds@prePostCoords)/2, ncol=2),
                                             rowData=preCoords, 
                                             colData=DataFrame(row.names=c("pre","post"))
                                             )

         for (i in 1:length(rds@rightBams)) {
            rightSE <- summOv(rds@rightBams[[i]])
            rightSEpost <- summOvPost(rds@rightBams[[i]])
            assay(countsRight,i) <- matrix(nrow=length(rds@prePostCoords)/2, ncol=2)
            assay(countsRight,i)[,"pre"] <- assays(rightSE)$counts[preElems,]
            assay(countsRight,i)[,"post"] <- assays(rightSEpost)$counts 
         }
         for (i in 1:length(rds@leftBams)) {
            leftSE <- summOv(rds@leftBams[[i]])
            leftSEpost <- summOvPost(rds@leftBams[[i]])
            assay(countsLeft,i) <- matrix(nrow=length(rds@prePostCoords)/2, ncol=2)
            assay(countsLeft,i)[,"pre"] <- assays(leftSE)$counts[preElems,]
            assay(countsLeft,i)[,"post"] <- assays(leftSEpost)$counts 
         }
         rds@countsRight <- countsRight
         rds@countsLeft <- countsLeft
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
   function(rds) {
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
      preLen <- end(rowData(rds)) - start(rowData(rds)) + 1
      postLen <- end(rds@postCoords) - start(rds@postCoords) + 1
      # I had to add "+1" as long as the end coords are not inclusive.
      # Then the bam lengths to correct our lengths, ie: postLen+ReadLength-1
      if (length(rds@rightBams) > 1 || length(rds@leftBams) > 1) {
         # Compute means and put them in place for roar calculations (in the rds-se object).
         # meanAcrossAssays given a list/SimpleList of assays performs the mean on a 
         # given col name (ie. "pre").
         # lapply(assays(testse), function(x) { x[,"a"]})
         # rowMeans(as.data.frame(a))
         assay(rds,1)[,"right_pre"] <- meanAcrossAssays(assays(rds@countsRight), "pre") # here peak of memory usage?
         assay(rds,1)[,"right_post"] <- meanAcrossAssays(assays(rds@countsRight),"post")
         assay(rds,1)[,"left_pre"] <- meanAcrossAssays(assays(rds@countsLeft), "pre")
         assay(rds,1)[,"left_post"] <- meanAcrossAssays(assays(rds@countsLeft), "post")
         # Also the length correction should consider all the samples!
         lenRight <- unlist(lapply(rds@rightBams, qwidth))
         lenLeft <- unlist(lapply(rds@leftBams, qwidth))
         corrRight <- mean(lenRight)
         corrLeft <- mean(lenLeft)
      } else {
         corrRight <- mean(qwidth(rds@rightBams[[1]]))
         # qwidth(x): Returns an integer vector of length length(x) containing the length 
         # of the query *after* hard clipping (i.e. the length of the query sequence 
         # that is stored in the corresponding SAM/BAM record).
         corrLeft <- mean(qwidth(rds@leftBams[[1]])) 
      }
      # Ok, now if we had a single sample for both conditions we had the data charged in
      # countPrePost, otherwise we have the means (in the same SE/RDS object).
      # If there is a single sample for one condition and more than one for the other there is
      # a little (I hope) unuseful overload to get the mean for the single sample. 
      # The countsRight/Left slot are still kept as long as we will need them in computePvals.
      postLenRight <- postLen + corrRight - 1
      postLenLeft <- postLen + corrLeft - 1
      mMright <- (assay(rds,1)[,"right_pre"]*postLenRight)/(assay(rds,1)[,"right_post"]*preLen)-1
      mMleft <- (assay(rds,1)[,"left_pre"]*postLenLeft)/(assay(rds,1)[,"left_post"]*preLen)-1
      roar <- mMright / mMleft
      pVal <- rep(NA, length(roar))
      assay(rds,2) <- as.matrix(data.frame(right_pre=mMright, right_post=mMleft, left_pre=roar, left_post=pVal))
      names(assays(rds)) <- c("counts", "stats")
      rds@step <- 2
      return(rds)
   }
)

setMethod("computePvals", signature(rds="RoarDataset"),
   function(rds) {
      goOn <- checkStep(rds, 2)
      if (!goOn[[1]]) {
         return(rds)
      }
      rds <- goOn[[2]]
      if (length(rds@rightBams) == 1 && length(rds@leftBams) == 1) {
         if (cores(rds) == 1) {
            assay(rds,2)[,"left_post"] <- apply(assay(rds,1), 1, getFisher)
         } else {
            stop("TODO")
         }
      } else {      
         # We have raw counts in two SE slots (countsRight/Left) and need to
         # compute pvalues for every combination of right/left samples.
         # We need a function that given two assays returns the fisher pvalue
         # and we need to pass every combination there. The results will be but
         # in still another SE slot with a number of columns in the matrix equal to
         # the number of combinations. The product of all the pvalues will be put in the
         # rds/SE object in place of the pvalue for the single sample case.
         countsRightAssays <- assays(rds@countsRight)
         countsLeftAssays <- assays(rds@countsLeft)
         nRight <- length(countsRightAssays)
         nLeft <- length(countsLeftAssays)
         comparisons <- nRight*nLeft
         rds@pVals <- SummarizedExperiment(assays = matrix(nrow=length(rds@prePostCoords)/2, ncol=comparisons),
                                           rowData=rowData(rds), 
                                           # To obtain all combination of two vectors (x,y) in the right order:
                                           # as.vector(t(outer(x,y,paste,sep=""))
                                           colData=DataFrame(row.names=paste("pvalue_", 
                                                                             as.vector(t(outer(seq(1,nRight), seq(1,nLeft), paste, sep="_"))),
                                                                             sep=""))
                                          )
         # Ok, I know that we are in R, but these two for seems straightforward to me.
         for (i in 1:nRight) { # the y
            for (j in 1:nLeft) { # the x
               mat <- cbind(countsRightAssays[[i]], countsLeftAssays[[j]])
               assay(rds@pVals,1)[,nLeft*(i-1)+j] <- apply(mat, 1, getFisher)
            }
         }
         assay(rds, 2)[,"left_post"] <- apply(assay(rds@pVals,1), 1, prod)
         # Here in theory we could remove countsRight/Left slots, TODO check memory footprint and decide.
      }
      rds@step <- 3
      return(rds)
   }
)

setMethod("totalResults", signature(rds="RoarDataset"),
   function(rds) {
      goOn <- checkStep(rds, 3)
      rds <- goOn[[2]]
      res <- data.frame(row.names=sub("^\\s+","",sub("_POST","",elementMetadata(rds@postCoords)$gene_id)), 
                        mM_right=assay(rds,2)[,"right_pre"], 
                        mM_left=assay(rds,2)[,"right_post"],
                        roar=assay(rds,2)[,"left_pre"],
                        pval=assay(rds,2)[,"left_post"])
      if (length(rds@rightBams) != 1 || length(rds@leftBams) != 1) {
         pvals <- data.frame(assay(rds@pVals,1))
         colnames(pvals) <- rownames(colData(rds@pVals))
         res <- cbind(res, pvals)
      }
      return(res)
   }
)

# This function will add to the totalResults dataframe RPKM gotten on the pre portions, then
# the user will be able to apply its preferred filtering/pvalue correction strategy.
# RPKM are gotten from counts over the PRE portions working on means across replicates.
# As total number of mapped reads we use the total number of reads mapped over all PRE portions.
setMethod("fpkmResults", signature(rds="RoarDataset"),
   function(rds) {
      df <- totalResults(rds)
      preLen <- end(rowData(rds)) - start(rowData(rds)) + 1
      sumPreRight <- sum(assay(rds, 1)[,"right_pre"])
      sumPreLeft <- sum(assay(rds, 1)[,"left_pre"])
      df$rightValue <- (assay(rds, 1)[,"right_pre"]*1000000000)/(preLen*sumPreRight)
      df$leftValue <- (assay(rds, 1)[,"left_pre"]*1000000000)/(preLen*sumPreLeft)
      return(df)
   }
)

setMethod("countResults", signature(rds="RoarDataset"),
          function(rds) {
             df <- totalResults(rds)
             df$rightValue <- assay(rds, 1)[,"right_pre"]
             df$leftValue <- assay(rds, 1)[,"left_pre"]
             return(df)
          }
)

setMethod("standardFilter", signature(rds="RoarDataset", fpkmCutoff="numeric"),
   function(rds, fpkmCutoff) {
      # Here we need to: remove all genes with a mean FPKM <= fpkmCutoff, 
      # a negative/NA m/M-roar.
      # P-value correction? In the single samples case it seems sensible to do that,
      # otherwise we will report all pvalues (and correct their product.)
      # Due to chr by chr scanning bonferroni correction has been removed.
      df <- fpkmResults(rds)
      # mM_right, mM_left , roar columns filtering (< 0 / NA)
      # df <- subset(df, mM_right >= 0) # subset is ok for interactive use only
      df <- df[df$mM_right >= 0,]
      #df <- subset(df, mM_left >= 0)
      df <- df[df$mM_left >= 0,]
      #df <- subset(df, !is.na(roar))
      # Changed is.na to is.finite to avoid Inf/-Inf, did not add a unitTest as long as it's trivial.
      df <- df[is.finite(df$roar),]
      # rightValue/leftValue filtering (<= fpkmCutoff)
      #df <- subset(df, rightValue > fpkmCutoff)
      #df <- subset(df, leftValue > fpkmCutoff)
      df <- df[df$rightValue > fpkmCutoff,]
      df <- df[df$leftValue > fpkmCutoff,]
      #df$bonferroniPval <- p.adjust(df$pval, method="bonferroni")
      return(df)
   }                  
)

setMethod("pvalueFilter", signature(rds="RoarDataset", fpkmCutoff="numeric", pvalCutoff="numeric"),
   function(rds, fpkmCutoff, pvalCutoff) {
      df <- standardFilter(rds, fpkmCutoff)
      if (length(rds@rightBams) != 1 || length(rds@leftBams) != 1) {
         # In this case we add to df a col that says how many comparisons yielded
         # a pvalue < pvalCutoff.
         # esany <- apply(data, 1, function(x) {any(x[seq(1,12)] < 0.05)})
         cols <- grep("^pvalue_", colnames(df))
         sel <- apply(df, 1, function(x) {x[cols] < pvalCutoff})
         # This yields a transposed df with cols rows and TRUE/FALSE. ncol = nrows of df
         df$nUnderCutoff <- apply(sel, 2, function(x){length(x[x==TRUE])})
      } else {
         #df <- subset(df, bonferroniPval < pvalCutoff)  
         #df <- df[df$bonferroniPval < pvalCutoff,]
         df <- df[df$pval < pvalCutoff,]
      }   
      return(df)
   }                  
)



# Simple getters and setters. Arf Arf!
# XXX TODO
setMethod("cores",  signature(rds="RoarDataset"),
 function(rds) {
    return(rds@cores)
 }
)
 
# setReplaceMethod("cores",  signature(rds="RoarDataset", value="numeric"),
#  function(rds, value) {
#     if (value < 1) {
#        stop("You can specify only a positive number of cores")
#     }
#     rds@cores <- value
#     return(rds)
#  }
# )