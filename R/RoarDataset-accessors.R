# Methods for the RoarDataset class.
# R/AllGenerics.R contains the signatures of the not overriding methods.

RoarDatasetFromFiles <- function(treatmentBams, controlBams, gtf) {
   # The format will be assumed using the file extension. Will work everytime?
   # Do we need to force a genome(eg. hg19)? It doesn't seem so.
   gtfGRanges<- import(gtf, asRangedData=FALSE)
   ordered <- order(elementMetadata(gtfGRanges)$gene_id)
   gtfGRanges <- gtfGRanges[ordered]
   treatmentBamsGenomicAlignments <- lapply(treatmentBams, readGAlignments)
   controlBamsGenomicAlignments <- lapply(controlBams, readGAlignments)
   new("RoarDataset", treatmentBams=treatmentBamsGenomicAlignments, controlBams=controlBamsGenomicAlignments, 
       prePostCoords=gtfGRanges, step = 0, paired=FALSE, cores=1)
}

RoarDataset <- function(treatmentGappedAlign, controlGappedAlign, gtfGRanges) {
   if (length(treatmentGappedAlign) == 0 || length(controlGappedAlign) == 0) {
      stop("Lists of GAlignments could not be empty")
   }
   ordered <- order(elementMetadata(gtfGRanges)$gene_id)
   gtfGRanges <- gtfGRanges[ordered]
   new("RoarDataset", treatmentBams=treatmentGappedAlign, controlBams=controlGappedAlign, 
       prePostCoords=gtfGRanges, step = 0, paired=FALSE, cores=1)
}

# Could have used setMethod("initialize", "xx",) but in this way should have had a gtf filename slot.
# setValidity instead of checks inside countPrePost? But this would be inefficient/add slots.

# I removed stranded="logical" from the signature because it has a default value, which
# is also set in setGeneric.
setMethod("countPrePost", signature(rds="RoarDataset"),
   function(rds, stranded=FALSE) {
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
         stop("The prePostCoords given for this RoarDataset are wrong: some of the gene_id
              does not end in _PRE/_POST.")
      }
      if (length(preElems) != length(postElems)) {
         stop("The prePostCoords given for this RoarDataset are wrong: the number of PRE is different
              from the number of POST.")
      }
      # The number of elements has to be checked because some horrible special case with recycling
      # of vector elements in the comparison are possible.
      if (!all(sub("_PRE", "", elementMetadata(rds@prePostCoords)$gene_id[preElems]) ==
                     sub("_POST", "", elementMetadata(rds@prePostCoords)$gene_id[postElems]))) {
         stop("The prePostCoords given for this RoarDataset are wrong: not all prefixes of PRE-POST
              correspond.")
      }
      # Check uniqueness of gene_id.
      geneIds <- sub("_PRE", "", elementMetadata(rds@prePostCoords)$gene_id[preElems])
      if (!all(geneIds==make.unique(geneIds))) {
         stop("The prePostCoords given for this RoarDataset are wrong: gene_ids (prefixes of PRE-POST)
               are not unique.")
      }
      preCoords <- rds@prePostCoords[preElems,]
      rds@postCoords <- rds@prePostCoords[postElems,]
      se <- SummarizedExperiment(assays = matrix(nrow=length(rds@prePostCoords)/2, ncol=4),
                                 rowData=preCoords, 
                                 colData=DataFrame(row.names=c("treatment_pre","treatment_post","control_pre", "control_post"))
                                 )
      if (length(rds@treatmentBams) == 1 && length(rds@controlBams) == 1) {
         # We obtain counts for both conditions on PRE and POST coords.
         treatmentSE <- summOv(rds@treatmentBams[[1]])
         controlSE <- summOv(rds@controlBams[[1]])    
         # We obtain counts for both conditions only on POST coords (notwithstanding overlap with PRE) -
         # this is needed to implement the prefer-POST policy where a read overlapping PRE and POST
         # is counted on the POST portion (it has to be derived from a long isoform).
         # Note that in the first case we need to use the whole gtf with PRE and POST to avoid
         # assigning to PRE the overlapping reads, even if the counts on POST gotten afterwards
         # will be discarded.
         # This works but if a PRE portion of a gene overlaps with a POST portion of another one 
         # reads falling there are counted two times.
         treatmentSEpost <- summOvPost(rds@treatmentBams[[1]])
         controlSEpost <- summOvPost(rds@controlBams[[1]])
         assay(se,1)[,"treatment_pre"] <- assays(treatmentSE)$counts[preElems,]
         assay(se,1)[,"treatment_post"] <- assays(treatmentSEpost)$counts 
         assay(se,1)[,"control_pre"] <- assays(controlSE)$counts[preElems,]
         assay(se,1)[,"control_post"] <- assays(controlSEpost)$counts
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
         len <- length(rds@treatmentBams) + length(rds@controlBams)
         # Do we need preallocation for the list of matrixes?
         # x <- vector(mode = "list", length = 10)
         # ma <- matrix(nrow=5, ncol=2)
         # test <- list(ma, ma)
         # testse <- SummarizedExperiment(assays=test, rowData=gtfGRanges[c(1,2,3,4,5)], colData=DataFrame(row.names=c("a","b")))
         countsControl <- SummarizedExperiment(assays = matrix(nrow=length(rds@prePostCoords)/2, ncol=2),
                                            rowData=preCoords, 
                                            colData=DataFrame(row.names=c("pre","post"))
                                            )
         countsTreatment <- SummarizedExperiment(assays = matrix(nrow=length(rds@prePostCoords)/2, ncol=2),
                                             rowData=preCoords, 
                                             colData=DataFrame(row.names=c("pre","post"))
                                             )

         for (i in 1:length(rds@treatmentBams)) {
            treatmentSE <- summOv(rds@treatmentBams[[i]])
            treatmentSEpost <- summOvPost(rds@treatmentBams[[i]])
            assay(countsTreatment,i) <- matrix(nrow=length(rds@prePostCoords)/2, ncol=2)
            assay(countsTreatment,i)[,"pre"] <- assays(treatmentSE)$counts[preElems,]
            assay(countsTreatment,i)[,"post"] <- assays(treatmentSEpost)$counts 
         }
         for (i in 1:length(rds@controlBams)) {
            controlSE <- summOv(rds@controlBams[[i]])
            controlSEpost <- summOvPost(rds@controlBams[[i]])
            assay(countsControl,i) <- matrix(nrow=length(rds@prePostCoords)/2, ncol=2)
            assay(countsControl,i)[,"pre"] <- assays(controlSE)$counts[preElems,]
            assay(countsControl,i)[,"post"] <- assays(controlSEpost)$counts 
         }
         rds@countsTreatment <- countsTreatment
         rds@countsControl <- countsControl
      }
      rowData(rds) <- rowData(se)
      colData(rds) <- colData(se)
      assays(rds) <- assays(se)
      names(assays(rds)) <- "counts"
      # We keep as "our" rds objects this SE, while those with counts are just slots.
      # In the case of multiple samples the primary object will be empty treatment now and will be filled
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
      # roar is the m/M of treatment condition divided by the m/M of the control one.
      # m/M = ((Lpost*Cpre)/(Lpre*Cpost))-1
      # Negative m/M are discarded.
      #roar = (m/M_treatment)/(m/M_control)
      # We must obtain the list of lengths from rds@postCoords and rowData(rds) (which is pre).
      preLen <- end(rowData(rds)) - start(rowData(rds)) + 1
      postLen <- end(rds@postCoords) - start(rds@postCoords) + 1
      # I had to add "+1" as long as the end coords are not inclusive.
      # Then the bam lengths to correct our lengths, ie: postLen+ReadLength-1
      if (length(rds@treatmentBams) > 1 || length(rds@controlBams) > 1) {
         # Compute means and put them in place for roar calculations (in the rds-se object).
         # meanAcrossAssays given a list/SimpleList of assays performs the mean on a 
         # given col name (ie. "pre").
         # lapply(assays(testse), function(x) { x[,"a"]})
         # rowMeans(as.data.frame(a))
         assay(rds,1)[,"treatment_pre"] <- meanAcrossAssays(assays(rds@countsTreatment), "pre") # here peak of memory usage?
         assay(rds,1)[,"treatment_post"] <- meanAcrossAssays(assays(rds@countsTreatment),"post")
         assay(rds,1)[,"control_pre"] <- meanAcrossAssays(assays(rds@countsControl), "pre")
         assay(rds,1)[,"control_post"] <- meanAcrossAssays(assays(rds@countsControl), "post")
         # Also the length correction should consider all the samples!
         lenTreatment <- unlist(lapply(rds@treatmentBams, qwidth))
         lenControl <- unlist(lapply(rds@controlBams, qwidth))
         corrTreatment <- mean(lenTreatment)
         corrControl <- mean(lenControl)
      } else {
         corrTreatment <- mean(qwidth(rds@treatmentBams[[1]]))
         # qwidth(x): Returns an integer vector of length length(x) containing the length 
         # of the query *after* hard clipping (i.e. the length of the query sequence 
         # that is stored in the corresponding SAM/BAM record).
         corrControl <- mean(qwidth(rds@controlBams[[1]])) 
      }
      # Ok, now if we had a single sample for both conditions we had the data charged in
      # countPrePost, otherwise we have the means (in the same SE/RDS object).
      # If there is a single sample for one condition and more than one for the other there is
      # a little (I hope) unuseful overload to get the mean for the single sample. 
      # The countsTreatment/control slot are still kept as long as we will need them in computePvals.
      postLenTreatment <- postLen + corrTreatment - 1
      postLenControl <- postLen + corrControl - 1
      mMtreatment <- (assay(rds,1)[,"treatment_pre"]*postLenTreatment)/(assay(rds,1)[,"treatment_post"]*preLen)-1
      mMcontrol <- (assay(rds,1)[,"control_pre"]*postLenControl)/(assay(rds,1)[,"control_post"]*preLen)-1
      roar <- mMtreatment / mMcontrol
      pVal <- rep(NA, length(roar))
      assay(rds,2) <- as.matrix(data.frame(treatment_pre=mMtreatment, treatment_post=mMcontrol, control_pre=roar, control_post=pVal))
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
      if (length(rds@treatmentBams) == 1 && length(rds@controlBams) == 1) {
         if (cores(rds) == 1) {
            assay(rds,2)[,"control_post"] <- apply(assay(rds,1), 1, getFisher)
         } else {
            stop("TODO")
         }
      } else {      
         # We have raw counts in two SE slots (countsTreatment/Control) and need to
         # compute pvalues for every combination of treatment/control samples.
         # We need a function that given two assays returns the fisher pvalue
         # and we need to pass every combination there. The results will be but
         # in still another SE slot with a number of columns in the matrix equal to
         # the number of combinations. The product of all the pvalues will be put in the
         # rds/SE object in place of the pvalue for the single sample case.
         countsTreatmentAssays <- assays(rds@countsTreatment)
         countsControlAssays <- assays(rds@countsControl)
         nTreatment <- length(countsTreatmentAssays)
         nControl <- length(countsControlAssays)
         comparisons <- nTreatment*nControl
         rds@pVals <- SummarizedExperiment(assays = matrix(nrow=length(rds@prePostCoords)/2, ncol=comparisons),
                                           rowData=rowData(rds), 
                                           # To obtain all combination of two vectors (x,y) in the treatment order:
                                           # as.vector(t(outer(x,y,paste,sep=""))
                                           colData=DataFrame(row.names=paste("pvalue_", 
                                                                             as.vector(t(outer(seq(1,nTreatment), seq(1,nControl), paste, sep="_"))),
                                                                             sep=""))
                                          )
         # Ok, I know that we are in R, but these two for seems straightforward to me.
         for (i in 1:nTreatment) { # the y
            for (j in 1:nControl) { # the x
               mat <- cbind(countsTreatmentAssays[[i]], countsControlAssays[[j]])
               assay(rds@pVals,1)[,nControl*(i-1)+j] <- apply(mat, 1, getFisher)
            }
         }
         assay(rds, 2)[,"control_post"] <- apply(assay(rds@pVals,1), 1, prod)
         # Here in theory we could remove countsTreatment/control slots, TODO check memory footprint and decide.
      }
      rds@step <- 3
      return(rds)
   }
)

setMethod("computePairedPvals", signature(rds="RoarDataset", treatmentSamples="list", controlSamples="list"),
   function(rds, treatmentSamples, controlSamples) {
      goOn <- checkStep(rds, 2)
      if (!goOn[[1]]) {
         return(rds)
      }
      rds <- goOn[[2]]
      if (length(rds@treatmentBams) == 1 && length(rds@controlBams) == 1) {
         # ERROR conditions
      } else {      
         # We have raw counts in two SE slots (countsTreatment/Control) and need to
         # compute pvalues for every combination of treatment/control samples.
         # We need a function that given two assays returns the fisher pvalue
         # and we need to pass every combination there. The results will be but
         # in still another SE slot with a number of columns in the matrix equal to
         # the number of combinations. The product of all the pvalues will be put in the
         # rds/SE object in place of the pvalue for the single sample case.
         countsTreatmentAssays <- assays(rds@countsTreatment)
         countsControlAssays <- assays(rds@countsControl)
         nTreatment <- length(countsTreatmentAssays)
         nControl <- length(countsControlAssays)
         comparisons <- nTreatment*nControl
         rds@pVals <- SummarizedExperiment(assays = matrix(nrow=length(rds@prePostCoords)/2, ncol=comparisons),
                                           rowData=rowData(rds), 
                                           # To obtain all combination of two vectors (x,y) in the treatment order:
                                           # as.vector(t(outer(x,y,paste,sep=""))
                                           colData=DataFrame(row.names=paste("pvalue_", 
                                                                             as.vector(t(outer(seq(1,nTreatment), seq(1,nControl), paste, sep="_"))),
                                                                             sep=""))
                                          )
         # Ok, I know that we are in R, but these two for seems straightforward to me.
         for (i in 1:nTreatment) { # the y
            for (j in 1:nControl) { # the x
               mat <- cbind(countsTreatmentAssays[[i]], countsControlAssays[[j]])
               assay(rds@pVals,1)[,nControl*(i-1)+j] <- apply(mat, 1, getFisher)
            }
         }
         assay(rds, 2)[,"control_post"] <- apply(assay(rds@pVals,1), 1, prod)
         # Here in theory we could remove countsTreatment/control slots, TODO check memory footprint and decide.
      }
      rds@paired <- TRUE
      rds@step <- 3
      return(rds)
   }
)

setMethod("totalResults", signature(rds="RoarDataset"),
   function(rds) {
      goOn <- checkStep(rds, 3)
      rds <- goOn[[2]]
      res <- data.frame(row.names=sub("^\\s+","",sub("_POST","",elementMetadata(rds@postCoords)$gene_id)), 
                        mM_treatment=assay(rds,2)[,"treatment_pre"], 
                        mM_control=assay(rds,2)[,"treatment_post"],
                        roar=assay(rds,2)[,"control_pre"],
                        pval=assay(rds,2)[,"control_post"])
      if (length(rds@treatmentBams) != 1 || length(rds@controlBams) != 1) {
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
      dat <- totalResults(rds)
      preLen <- end(rowData(rds)) - start(rowData(rds)) + 1
      sumPreTreatment <- sum(assay(rds, 1)[,"treatment_pre"])
      sumPreControl <- sum(assay(rds, 1)[,"control_pre"])
      dat$treatmentValue <- (assay(rds, 1)[,"treatment_pre"]*1000000000)/(preLen*sumPreTreatment)
      dat$controlValue <- (assay(rds, 1)[,"control_pre"]*1000000000)/(preLen*sumPreControl)
      return(dat)
   }
)

setMethod("countResults", signature(rds="RoarDataset"),
   function(rds) {
      dat <- totalResults(rds)
      dat$treatmentValue <- assay(rds, 1)[,"treatment_pre"]
      dat$controlValue <- assay(rds, 1)[,"control_pre"]
      return(dat)
   }
)

setMethod("standardFilter", signature(rds="RoarDataset", fpkmCutoff="numeric"),
   function(rds, fpkmCutoff) {
      # Here we need to: remove all genes with a mean FPKM <= fpkmCutoff, 
      # a negative/NA m/M-roar.
      # P-value correction? In the single samples case it seems sensible to do that,
      # otherwise we will report all pvalues (and correct their product.)
      # Due to chr by chr scanning bonferroni correction has been removed.
      dat <- fpkmResults(rds)
      # mM_treatment, mM_control , roar columns filtering (< 0 / NA)
      # dat <- subset(dat, mM_treatment >= 0) # subset is ok for interactive use only
      dat <- dat[is.finite(dat$mM_treatment) & dat$mM_treatment >= 0,]
      #dat <- subset(dat, mM_control >= 0)
      dat <- dat[is.finite(dat$mM_control) & dat$mM_control >= 0,]
      #dat <- subset(dat, !is.na(roar))
      # Changed is.na to is.finite to avoid Inf/-Inf, did not add a unitTest as long as it's trivial.
      dat <- dat[is.finite(dat$roar),]
      # treatmentValue/controlValue filtering (<= fpkmCutoff)
      #dat <- subset(dat, treatmentValue > fpkmCutoff)
      #dat <- subset(dat, controlValue > fpkmCutoff)
      dat <- dat[dat$treatmentValue > fpkmCutoff,]
      dat <- dat[dat$controlValue > fpkmCutoff,]
      return(dat)
   }                  
)

setMethod("pvalueFilter", signature(rds="RoarDataset", fpkmCutoff="numeric", pvalCutoff="numeric"),
   function(rds, fpkmCutoff, pvalCutoff) {
      dat <- standardFilter(rds, fpkmCutoff)
      if ((length(rds@treatmentBams) != 1 || length(rds@controlBams) != 1) && !rds@paired) {
         # In this case we add to dat a col that says how many comparisons yielded
         # a pvalue < pvalCutoff.
         # esany <- apply(data, 1, function(x) {any(x[seq(1,12)] < 0.05)})
         if(nrow(dat) != 0) {
            cols <- grep("^pvalue_", colnames(dat))
            sel <- apply(dat, 1, function(x) {x[cols] < pvalCutoff})
            # This yields a transposed dat with cols rows and TRUE/FALSE. ncol = nrows of dat
            dat$nUnderCutoff <- apply(sel, 2, function(x){length(x[x==TRUE])})
         }
      } else {
         dat <- dat[dat$pval < pvalCutoff,]
      }   
      return(dat)
   }                  
)

setMethod("pvalueCorrectFilter", signature(rds="RoarDataset", fpkmCutoff="numeric", pvalCutoff="numeric", method="character"),
   function(rds, fpkmCutoff, pvalCutoff, method) {
      dat <- standardFilter(rds, fpkmCutoff)
      if ((length(rds@treatmentBams) != 1 || length(rds@controlBams) != 1) && !rds@paired) {
         # In this case we add to dat a col that says how many comparisons yielded
         # a pvalue < pvalCutoff.
         # esany <- apply(data, 1, function(x) {any(x[seq(1,12)] < 0.05)})
         if(nrow(dat) != 0) {
            cols <- grep("^pvalue_", colnames(dat))
            sel <- apply(dat, 1, function(x) {x[cols] < pvalCutoff})
            # This yields a transposed dat with cols rows and TRUE/FALSE. ncol = nrows of dat
            dat$nUnderCutoff <- apply(sel, 2, function(x){length(x[x==TRUE])})
         }
      } else {
         dat$pval <- p.adjust(dat$pval, method=method)
         dat <- dat[dat$pval < pvalCutoff,]
      }   
      return(dat)
   }                 
)

# Simple getters and setters.
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

# Here signature does not work. We are overriding show, that's why?
setMethod("show", "RoarDataset",
   function(object) {
      cat("RoarDataset object\n")
      cat("N. of treatment alignments:", length(object@treatmentBams), "\n")
      cat("N. of control alignments:", length(object@controlBams), "\n")
      cat("N. of genes in study:", length(object@prePostCoords)/2 , "\n")
      cat("N. of cores:", object@cores, "\n")
      cat("Analysis step reached [0-3]:", object@step, "\n")
      cat("\n")     
   }
)