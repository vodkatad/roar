RoarDatasetMultipleAPAFromFiles <- function(treatmentBams, controlBams, gtf) {
   gtfGRanges <- import(gtf)
   treatmentBamsGenomicAlignments <- lapply(treatmentBams, readGAlignments)
   controlBamsGenomicAlignments <- lapply(controlBams, readGAlignments)
   RoarDatasetMultipleAPA(treatmentBamsGenomicAlignments, controlBamsGenomicAlignments, gtfGRanges)
}

RoarDatasetMultipleAPA <- function(treatmentBamsGenomicAlignments, controlBamsGenomicAlignments, gtfGRanges) {
   apas_melted <- gtfGRanges[mcols(gtfGRanges)$type=="apa"]
   genes_melted <- gtfGRanges[mcols(gtfGRanges)$type=="gene"]
   mcols(apas_melted)$gene <- sapply(strsplit(mcols(apas_melted)$apa, '_', 
                                              fixed=TRUE),
                                     function(x) { x[length(x)]})
   genes_ids <- sort(unique(as.character(mcols(genes_melted)$gene)))
   genes_ids_apas <- sort(unique(mcols(apas_melted)$gene))
   if (!all(genes_ids==genes_ids_apas)) {
      stop("All the genes in the gtf should have at least one apa")
   }
   genes <- do.call(GRangesList, sapply(genes_ids, 
                                        function(x) {genes_melted[mcols(genes_melted)$gene==x]}))
   apas <- do.call(GRangesList, sapply(genes_ids, 
                                       function(x) {apas_melted[mcols(apas_melted)$gene==x]}))
   names(apas) <- genes_ids
   names(genes) <- genes_ids
   new("RoarDatasetMultipleAPA", treatmentBams=treatmentBamsGenomicAlignments, 
       controlBams=controlBamsGenomicAlignments, 
       geneCoords=genes, apaCoords=apas, step=0, paired=FALSE, cores=1)
}

setMethod("countPrePost", signature(rds="RoarDatasetMultipleAPA"),
   function(rds, stranded=FALSE) {
      allFragmentsAndPrePostDef <- mapply(getApaGenesFractions,
                              rds@geneCoords, 
                              rds@apaCoords)
      rds@fragments <- GRangesList(allFragmentsAndPrePostDef[1,]) 
      #Right subsetting? Yes because in this situation mapply gives us a matrix with 
      #columns == gene names and rows that are fragments and prePostDef (SIMPLIFY=TRUE).
      rds@prePostDef <- allFragmentsAndPrePostDef[2,]
      # A GRangesList with GRanges for all fragments defining pre/post
      # in a gene and a list with info on APA choices and which
      # fragments are to be summed.
      ResizeReadsPlus <- function(reads, width=1, fix="end", ...) {
         reads <- as(reads, "GRanges")
         #stopifnot(all(strand(reads) != "*"))
         # Bad and ugly and will do horrible things for stranded data but still...
         # FIXME: if stranded=T different resize functions that do not change read strand but
         # move them always at their ends (should be the right call because if they are on a given
         # strand they align on a feature on that strand).
         strand(reads) <- rep("+", length(reads))
         resize(reads, width=width, fix=fix, ...)
      }
      ResizeReadsMinus <- function(reads, width=1, fix="start", ...) {
         reads <- as(reads, "GRanges")
         #stopifnot(all(strand(reads) != "*"))
         strand(reads) <- rep("+", length(reads))
         resize(reads, width=width, fix=fix, ...)
      }
      # Does not work as expected as long as for GRangesList counts
      # are collapsed: need to unlist.
      summOv <- function(x) {
         frag <-  unlist(rds@fragments)
         # We unlist because 
         # When a GRanges is supplied, each row is considered a feature.
         # When a GRangesList is supplied, each higher list-level 
         # is considered a feature. 
         # This distinction is important when defining overlaps.
         # Having genes as feature would help? 
         # No, multiple counts would happen I believe.
         featPlus <- frag[strand(frag)=="+"]
         featMinus <- frag[strand(frag)=="-"]
         plus <- summarizeOverlaps(features=featPlus, reads=x, 
                                    ignore.strand=!stranded,
                                    preprocess.reads=ResizeReadsPlus)
         minus <- summarizeOverlaps(features=featMinus, reads=x, 
                                    ignore.strand=!stranded, 
                                    preprocess.reads=ResizeReadsMinus)
         return(rbind(plus, minus)) # Does this work as expected? Yes. 
         # A more insightful comment would be what I expected...fuck you past Elena.
         # It's like an rbind for counts as long as we have only counts as assays columns 
         # and rowRanges are our genes. The order of fragments inside genes is kept while
         # not their relative order but we will get them from names later on so it's ok. 
      }
      if (length(rds@treatmentBams) == 1 && length(rds@controlBams) == 1) {
         # We obtain counts for both conditions on gene fragments and them sum 
         # them to obtain PRE/POST counts.
         treatmentSE <- summOv(rds@treatmentBams[[1]])
         controlSE <- summOv(rds@controlBams[[1]])    
         rds <- generateRoarsSingleBam(rds, treatmentSE, controlSE)
         rds@corrTreatment <- mean(qwidth(rds@treatmentBams[[1]]))
         rds@corrControl <- mean(qwidth(rds@controlBams[[1]]))
      } else {
         countsControl <- vector(mode = "list", length = length(rds@controlBams))
         countsTreatment <- vector(mode = "list", length = length(rds@treatmentBams))
         for (i in 1:length(rds@treatmentBams)) {
            countsTreatment[[i]] <- summOv(rds@treatmentBams[[i]])
         }
         for (i in 1:length(rds@controlBams)) {
            countsControl[[i]] <- summOv(rds@controlBams[[i]])
         }
         rds <- generateRoarsMultipleBam(rds, countsTreatment, countsControl)
         rds@corrTreatment <- mean(unlist(lapply(rds@treatmentBams, qwidth)))
         rds@corrControl <- mean(unlist(lapply(rds@controlBams, qwidth)))
      }
      return(rds)
   }       
)

setMethod("generateRoarsSingleBam", signature(rds="RoarDatasetMultipleAPA", 
                                              "RangedSummarizedExperiment",
                                              "RangedSummarizedExperiment"),
   function(rds, treatmentSE, controlSE)
   {
      # treatmentSE and controlSE are MoreArgs?
      #rds@roars <- mapply(createRoarSingleBAM, rds@fragments, rds@prePostDef,
      #                treatmentSE, controlSE)
      # Could be:
      rds@roars <- lapply(names(rds@fragments), createRoarSingleBam,
                          rds, treatmentSE, controlSE)
      names(rds@roars) <- names(rds@fragments)
      # to have the names...other ways? Store them in fragments or prePostDef
      # in an accessible way another time seems a waste of space.
      # Will have to compare times!
      return(rds)
   }
)

setMethod("generateRoarsMultipleBam", signature(rds="RoarDatasetMultipleAPA", 
                                              "list",
                                              "list"),
          function(rds, treatmentSE, controlSE)
          {
             # treatmentSE and controlSE are MoreArgs?
             #rds@roars <- mapply(createRoarSingleBAM, rds@fragments, rds@prePostDef,
             #                treatmentSE, controlSE)
             # Could be:
             rds@roars <- lapply(names(rds@fragments), createRoarMultipleBam,
                                 rds, treatmentSE, controlSE)
             names(rds@roars) <- names(rds@fragments)
             # to have the names...other ways? Store them in fragments or prePostDef
             # in an accessible way another time seems a waste of space.
             # Will have to compare times!
             return(rds)
          }
)

setMethod("computeRoars", signature(rds="RoarDatasetMultipleAPA"),
      function(rds) 
      {
         rds@roars <- lapply(rds@roars, computeRoars, 
                             rds@corrTreatment, rds@corrControl)
         return(rds)
      }
)

setMethod("computePvals", signature(rds="RoarDatasetMultipleAPA"),
      function(rds)
      {
         rds@roars <- lapply(rds@roars, computePvals)
         return(rds)
      }
)

setMethod("computePairedPvals", 
            signature(rds="RoarDatasetMultipleAPA",
            treatmentSamples="numeric", controlSamples="numeric"),
      function(rds, treatmentSamples, controlSamples) 
      {
         rds@roars <- lapply(rds@roars, computePairedPvals, 
                             treatmentSamples, controlSamples)
         return(rds)
      }
)

setMethod("totalResults", signature(rds="RoarDatasetMultipleAPA"),
      function(rds) 
      {
         totRes <- lapply(rds@roars, totalResults)
         r <- do.call(rbind, totRes)
         rownames(r) <- unlist(lapply(names(totRes), function(x) {
                                 paste(x, rownames(totRes[[x]]), sep="_") }))
         return(r)
         # For names do.call/rbind added APA ids only when there were multiple
         # choices for the same gene. There could be best ways to add names.
      }
)

setMethod("countResults", signature(rds="RoarDatasetMultipleAPA"),
      function(rds) 
      {
         # We cannot call fpkmResults on the roars objects because in this
         # case we want FPKM for whole genes. Therefore we
         # need to sum all the fragments counts. We work on every rds@roars.
         fpkmList <- lapply(rds@roars, sumRoarCounts)
         fpkmDf <- do.call(rbind, fpkmList)
         #t(sapply(rds@roars, s1)) # the same
         totRes <- totalResults(rds)
         totRes$geneName <- sapply(rownames(totRes), function(x) unlist(strsplit(x,"_", fixed=T))[1])
         res <- merge(totRes, fpkmDf, by.x="geneName", by.y="row.names")
         rownames(res) <- rownames(totRes)
         res$geneName <- NULL
         return(res)
      }
)

setMethod("fpkmResults", signature(rds="RoarDatasetMultipleAPA"),
      function(rds) 
      {
         counts <- countResults(rds)
         lengths <- sapply(rds@fragments, function(x) { sum(mcols(x)$length)})
         # We put in rds@fragments mcol length the lengths of the exon portions of 
         # the fragments that we are considering.
         lengthsdf <- data.frame(length=lengths)
         counts$genes <- sapply(rownames(counts), function(x) unlist(strsplit(x,"_", fixed=T))[1])
         countsonly <- counts[,c("counts_treatment", "counts_control","genes")]
         countsonly <- unique(countsonly)
         sumTreatment <- sum(countsonly[,"counts_treatment"])
         sumControl <- sum(countsonly[,"counts_control"])
         dat <- merge(countsonly, lengthsdf, by.x="genes", by.y="row.names")
         dat$treatmentFpkm <- (dat[,"counts_treatment"]*1000000000)/(dat[,"length"]*sumTreatment)
         dat$controlFpkm <- (dat[,"counts_control"]*1000000000)/(dat[,"length"]*sumControl)
         res <- merge(counts, dat, by="genes")
         res$genes <- NULL
         res$counts_control.x <- NULL
         res$counts_control.y <- NULL
         res$counts_treatment.x <- NULL
         res$counts_treatment.y <- NULL
         #res$length <- NULL
         rownames(res) <- rownames(counts)
         return(res)
      }
)

setMethod("pvalueFilter", signature(rds="RoarDatasetMultipleAPA", fpkmCutoff="numeric", pvalCutoff="numeric"),
      function(rds, fpkmCutoff, pvalCutoff) {
         dat <- standardFilter(rds, fpkmCutoff)
         if (nrow(dat) > 0) {
            resby <- by(dat, INDICES=as.factor(sapply(rownames(dat), function(x) unlist(strsplit(x,"_", fixed=T))[1])), 
                        FUN = function(x) {res<-x[with(x, order(x[,4])),]; res[1,]})
            dat <- do.call(rbind, resby)
            # For each gene we select the APA choice that is associated with the smallest p-value (after fpkm filtering)
            # then proceed exactly like roar [NDRY ALERT FIXME].
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
         } else {
            return(data.frame())
         }
      }                  
)

setMethod("pvalueCorrectFilter", signature(rds="RoarDatasetMultipleAPA", fpkmCutoff="numeric", pvalCutoff="numeric", method="character"),
      function(rds, fpkmCutoff, pvalCutoff, method) 
      {
         dat <- standardFilter(rds, fpkmCutoff)
         if (nrow(dat) > 0) {   
            resby <- by(dat, INDICES=as.factor(sapply(rownames(dat), function(x) unlist(strsplit(x,"_", fixed=T))[1])), 
                        FUN = function(x) {res<-x[with(x, order(x[,4])),]; res[1,]})
            dat <- do.call(rbind, resby)
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
         } else {
            return(data.frame())
         }    
      }
)