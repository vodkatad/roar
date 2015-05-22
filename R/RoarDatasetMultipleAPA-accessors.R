RoarDatasetMultipleAPAFromFiles <- function(treatmentBams, controlBams, gtf) {
   gtfGRanges<- import(gtf, asRangedData=FALSE)
   apas_melted <- gtfGRanges[mcols(gtfGRanges)$type=="apa"]
   genes_melted <- gtfGRanges[mcols(gtfGRanges)$type=="gene"]
   mcols(apas_melted)$gene <- sapply(strsplit(mcols(apas_melted)$apa, '_', 
                                              fixed=TRUE),
                              function(x) { x[length(x)]})
   genes_ids <- sort(unique(mcols(genes_melted)$gene))
   genes_ids_apas <- sort(unique(mcols(apas_melted)$gene))
   if (!all(genes_ids==genes_ids_apas)) {
      stop("Lists of GAlignments could not be empty")
   }
   genes <- do.call(GRangesList, sapply(genes_ids, 
                     function(x) {genes_melted[mcols(genes_melted)$gene==x]}))
   apas <- do.call(GRangesList, sapply(genes_ids, 
                     function(x) {apas_melted[mcols(apas_melted)$gene==x]}))
   names(apas) <- genes_ids
   names(genes) <- genes_ids
   treatmentBamsGenomicAlignments <- lapply(treatmentBams, readGAlignments)
   controlBamsGenomicAlignments <- lapply(controlBams, readGAlignments)
   new("RoarDatasetMultipleAPA", treatmentBams=treatmentBamsGenomicAlignments, 
       controlBams=controlBamsGenomicAlignments, 
       geneCoords=genes, apaCoords=apas, step=0, paired=FALSE, cores=1)
}

setMethod("countPrePost", signature(rds="RoarDatasetMultipleAPA"),
   function(rds, stranded=FALSE) {
      allFragmentsAndPrePostDef <- mapply(getApaGenesFractions,
                              rds@geneCoords, 
                              rds@apaCoords)
      rds@fragments <- allFragmentsAndPrePostDef[1,] #Right subsetting?
      rds@prePostDef <- allFragmentsAndPrePostDef[2,]
      # A GRangesList with GRanges for all fragments defining pre/post
      # in a gene and a list with info on APA choices and which
      # fragments are to be summed.
      ResizeReadsPlus <- function(reads, width=1, fix="end", ...) {
         reads <- as(reads, "GRanges")
         #stopifnot(all(strand(reads) != "*"))
         resize(reads, width=width, fix=fix, ...)
      }
      ResizeReadsMinus <- function(reads, width=1, fix="start", ...) {
         reads <- as(reads, "GRanges")
         #stopifnot(all(strand(reads) != "*"))
         resize(reads, width=width, fix=fix, ...)
      }
      summOv <- function(x) {
         featPlus <- rds@fragments[strand(rds@fragments)=="+"]
         featMinus <- rds@fragments[strand(rds@fragments)=="-"]
         plus <- summarizeOverlaps(features=featPlus, reads=x, 
                                    ignore.strand=!stranded, mc.cores=rds@cores
                                    preprocess.reads=ResizeReadsPlus)
         minus <- summarizeOverlaps(features=featMinus, reads=x, 
                                    ignore.strand=!stranded, mc.cores=rds@cores
                                    preprocess.reads=ResizeReadsMinus)
         return(rbind(plus, minus)) # Does this work as expected? XXX
      }
      if (length(rds@treatmentBams) == 1 && length(rds@controlBams) == 1) {
         # We obtain counts for both conditions on gene fragments and them sum 
         # them to obtain PRE/POST counts.
         treatmentSE <- summOv(rds@treatmentBams[[1]])
         controlSE <- summOv(rds@controlBams[[1]])    
         rds <- generateRoarsSingleBAM(rds, treatmentSE, controlSE)
      } else {
         # Still to be implemented.
      }
   }       
)

# BOf, why?
setMethod("generateRoarsSingleBAM", signature(rds="RoarDatasetMultipleAPA", 
                                              "RangedSummarizedExperiment",
                                              "RangedSummarizedExperiment"),
   function(rds, treatmentSE, controlSE)
   {
      # treatmentSE and controlSE are MoreArgs?
      rds@roars <- mapply(createRoarSingleBAM, rds@fragments, rds@prePostDef,
                       treatmentSE, controlSE)
      # set treatmentBams controlBams step paired (cores)
      return(rds)
   }
)

setMethod("computeRoars", signature(rds="RoarDatasetMultipleAPA"),
         function(rds) 
         {
            rds@roars <- lapply(rds@roars, computeRoars)
         }
)

setMethod("computePvals", signature(rds="RoarDatasetMultipleAPA"),
         function(rds)
         {
            rds@roars <- lapply(rds@roars, computePvals)
         }
)

setMethod("computePairedPvals", 
            signature(rds="RoarDatasetMultipleAPA",
            treatmentSamples="numeric", controlSamples="numeric"),
         function(rds, treatmentSamples, controlSamples) 
         {
            rds@roars <- lapply(rds@roars, computePairedPvals, 
                                treatmentSamples, controlSamples)
         }
)

setMethod("fpkmResults", signature(rds="RoarDatasetMultipleAPA"),
         function(rds) 
         {
            fpkmRes <- sapply(rds@roars, fpkmResults)
            # XXX TODO knit together results.
         }
)

