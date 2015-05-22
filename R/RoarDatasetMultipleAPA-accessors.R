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

setMethod("countPrePost", signature(rds="RoarDataset"),
         function(rds, stranded=FALSE) {
            allFragments <- mapply(getApaGenesFractions,
                                    rds@geneCoords, 
                                    rds@apaCoords)
            # A GRangesList with GRanges for all fragments defining pre/post
            # in a gene. 
         }
          
)

