#!/usr/bin/env Rscript
# Script to perform Roar analysis. 
# Requires a gtf with _PRE and _POST gene_ids and bam files from the two 
# conditions to be compared.

roarAnalysis <- function(gtf, treatmentBams, controlBams)
{
   ## This will need to be run foreach gene and foreach PRE/POST choice for them.
   # Get counts
   rds <- RoarDataset(treatmentBams, controlBams, gtf)
   rds <- countPrePost(rds, FALSE)
   
   # Get m/M and Roar
   rds <- computeRoars(rds)
   
   # Fisher test
   rds <- computePvals(rds)
   
   results <- fpkmResults(rds)
   apa_used <- unique(mcols(gtf)$apa_used)
   return(cbind(apa_used, results))
   
   #return(data.frame(matrix(c(rnorm(2)), nrow=1, ncol=2), row.names=apa_used))
   
   # filteredResults <- standardFilter(rds, fpkmCutoff=1)
   # write.table(filteredResults, sep="\t", quote=FALSE)
   
   # pvals <- pvalueFilter(rds, fpkmCutoff = 1, pvalCutoff = 0.05)
   # write.table(pvals, sep="\t", quote=FALSE)
}

# A function that puts together all roar results (or choose among them)
# calling them with every possible PRE/POST choice
callRoar <- function(gtf, treatmentBams, controlBams)
{
   return(lapply(gtf, roarAnalysis, treatmentBams, controlBams))
}

# A function that for a gene and its overlapping apas producesa
# a GRangesList object with all the choices of PRE/POST as GRanges. 
getAllPrePost <- function(geneGr, apaGr)
{
   apaGr <- sort(apaGr)
   # XXX check che ordini senza badare agli strand! (e idem start/end)
   strand <- unique(as.character(strand(geneGr)))
   chr <- unique(as.character(seqnames(geneGr)))
   gene_id <- unique(as.character(mcols(geneGr)$gene))
   if (length(chr) != 1) {
      stop("A gene is on different chrs")
   }
   if (length(gene_id) != 1) {
      stop("A gene has different names")
   }
   if (strand == '+') {
      apas <- head(apaGr, n=length(apaGr)-1)
      last <- tail(apaGr, n=1)
   } else if (strand == '-') {
      apas <- tail(apaGr, n=length(apaGr)-1)
      last <- head(apaGr, n=1)
   } else {
      stop("A gene has no strand info or both + and -")
   }
   res <- lapply(apas, FUN=definePrePost, last, geneGr, strand, chr, gene_id)
   return(GRangesList(res))
}

definePrePost <- function(firstApa, secondApa, geneGr, strand, chr, gene_id)
{
   introns <- gaps(geneGr)
   # We remove the first intron (from the beginning of the chr)
   introns <- tail(introns, n=length(introns)-1)
   mcols(geneGr) <- NULL
   mcols(geneGr)$type <- 'e'
   mcols(introns)$type <- 'i'   
   whole <- sort(c(geneGr, introns))
   hitsPre <- findOverlaps(whole, firstApa)
   #hitsPost <- findOverlaps(whole, secondApa)
   # We do not map the secondAPA: we cannot be sure that it has to map (it could
   # fall outside our gene in theory).
   foundPre <- NA
   if (!(countSubjectHits(hitsPre) == 1)) {
      # Then also the PRE is outside our gene. We use the last exon.
      #stop("Error: a given apa does not overlap its gene", 
      #     print(geneGr), print(firstApa))
      foundPre <- length(whole)
      if (strand == "-") {
         foundPre <- 1
      }
      if (whole[foundPre]$type != "e") {
         stop("Two apas outside the gene and I am unable to find the last exon", print(geneGr))
      }
   } else {
      foundPre <- queryHits(hitsPre)
   }
   # strand == "+"
   # XXX TODO reason about +-1 at begins and ends XXX TODO
   endPost <- start(secondApa)
   startPost <- end(firstApa) +1
   # We add 1 to avoid overlapping PRE/POST - APA are considered before the cut.
   # In the previous gtf we skipped the "cut" base altogether.
   startPre <- NA
   if (mcols(whole[foundPre])$type =="e") {
      startPre <- start(whole[foundPre])
      if (strand == "-") {
         startPre <- end(whole[foundPre])
      }
   } else { # if (mcols(foundPre)$type =="i") # true by costruction
      # We should not fall over boundaries if all preconditions are met.
      previousExon <- foundPre-1
      startPre <- start(whole[previousExon])
      if (strand == "-") {
         previousExon <- foundPre+1
         startPre <- end(whole[previousExon])
      }
   }
   endPre <- start(firstApa)
   if (strand == "-") {
      # We subtract 1 to avoid overlapping PRE/POST 
      # - APA are considered before the cut.
      # In the previous gtf we skipped the "cut" base altogether.
      endPost <- start(firstApa)-1
      startPost <- end(secondApa)
      sw <- startPre
      startPre <- endPre
      endPre <- sw
   }
   firstApa_name <- unlist(strsplit(mcols(firstApa)$apa, '_', fixed=TRUE))[1]
   names <- unlist(strsplit(mcols(secondApa)$apa, '_', fixed=TRUE))
   secondApa_name <- names[1]
   apa_names <- paste(firstApa_name, secondApa_name, sep="-")
   gene <- names[2]
   apa_names <- paste(gene, apa_names, sep="_")
   post <- GRanges(seqnames=chr, strand=strand, 
                   ranges=IRanges(start=startPost, end=endPost), 
                   gene_id=paste(gene_id, "POST", sep="_"),
                   apa_used=apa_names)
   pre <- GRanges(seqnames=chr, strand=strand, 
                   ranges=IRanges(start=startPre, end=endPre), 
                   gene_id=paste(gene_id, "PRE", sep="_"),
                   apa_used=apa_names)
   return(c(pre, post))
}

checkReadable <- function(filename) 
{
   res <- file.access(names=filename, mode=4) == 0
   if (!res) {
      warning(paste(filename, "is not readable", sep=" "))
   }
   res
}

printResults <- function(listRes)
{
   write.table(do.call(rbind, listRes), sep="\t", quote=FALSE, col.names=FALSE)
}

arguments <- matrix(c(
   'help', 'h', 0, "logical",
   'debug', 'd', 1, "character",
   'gtf' , 'a', 1, "character",
   'treatment'  , 't', 1, "character",
   'control'  , 'c', 1, "character"
), ncol=4, byrow=T)

library(getopt)
opt <- getopt(arguments)

if (!is.null(opt$help)) {
   stop(getopt(arguments, command=get_Rscript_filename(), usage=TRUE))
}

if (is.null(opt$gtf)) {
   stop("Missing gtf [-a filename] annotation option\n")
}

if (is.null(opt$treatment) | is.null(opt$control)) {
   stop("Missing treatment or control",
         "[-t, -c followed by comma separated bam files] param")
}

library(roar)
library(rtracklayer)
library(GenomicAlignments)
treatmentBams <- as.vector(unlist(strsplit(opt$treatment, ",")))
controlBams <- as.vector(unlist(strsplit(opt$control, ",")))
if (!all(sapply(c(treatmentBams, controlBams, opt$gtf), checkReadable))) {
   stop("One of the given files does not exist or is not readable")  
}
treatmentBamsGenomicAlignments <- lapply(treatmentBams, readGAlignments)
controlBamsGenomicAlignments <- lapply(controlBams, readGAlignments)

gtfGRanges<- import(opt$gtf, asRangedData=FALSE)
apas_melted <- gtfGRanges[mcols(gtfGRanges)$type=="apa"]
genes_melted <- gtfGRanges[mcols(gtfGRanges)$type=="gene"]
mcols(apas_melted)$gene <- as.numeric(sapply(
                                 strsplit(mcols(apas_melted)$apa, '_', 
                                          fixed=TRUE),
                                 function(x) { x[length(x)]}))
genes_ids <- sort(unique(mcols(genes_melted)$gene))
genes_ids_apas <- sort(unique(mcols(apas_melted)$gene))
if (!all(genes_ids==genes_ids_apas)) {
   stop("Lists of GAlignments could not be empty")
}
genes <- do.call(GRangesList, sapply(genes_ids, 
                                    function(x) {
                                       genes_melted[mcols(genes_melted)$gene==x]
                                    }))
apas <- do.call(GRangesList, sapply(genes_ids, 
                                    function(x) {
                                       apas_melted[mcols(apas_melted)$gene==x]
                                    }))
names(apas) <- genes_ids
names(genes) <- genes_ids

# We want a list of GRangesList: foreach gene a GRangesList 
# object with pre/post. 
# Are they sorted by name? They should be after the sapply. 
# XXX Check if that's true.
all_pre_post <- mapply(getAllPrePost, genes, apas)
# Foreach list we have a function that puts together all roar results 
# (or choose among them)
# calling them with every possible PRE/POST choice.
res <- lapply(all_pre_post, callRoar, 
               treatmentBamsGenomicAlignments, controlBamsGenomicAlignments)
# A list with n. elements == n. genes, elements are list of results for every
# PRE/POST choice for the gene.
garbage <- lapply(res, printResults)

if (!is.null(opt$debug)) {
   save.image(file=opt$debug)
}