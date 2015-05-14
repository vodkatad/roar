#!/usr/bin/env Rscript
# Script to perform Roar analysis. 
# Requires a gtf with _PRE and _POST gene_ids and bam files from the two 
# conditions to be compared.

roarAnalysis <- function(gtf, treatmentBams, controlBams)
{
   ## This will need to be run foreach gene and foreach PRE/POST choice for them.
   # Get counts
   rds <- RoarDatasetFromFiles(treatmentBams, controlBams, gtf)
   rds <- countPrePost(rds, FALSE)
   
   # Get m/M and Roar
   rds <- computeRoars(rds)
   
   # Fisher test
   rds <- computePvals(rds)
   
   results <- fpkmResults(rds)
   
   return(fpkmResults)
   # filteredResults <- standardFilter(rds, fpkmCutoff=1)
   # write.table(filteredResults, sep="\t", quote=FALSE)
   
   # pvals <- pvalueFilter(rds, fpkmCutoff = 1, pvalCutoff = 0.05)
   # write.table(pvals, sep="\t", quote=FALSE)
}

# A function that puts together all roar results (or choose among them)
# calling them with every possible PRE/POST choice
callRoar <- function(treatmentBams, controlBams, gtf)
{
   
}

# A function that for a gene and its overlapping apas producesa
# a GRangesList object with all the choices of PRE/POST as GRanges. 
getAllPrePost <- function(geneGr, apaGr)
{
   introns <- gaps(geneGr)
   mcols(geneGr) <- NULL
   mcols(geneGr)$type <- 'e'
   mcols(introns)$type <- 'i'   
   whole <- sort(c(geneGr, introns))
   hits <- findOverlaps(whole, apaGr)
   if (!all(countSubjectHits(hits) == 1)) {
      stop("Error: a given apa does not overlap its gene")
   }
   foundov <- length(unique(whole[queryHits(hits)]))
   mcols(whole)$overlap <- rep(FALSE, length(whole))   
   mcols(whole[unique(queryHits(hits))])$overlap <- rep(TRUE, foundov)
   if (whole[1]$type != 'i' || whole[1]$overlap) {
      stop("Error in the given gene structure or overlapping apas")
   }
}

checkReadable <- function(filename) {
   res <- file.access(names=filename, mode=4) == 0
   if (!res) {
      warning(paste(filename, "is not readable", sep=" "))
   }
   res
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
   stop("Missing treatment or control [-t, -c followed by comma separated bam files] param")
}

library(roar)
library(rtracklayer)
treatmentBams <- as.vector(unlist(strsplit(opt$treatment, ",")))
controlBams <- as.vector(unlist(strsplit(opt$control, ",")))

if (!all(sapply(c(treatmentBams, controlBams, opt$gtf), checkReadable))) {
   stop("One of the given files does not exist or is not readable")  
}


gtfGRanges<- import(gtf, asRangedData=FALSE)
apas_melted <- gtfGRanges[mcols(gtfGRanges)$type=="apa"]
genes_melted <- gtfGRanges[mcols(gtfGRanges)$type=="gene"]
mcols(apas_melted)$gene <- sapply(strsplit(mcols(apas_melted)$apa, '_', 
                                           fixed=TRUE),
                                  function(x) { x[length(x)]})
genes_ids <- sort(unique(mcols(genes_melted)$gene))
genes_ids_apas <- sort(unique(mcols(apas_melted)$gene))
names(apas) <- genes_ids
names(genes) <- genes_ids
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

# We want a list of GRangesList: foreach gene a GRangesList 
# object with pre/post. 
all_pre_post <- mapply(getAllPrePost, genes, apas)
# Foreach list we have a function that puts together all roar results 
# (or choose among them)
# calling them with every possible PRE/POST choice.
res <- qualkapply(callRoar, all_pre_post)



if (!is.null(opt$debug)) {
   save.image(file=opt$debug)
}