#!/usr/bin/env Rscript
# Script to perform Roar analysis. 
# Requires a gtf with _PRE and _POST gene_ids and bam files from the two 
# conditions to be compared.

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
options(error = function() traceback(2))
library(roar)
library(rtracklayer)
library(GenomicAlignments)
treatmentBams <- as.vector(unlist(strsplit(opt$treatment, ",")))
controlBams <- as.vector(unlist(strsplit(opt$control, ",")))
stranded <- FALSE

if (!all(sapply(c(treatmentBams, controlBams, opt$gtf), checkReadable))) {
   stop("One of the given files does not exist or is not readable")  
}


gtfGRanges<- import(opt$gtf, asRangedData=FALSE)
chrs <- seqlevels(gtfGRanges)

orderBam <- function(bam) {
   tmp <- tempfile()
   ordered <- sortBam(bam, tmp, byQname=FALSE, maxMemory=512)
   garbage <- indexBam(ordered)
   return(ordered)
}

orderedRightBams <- lapply(treatmentBams, orderBam)
orderedLeftBams <- lapply(controlBams, orderBam)

workOnChr <- function(chr) {
   write(paste("Working on", chr), stderr())
   reduced <- keepSeqlevels(gtfGRanges, chr) 
   coords <- c(start(reduced), end(reduced))  # To keep strandness in consideration! Is this needed?
   begin <- min(coords)
   end <- max(coords)
   spanChr <- GRanges(seqnames=chr,ranges=IRanges(start=begin,width=end-begin+1))
   loadBam <- function(bam) {
      param <- ScanBamParam(which=spanChr)
      res <- readGAlignments(file=bam, param = param)
      return(res)
   } 
   
   treatmentBamsGenomicAlignments <- lapply(orderedRightBams, loadBam)
   controlBamsGenomicAlignments <- lapply(orderedLeftBams, loadBam)
   
   rds <- RoarDatasetMultipleAPA(treatmentBamsGenomicAlignments, controlBamsGenomicAlignments, reduced)
   
   # Get counts
   rds <- countPrePost(rds, stranded)
   # Get m/M and Roar
   rds <- computeRoars(rds)
   # Fisher test
   rds <- computePvals(rds)
   res <- countResults(rds)
   return(res)
}

allRes <- lapply(chrs, workOnChr)
meltedRes <- do.call("rbind", allRes)

sumTreatment <- sum(meltedRes[,"counts_treatment"])
sumControl <- sum(meltedRes[,"counts_control"])
meltedRes$treatmentFpkm <- (meltedRes[,"counts_treatment"]*1000000000)/(meltedRes[,"length"]*sumTreatment)
meltedRes$controlFpkm <- (meltedRes[,"counts_control"]*1000000000)/(meltedRes[,"length"]*sumControl)
meltedRes$treatmentValue <- NULL
meltedRes$controlValue <- NULL
write.table(meltedRes, sep="\t", quote=FALSE)

unlink(orderedRightBams)
unlink(orderedLeftBams)
treatmentBai <- sub(".bam", ".bai", orderedRightBams)
controlBai <- sub(".bam", ".bai", orderedLeftBams)
unlink(treatmentBai)
unlink(controlBai)

if (!is.null(opt$debug)) {
   save.image(file=opt$debug)
}
