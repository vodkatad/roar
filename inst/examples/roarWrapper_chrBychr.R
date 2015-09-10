#!/usr/bin/env Rscript
# Script to perform stepwise (chr by chr) Roar analysis. 
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

library(roar)
library(rtracklayer)
library(Rsamtools)
library(GenomicAlignments)

treatmentBams <- as.vector(unlist(strsplit(opt$treatment, ",")))
controlBams <- as.vector(unlist(strsplit(opt$control, ",")))

if (!all(sapply(c(treatmentBams, controlBams, opt$gtf), checkReadable))) {
   stop("One of the given files does not exist or is not readable")  
}

gtfGRanges<- import(opt$gtf)
chrs <- seqlevels(gtfGRanges)

orderBam <- function(bam) {
      tmp <- tempfile()
      ordered <- sortBam(bam, tmp, byQname=FALSE, maxMemory=512)
      garbage <- indexBam(ordered)
      return(ordered)
}

orderedTreatmentBams <- lapply(treatmentBams, orderBam)
orderedControlBams <- lapply(controlBams, orderBam)

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

   treatmentBamsGenomicAlignments <- lapply(orderedTreatmentBams, loadBam)
   controlBamsGenomicAlignments <- lapply(orderedControlBams, loadBam)
   
   rds <- RoarDataset(treatmentBamsGenomicAlignments, controlBamsGenomicAlignments, reduced)
   
   # Get counts
   rds <- countPrePost(rds, FALSE)
   # Get m/M and Roar
   rds <- computeRoars(rds)
   # Fisher test
   rds <- computePvals(rds)
   res <- countResults(rds)
   return(res)
}

allRes <- lapply(chrs, workOnChr)
meltedRes <- do.call("rbind", allRes)
preElems <- grep("_PRE$", mcols(gtfGRanges)$gene_id)
pre <- gtfGRanges[preElems,]
preLen <- end(pre) - start(pre) + 1
names <- sub("^\\s+","",sub("_PRE", "",mcols(pre)$gene_id))
meltedRes <- meltedRes[match(names, rownames(meltedRes)),]
sumPreTreatment <- sum(meltedRes[,"treatmentValue"])
sumPreControl <- sum(meltedRes[,"controlValue"])
meltedRes$treatmentFpkm <- (meltedRes[,"treatmentValue"]*1000000000)/(preLen*sumPreTreatment)
meltedRes$controlFpkm <- (meltedRes[,"controlValue"]*1000000000)/(preLen*sumPreControl)
write.table(meltedRes, sep="\t", quote=FALSE)

unlink(orderedTreatmentBams)
unlink(orderedControlBams)
treatmentBai <- sub(".bam", ".bai", orderedTreatmentBams)
controlBai <- sub(".bam", ".bai", orderedControlBams)
unlink(treatmentBai)
unlink(controlBai)

if (!is.null(opt$debug)) {
   save.image(file=opt$debug)
}
