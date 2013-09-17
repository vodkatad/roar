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
   'right'  , 'r', 1, "character",
   'left'  , 'l', 1, "character"
), ncol=4, byrow=T)

library(getopt)
opt <- getopt(arguments)

if (!is.null(opt$help)) {
   stop(getopt(arguments, command=get_Rscript_filename(), usage=TRUE))
}

if (is.null(opt$gtf)) {
   stop("Missing gtf [-a filename] annotation option\n")
}

if (is.null(opt$right) | is.null(opt$left)) {
   stop("Missing right or left [-r, -l followed by comma separated bam files] param")
}

library(roar)
rightBams = as.vector(unlist(strsplit(opt$right, ",")))
leftBams = as.vector(unlist(strsplit(opt$left, ",")))

if (!all(sapply(c(rightBams, leftBams, opt$gtf), checkReadable))) {
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

orderedRightBams <- lapply(rightBams, orderBam)
orderedLeftBams <- lapply(leftBams, orderBam)

workOnChr <- function(chr) {
   write(paste("Working on", chr), stderr())
   reduced <- keepSeqlevels(gtfGRanges, chr) 
   coords <- c(start(reduced), end(reduced))  # To keep strandness in consideration! Is this needed?
   begin <- min(coords)
   end <- max(coords)
   spanChr <- GRanges(seqnames=chr,ranges=IRanges(start=begin,width=end-begin+1))
   loadBam <- function(bam) {
      param <- ScanBamParam(which=spanChr)
      res <- readGappedAlignments(file=bam, param = param)
      return(res)
   } 

   rightBamsGenomicAlignments <- lapply(orderedRightBams, loadBam)
   leftBamsGenomicAlignments <- lapply(orderedLeftBams, loadBam)
   
   rds <- RoarDataset(rightBamsGenomicAlignments, leftBamsGenomicAlignments, reduced)
   
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
write.table(meltedRes, sep="\t", quote=FALSE)
# XXX TODO ADD knit together rpkm values

# unlink all bam now I don't care

if (!is.null(opt$debug)) {
   save.image(file=opt$debug)
}
