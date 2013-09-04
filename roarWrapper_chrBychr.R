#!/rogue_bis/data/R/bin/Rscript
# First script to perform Roar analysis. 
# Requires a gtf with _PRE and _POST gene_ids and bam files from the two 
# conditions to be compared.
# Will become a wrapper for our Bioconductor library?

checkReadable <- function(filename) {
   res <- file.access(names=filename, mode=4) == 0
   if (!res) {
      warning(paste(filename, "is not readable", sep=" "))
   }
   res
}

# To understand at least something about occurred errors.
# In this way it goes further on after the first error. TODO avoid this!
options(error=traceback) 

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

#if (length(rightBams) > 1 || length(leftBams) > 1) {
#   warning("Statistical analysis for replicates right now is not fully implemented.")
#}

if (!all(sapply(c(rightBams, leftBams, opt$gtf), checkReadable))) {
   stop("One of the given files does not exist or is not readable")  
}

gtfGRanges<- import(opt$gtf, asRangedData=FALSE)
ordered <- order(elementMetadata(gtfGRanges)$gene_id)
gtfGRanges <- gtfGRanges[ordered]
chrs <- seqlevels(gtfGRanges)

workOnChr <- function(chr) {
   reduced <- keepSeqlevels(gtfGRanges, chr) # just a try
 
   loadBam <- function(bam) {
      tmp <- tempfile()
      garbage <- sortBam(bam, tmp, byQname=FALSE, maxMemory=512)
      garbage <- indexBam(paste(tmp,"bam",sep="."))
      
      param <- ScanBamParam(what=c("rname", "strand", "pos", "qwidth"), which=reduced)
      res <- readGappedAlignments(file=paste(tmp,"bam",sep="."), param = param)
      unlink(paste(tmp,"bam",sep="."))
      return(res)
   } 

   rightBamsGenomicAlignments <- lapply(rightBams, loadBam)
   leftBamsGenomicAlignments <- lapply(leftBams, loadBam)
   
   rds <- RoarDataset(rightBamsGenomicAlignments, leftBamsGenomicAlignments, reduced)
   # Get counts
   
   rds <- countPrePost(rds, FALSE)
   
   # Get m/M and Roar
   rds <- computeRoars(rds)
   
   # Fisher test
   rds <- computePvals(rds)
   
   # results <- filteringInfoResults(rds)
   # write.table(results, sep="\t", quote=FALSE)
            
   # filteredResults <- standardFilter(rds, fpkmCutoff=1)
   # write.table(filteredResults, sep="\t", quote=FALSE)
   
   res <- standardFilter(rds, fpkmCutoff = 1)
}

allRes <- lapply(chrs, workOnChr)
meltedRes <- do.call("rbind", allRes)
write.table(meltedRes, sep="\t", quote=FALSE)

if (!is.null(opt$debug)) {
   save.image(file=opt$debug)
}
