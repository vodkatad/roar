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

Rprof(tf <- "rprof.log", memory.profiling=TRUE)

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
   start.time <- Sys.time()
   reduced <- keepSeqlevels(gtfGRanges, chr) 
   coords <- c(start(reduced), end(reduced))  # Keep strandness in consideration! Is this needed?
   begin <- min(coords)
   end <- max(coords)
   spanChr <- GRanges(seqnames=chr,ranges=IRanges(start=begin,width=end-begin+1))
   loadBam <- function(bam) {
      #param <- ScanBamParam(what=c("rname", "strand", "pos", "qwidth"), which=reduced)
      param <- ScanBamParam(which=spanChr)
      res <- readGappedAlignments(file=bam, param = param)
      return(res)
   } 

   rightBamsGenomicAlignments <- lapply(orderedRightBams, loadBam)
   leftBamsGenomicAlignments <- lapply(orderedLeftBams, loadBam)
   end.time <- Sys.time()
   time.taken <- end.time - start.time
   write(paste("Loading took", time.taken), stderr())
   
   start.time <- Sys.time()
   rds <- RoarDataset(rightBamsGenomicAlignments, leftBamsGenomicAlignments, reduced)
   end.time <- Sys.time()
   time.taken <- end.time - start.time
   write(paste("Building object took", time.taken), stderr())
   # Get counts
   
   start.time <- Sys.time()
   rds <- countPrePost(rds, FALSE)
   end.time <- Sys.time()
   time.taken <- end.time - start.time
   write(paste("Counting took", time.taken),  stderr())
   # Get m/M and Roar
   start.time <- Sys.time()
   rds <- computeRoars(rds)
   end.time <- Sys.time()
   time.taken <- end.time - start.time
   write(paste("Roars took", time.taken),  stderr())
   
   # Fisher test
   start.time <- Sys.time()
   rds <- computePvals(rds)
   end.time <- Sys.time()
   time.taken <- end.time - start.time
   write(paste("Pvalues took", time.taken), stderr())
   size <- object.size(rds)
   write("Size:",  stderr())
   write(size, stderr())
   # results <- filteringInfoResults(rds)
   # write.table(results, sep="\t", quote=FALSE)
            
   # filteredResults <- standardFilter(rds, fpkmCutoff=1)
   # write.table(filteredResults, sep="\t", quote=FALSE)
   start.time <- Sys.time()
   res <- filteringInfoResults(rds)
   end.time <- Sys.time()
   time.taken <- end.time - start.time
   write(paste("results took", time.taken), stderr())
   return(res)
}

allRes <- lapply(chrs, workOnChr)
meltedRes <- do.call("rbind", allRes)
write.table(meltedRes, sep="\t", quote=FALSE)
# unlink all bam now I don't care
Rprof(NULL)

if (!is.null(opt$debug)) {
   save.image(file=opt$debug)
}
