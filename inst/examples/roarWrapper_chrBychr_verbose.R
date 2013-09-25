#!/rogue_bis/data/R/bin/Rscript
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

start.time <- Sys.time()
orderedRightBams <- lapply(rightBams, orderBam)
orderedLeftBams <- lapply(leftBams, orderBam)
end.time <- Sys.time()
time.taken <- end.time - start.time
print("# order bam took")
print(time.taken)


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
   start.time <- Sys.time()
   rightBamsGenomicAlignments <- lapply(orderedRightBams, loadBam)
   leftBamsGenomicAlignments <- lapply(orderedLeftBams, loadBam)
   end.time <- Sys.time()
   time.taken <- end.time - start.time
   print("# loading bam took")
   print(time.taken)
   
   start.time <- Sys.time()
   rds <- RoarDataset(rightBamsGenomicAlignments, leftBamsGenomicAlignments, reduced)
   # Get counts
   rds <- countPrePost(rds, FALSE)
   end.time <- Sys.time()
   time.taken <- end.time - start.time
   print("# creating rds and countPrePost took")
   print(time.taken)
   # Get m/M and Roar
   start.time <- Sys.time()
   rds <- computeRoars(rds)
   end.time <- Sys.time()
   time.taken <- end.time - start.time
   print("# computeRoars took")
   print(time.taken)
   # Fisher test
   start.time <- Sys.time()
   rds <- computePvals(rds)
   end.time <- Sys.time()
   time.taken <- end.time - start.time
   print("# computePvals took")
   print(time.taken)
   start.time <- Sys.time()
   res <- countResults(rds)
   end.time <- Sys.time()
   time.taken <- end.time - start.time
   print("# countResults took")
   print(time.taken)
   return(res)
}

start.time <- Sys.time()
allRes <- lapply(chrs, workOnChr)
meltedRes <- do.call("rbind", allRes)
preElems <- grep("_PRE$", elementMetadata(gtfGRanges)$gene_id)
pre <- gtfGRanges[preElems,]
preLen <- end(pre) - start(pre) + 1
names <- sub("^\\s+","",sub("_PRE", "",elementMetadata(pre)$gene_id))
meltedRes <- meltedRes[match(names, rownames(meltedRes)),]
sumPreRight <- sum(meltedRes[,"rightValue"])
sumPreLeft <- sum(meltedRes[,"leftValue"])
meltedRes$rightFpkm <- (meltedRes[,"rightValue"]*1000000000)/(preLen*sumPreRight)
meltedRes$leftFpkm <- (meltedRes[,"leftValue"]*1000000000)/(preLen*sumPreLeft)
write.table(meltedRes, sep="\t", quote=FALSE)
end.time <- Sys.time()
time.taken <- end.time - start.time
print("# put together took")
print(time.taken)
# unlink all bam now I don't care

if (!is.null(opt$debug)) {
   save.image(file=opt$debug)
}
