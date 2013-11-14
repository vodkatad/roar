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
rightBams <- as.vector(unlist(strsplit(opt$right, ",")))
leftBams <- as.vector(unlist(strsplit(opt$left, ",")))

if (!all(sapply(c(rightBams, leftBams, opt$gtf), checkReadable))) {
   stop("One of the given files does not exist or is not readable")  
}

# Get counts
rds <- RoarDatasetFromFiles(rightBams, leftBams, opt$gtf)
rds <- countPrePost(rds, FALSE)

# Get m/M and Roar
rds <- computeRoars(rds)

# Fisher test
rds <- computePvals(rds)

# results <- fpkmResults(rds)
# write.table(results, sep="\t", quote=FALSE)
         
# filteredResults <- standardFilter(rds, fpkmCutoff=1)
# write.table(filteredResults, sep="\t", quote=FALSE)

pvals <- pvalueFilter(rds, fpkmCutoff = 1, pvalCutoff = 0.05)
write.table(pvals, sep="\t", quote=FALSE)

if (!is.null(opt$debug)) {
   save.image(file=opt$debug)
}