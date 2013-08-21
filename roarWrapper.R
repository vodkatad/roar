#!/usr/bin/env Rscript
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
# options(error=traceback) 

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
   stop(getopt(arguments, command=get_Rscript_filename(), usage=T))
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

if (length(rightBams) > 1 || length(leftBams) > 1) {
   warning("Statistical analysis for replicates right now is not fully implemented.")
}

if (!all(sapply(c(rightBams, leftBams, opt$gtf), checkReadable))) {
   stop("One of the given files does not exist or is not readable")  
}

# Get counts
roar <- RoarDataset(rightBams, leftBams, opt$gtf)
roar <- countPrePost(roar, FALSE)

# Get m/M and Roar
roar <- computeRoars(roar)

# Fisher test
roar <- computePvals(roar)

# Obtain results with FPKM info for filtering and p-value correction
results <- filteringInfoResults(roar)
write.table(results, sep="\t", quote=F)
         
# XXX TODO FILTER AND CORRECTION.

if (!is.null(opt$debug)) {
   save.image(file=opt$debug)
}