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
options(error=traceback) 

library(getopt)

arguments <- matrix(c(
   'help', 'h', 0, "logical",
   'gtf' , 'a', 1, "character",
   'right'  , 'r', 1, "character",
   'left'  , 'l', 1, "character"
), ncol=4, byrow=T)

opt <- getopt(arguments)

if (!is.null(opt$help)) {
   stop(getopt(arguments, command=get_Rscript_filename(), usage=T))
}

if (is.null(opt$gtf)) {
   stop("Missing gtf [-g filename] annotation option\n")
}

if (is.null(opt$right) | is.null(opt$left)) {
   stop("Missing right or left [-r, -l followed by comma separated bam files] param")
}

rightBams = as.vector(unlist(strsplit(opt$right, ",")))
leftBams = as.vector(unlist(strsplit(opt$left, ",")))

if (length(rightBams) > 1 || length(leftBams) > 1) {
   warning(c("Statistical analysis for replicates right now is not fully implemented.\n", 
            "Simple sum of counts will be used to determine the p-value for the Fisher test."))
}

if (!all(sapply(c(rightBams, leftBams, opt$gtf), checkReadable))) {
   stop("One of the given files does not exist or is not readable")  
}

# Get counts
roar <- RoarDataset(rightBams, leftBams, opt$gtf)
roar <- CountPrePost(roar)

# Get m/M and Roar

# Fisher test

# Filter results based on PRE counts and bonferroni p-value correction
