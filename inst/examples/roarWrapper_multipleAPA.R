#!/home/data/work/R-devel/bin/Rscript
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

library(roar)
library(rtracklayer)
library(GenomicAlignments)
treatmentBams <- as.vector(unlist(strsplit(opt$treatment, ",")))
controlBams <- as.vector(unlist(strsplit(opt$control, ",")))

if (!all(sapply(c(treatmentBams, controlBams, opt$gtf), checkReadable))) {
   stop("One of the given files does not exist or is not readable")  
}

# Get counts
ptm <- proc.time()
rds <- RoarDatasetMultipleAPAFromFiles(treatmentBams, controlBams, opt$gtf)
gap <- proc.time() - ptm
cat("create obj", "\n", file=stderr())
cat(gap,"\n", file=stderr())
ptm <- proc.time()

rds <- countPrePost(rds, FALSE)
gap <- proc.time() - ptm
cat("count", "\n", file=stderr())
cat(gap,"\n", file=stderr())
ptm <- proc.time()

# Get m/M and Roar
rds <- computeRoars(rds)
gap <- proc.time() - ptm
cat("computeRoar", "\n", file=stderr())
cat(gap,"\n", file=stderr())
ptm <- proc.time()

# Fisher test
rds <- computePvals(rds)
gap <- proc.time() - ptm
cat("computePvals", "\n", file=stderr())
cat(gap,"\n", file=stderr())
ptm <- proc.time()


results <- totalResults(rds)
gap <- proc.time() - ptm
cat("results", "\n", file=stderr())
cat(gap,"\n", file=stderr())
ptm <- proc.time()

write.table(results, sep="\t", quote=FALSE)

garbage <- lapply(rds@roars, function(x) { export(x@prePostCoords, con = "prova2.gtf", append=TRUE)})

# filteredResults <- standardFilter(rds, fpkmCutoff=1)
# write.table(filteredResults, sep="\t", quote=FALSE)

# pvals <- pvalueFilter(rds, fpkmCutoff = 1, pvalCutoff = 0.05)
# write.table(pvals, sep="\t", quote=FALSE)

if (!is.null(opt$debug)) {
   save.image(file=opt$debug)
}
