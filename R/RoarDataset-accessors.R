# Methods for the RoarDataset class.
# R/AllGenerics.R contains the signature of the not overriding methods.

RoarDataset <- function(rightBams, leftBams, gtf) {
   genome = "hg19" # FIXME
   # The format will be assumed using the file extension. Will work everytime?
   # Do we force hg19?
   gtfGRanges<- import(gtf, genome=genome, asRangedData=F)
   rightBamsGenomicAlignments <- lapply(rightBams, readGappedAlignments)
   leftBamsGenomicAlignments <- lapply(leftBams, readGappedAlignments)
   new("RoarDataset", rightBams=rightBamsGenomicAlignments, leftBams=leftBamsGenomicAlignments, prePostCoords=gtfGRanges)
}

# Could have used setMethod("initialize", "xx",) but in this way should have had a gtf filename slot.

setMethod("show","ROC", # To override "standard methods"
          function(object){
             cat("ROC curve: ")
             print(object@call)
          })
