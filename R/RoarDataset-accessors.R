# Methods for the RoarDataset class.
# R/AllGenerics.R contains the signature of the not overriding methods.

CreateRoarDataset <- function(rightBams, leftBams, gtf) {
   # TODO
   gtfGenomicAnnot <- qualcosa(gtf)
   new("RoarDataset", rightBams=rightBams, leftBams=leftBams, prePostCoords=gtfGenomicAnnot)
}

# Could have used setMethod("initialize", "xx",) but in this way should have had a gtf filename slot.

setMethod("show","ROC", # To override "standard methods"
          function(object){
             cat("ROC curve: ")
             print(object@call)
          })
