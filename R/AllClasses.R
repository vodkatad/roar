# This class will hold a GenomicAnnotation (pre-post GTF) and count results in a RangedSummarizedExperiment
setClass( "RoarDataset",
          contains = "RangedSummarizedExperiment", # Superclasses
          representation = representation( # Slots
             treatmentBams = "list",  #list of GAlignments 
             controlBams = "list",
             prePostCoords = "GRanges",
             postCoords = "GRanges",
             countsTreatment = "RangedSummarizedExperiment",
             countsControl = "RangedSummarizedExperiment",
             pVals = "RangedSummarizedExperiment",
             step = "numeric",
             paired = "logical",
             cores = "numeric"
          )
#validity=function(roc) length(roc@sens)==length(roc@mspec)
#          && length(roc@sens)==length(roc@test)
# check validity of kind of alignments, existence of PRE/POST ids -> moved in other methods.
)
# prototype -> default values
