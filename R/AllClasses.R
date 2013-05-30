# This class will old a GenomicAnnotation (pre-post GTF) and count results in a SummarizedExperiment
setClass( "RoarDataset",
          contains = "SummarizedExperiment", # Superclasses
          representation = representation( # Slots
             rightBams = "list", #list of GappedAlignments (were)
             leftBams = "list",
             prePostCoords = "GRanges",
             cores = "numeric"
          )
#validity=function(roc) length(roc@sens)==length(roc@mspec)
#          && length(roc@sens)==length(roc@test)
# TODO check validity of kind of alignments, existence of PRE/POST ids
)
# prototype -> default values