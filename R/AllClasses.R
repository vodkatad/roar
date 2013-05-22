# This class will old a GenomicAnnotation (pre-post GTF) and count results in a SummarizedExperiment
setClass( "RoarDataset",
          contains = "SummarizedExperiment" # Superclasses
          representation = representation( # Slots
             rightBams = "character",
             leftBams = "character",
             prePostCoords = "GenomicAnnotation"
          ),
#validity=function(roc) length(roc@sens)==length(roc@mspec)
#          && length(roc@sens)==length(roc@test)
          
)
# prototype -> default values