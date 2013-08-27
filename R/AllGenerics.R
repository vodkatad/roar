
setGeneric("countPrePost", def=function(rds, stranded) { standardGeneric("countPrePost") }, valueClass="RoarDataset") 
setGeneric("computeRoars", def=function(rds) { standardGeneric("computeRoars") }, valueClass="RoarDataset")
setGeneric("computePvals", def=function(rds) { standardGeneric("computePvals") }, valueClass="RoarDataset")
setGeneric("totalResults", def=function(rds) { standardGeneric("totalResults") }, valueClass="data.frame")
setGeneric("filteringInfoResults", def=function(rds) { standardGeneric("filteringInfoResults") }, valueClass="data.frame")
setGeneric("standardFilter", def=function(rds, fpkmCutoff) { standardGeneric("standardFilter") }, valueClass="data.frame")
setGeneric("pvalueFilter", def=function(rds, fpkmCutoff, pvalCutoff) { standardGeneric("pvalueFilter") }, valueClass="data.frame")


setGeneric("cores", def=function(rds) { standardGeneric("cores") }, valueClass="numeric")
#setGeneric("cores<-", def=function(rds, value) { standardGeneric("cores<-") })
