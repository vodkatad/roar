
setGeneric("countPrePost", def=function(rds, stranded=FALSE) { standardGeneric("countPrePost") }, valueClass="RoarDataset")  
setGeneric("computeRoars", def=function(rds) { standardGeneric("computeRoars") }, valueClass="RoarDataset")
setGeneric("computePvals", def=function(rds) { standardGeneric("computePvals") }, valueClass="RoarDataset")
setGeneric("computePairedPvals", def=function(rds, treatmentSamples, controlSamples) { standardGeneric("computePairedPvals") }, valueClass="RoarDataset")
setGeneric("totalResults", def=function(rds) { standardGeneric("totalResults") }, valueClass="data.frame")
setGeneric("fpkmResults", def=function(rds) { standardGeneric("fpkmResults") }, valueClass="data.frame")
setGeneric("countResults", def=function(rds) { standardGeneric("countResults") }, valueClass="data.frame")
setGeneric("standardFilter", def=function(rds, fpkmCutoff) { standardGeneric("standardFilter") }, valueClass="data.frame")
setGeneric("pvalueFilter", def=function(rds, fpkmCutoff, pvalCutoff) { standardGeneric("pvalueFilter") }, valueClass="data.frame")
setGeneric("pvalueCorrectFilter", def=function(rds, fpkmCutoff, pvalCutoff, method) { standardGeneric("pvalueCorrectFilter") }, valueClass="data.frame")

setGeneric("generateRoarsSingleBam", def=function(rds, treatmentSE, controlSE) {
   standardGeneric("generateRoarsSingleBam") }, valueClass="RoarDatasetMultipleAPA")
setGeneric("mcountPrePost", def=function(rds, stranded=FALSE) { standardGeneric("mcountPrePost") }, valueClass="RoarDatasetMultipleAPA")  
setGeneric("mcomputeRoars", def=function(rds) { standardGeneric("mcomputeRoars") }, valueClass="RoarDatasetMultipleAPA")
setGeneric("mcomputePvals", def=function(rds) { standardGeneric("mcomputePvals") }, valueClass="RoarDatasetMultipleAPA")
setGeneric("mcomputePairedPvals", def=function(rds, treatmentSamples, controlSamples) { standardGeneric("mcomputePairedPvals") }, valueClass="RoarDatasetMultipleAPA")
setGeneric("mfpkmResults", def=function(rds) { standardGeneric("mfpkmResults") }, valueClass="data.frame")


setGeneric("cores", def=function(rds) { standardGeneric("cores") }, valueClass="numeric")
#setGeneric("cores<-", def=function(rds, value) { standardGeneric("cores<-") })
