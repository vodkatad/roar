setGeneric("countPrePost", def=function(rds, stranded=FALSE) { standardGeneric("countPrePost") }, valueClass=c("RoarDataset", "RoarDatasetMultipleAPA"))
setGeneric("computeRoars", def=function(rds, qwidthTreatment = NA, qwidthControl = NA) { standardGeneric("computeRoars") }, valueClass=c("RoarDataset", "RoarDatasetMultipleAPA"))
setGeneric("computePvals", def=function(rds) { standardGeneric("computePvals") }, valueClass=c("RoarDataset", "RoarDatasetMultipleAPA"))
setGeneric("computePairedPvals", def=function(rds, treatmentSamples, controlSamples) { standardGeneric("computePairedPvals") }, valueClass=c("RoarDataset", "RoarDatasetMultipleAPA"))
setGeneric("totalResults", def=function(rds) { standardGeneric("totalResults") }, valueClass="data.frame")
setGeneric("fpkmResults", def=function(rds) { standardGeneric("fpkmResults") }, valueClass="data.frame")
setGeneric("countResults", def=function(rds) { standardGeneric("countResults") }, valueClass="data.frame")
setGeneric("standardFilter", def=function(rds, fpkmCutoff) { standardGeneric("standardFilter") }, valueClass="data.frame")
setGeneric("pvalueFilter", def=function(rds, fpkmCutoff, pvalCutoff) { standardGeneric("pvalueFilter") }, valueClass="data.frame")
setGeneric("pvalueCorrectFilter", def=function(rds, fpkmCutoff, pvalCutoff, method) { standardGeneric("pvalueCorrectFilter") }, valueClass="data.frame")

setGeneric("generateRoarsSingleBam", def=function(rds, treatmentSE, controlSE) {
   standardGeneric("generateRoarsSingleBam") }, valueClass="RoarDatasetMultipleAPA")
setGeneric("generateRoarsMultipleBam", def=function(rds, treatmentSE, controlSE) {
   standardGeneric("generateRoarsMultipleBam") }, valueClass="RoarDatasetMultipleAPA")

setGeneric("cores", def=function(rds) { standardGeneric("cores") }, valueClass="numeric")
#setGeneric("cores<-", def=function(rds, value) { standardGeneric("cores<-") })
