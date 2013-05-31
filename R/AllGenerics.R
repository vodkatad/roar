
setGeneric("countPrePost", def=function(rds) { standardGeneric("countPrePost") }, valueClass="RoarDataset") 
setGeneric("computeRoars", def=function(rds) { standardGeneric("computeRoars") }, valueClass="RoarDataset")
setGeneric("computePvals", def=function(rds) { standardGeneric("computePvals") }, valueClass="RoarDataset")
setGeneric("formatResults", def=function(rds) { standardGeneric("formatResults") }, valueClass="data.frame")

setGeneric("cores", def=function(rds) { standardGeneric("cores") }, valueClass="numeric")
setGeneric("cores<-", def=function(rds, value) { standardGeneric("cores<-") })
#setGeneric("cores", def=function(rds, value) { standardGeneric("cores<-") }) # FOAD
# HOWTO USE INTERFACE?