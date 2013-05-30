
# Do we need a generic?
setGeneric("countPrePost", def=function(rds) { standardGeneric("countPrePost") }, valueClass="RoarDataset") 

setGeneric("cores", def=function(rds) { standardGeneric("cores") }, valueClass="numeric")
setGeneric("cores<-", def=function(rds, value) { standardGeneric("cores<-") })
#setGeneric("cores", def=function(rds, value) { standardGeneric("cores<-") }) # FOAD
# HOWTO USE INTERFACE?