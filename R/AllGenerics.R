
setGeneric("lines", signature=c("x","y")) # To define new method of a class
setMethod("lines",c("ROC","missing"), function(x,y,...){
   lines(x@mspec, x@sens,...)
})
