setGeneric("quality", function(object) standardGeneric("quality"))
setGeneric("quality<-", function(object, value) standardGeneric("quality<-"))

setMethod("quality", signature(object="qPCRset"), definition = 
    function (object) {x <- assayDataElement(object, "quality"); rownames(x) <- featureNames(object); x}
)

setReplaceMethod("quality", signature(object="qPCRset"), definition = 
    function (object, value) assayDataElementReplace(object, "quality", value)
)
