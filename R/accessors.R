#######
# nLoc
#######
setGeneric("nLoc", function(x,...){
    standardGeneric("nLoc")
})



# setMethod("nLoc","genind", function(x,...){
#   if (x@type == "PA"){
#     return(ncol(x@tab))
#   } else {
#     return(length(levels(x@loc.fac)))
#   }
# })
# 
# 
# 
# setMethod("nLoc","genpop", function(x,...){
#     if (x@type == "PA"){
#       return(ncol(x@tab))
#     } else {
#       return(length(levels(x@loc.fac)))
#     }
# 
# })

setMethod("nLoc","gen", function(x,...){
    if (x@type == "PA"){
      return(ncol(x@tab))
    } else {
      return(length(levels(x@loc.fac)))
    }
})

setMethod("nLoc", "genind", function(x, ...){
  callNextMethod()
})

setMethod("nLoc", "genpop", function(x, ...){
  callNextMethod()
})

#########
# locFac
#########
setGeneric("locFac", function(x, ...){
    standardGeneric("locFac")
})

setMethod("locFac","gen", function(x,...){
    return(x@loc.fac)
})

setMethod("locFac","genind", function(x,...){
  callNextMethod()
})

setMethod("locFac","genpop", function(x,...){
    callNextMethod()
})


#######
# nAll
#######
setGeneric("nAll", function(x,...){
    standardGeneric("nAll")
})


setMethod("nAll","gen", function(x,...){
  if (x@type == "PA"){
    return(ncol(x@tab))
  } else {
    return(x@loc.n.all)
  }
})

setMethod("nAll","genind", function(x,...){
  callNextMethod()
})

setMethod("nAll","genpop", function(x,...){
    callNextMethod()
})


#######
# nPop (no gen method)
#######
setGeneric("nPop", function(x,...){
    standardGeneric("nPop")
})

setMethod("nPop","genind", function(x,...){
    return(length(levels(x@pop)))
})

setMethod("nPop","genpop", function(x,...){
     return(nrow(x@tab))
})


#######
# nInd (no gen method)
#######
setGeneric("nInd", function(x,...){
    standardGeneric("nInd")
})



setMethod("nInd","genind", function(x,...){
    return(nrow(x@tab))
})





######
# pop (no gen method)
######
setGeneric("pop", function(x) {
  standardGeneric("pop")
})



setGeneric("pop<-",
           function(x, value) {
               standardGeneric("pop<-")
           })



setMethod("pop","genind", function(x){
    return(x@pop)
})



setReplaceMethod("pop", "genind", function(x, value) {
    if(is.null(value)){
        x@pop <- NULL
        return(x)
    }

    if(length(value) != nrow(x$tab)) stop("wrong length for population factor")

    ## coerce to factor (put levels in their order of appearance)
    newPop <- as.character(value)
    newPop <- factor(newPop, levels=unique(newPop))

    ## construct output
    x$pop <- newPop

    return(x)
})





###########
# locNames
###########
setGeneric("locNames", function(x,...){
    standardGeneric("locNames")
})

setGeneric("locNames<-", function(x, value) {
    standardGeneric("locNames<-")
})


setMethod("locNames","gen", function(x, withAlleles=FALSE, ...){
    if (withAlleles){
        res <- colnames(x@tab)
    } else {
        res <- unique(sub("[.][^.]*$","",colnames(x@tab)))
    }
    return(res)
})


setReplaceMethod("locNames","gen",function(x,value) {
    ## check input
    value <- as.character(value)
    if(length(value) != nLoc(x)) stop("Vector length does no match number of loci")

    ## make changes in the object
    names(x@all.names) <- value
    levels(x@loc.fac) <- value
    names(x@loc.n.all) <- value
    newColNames <- paste(rep(value, x@loc.n.all), unlist(x@all.names), sep=".")
    colnames(x@tab) <- newColNames

    ## return
    return(x)
})


setMethod("locNames","genpop", function(x, withAlleles=FALSE, ...){
    callNextMethod()
})
setReplaceMethod("locNames", "genpop", function(x, value) {
    callNextMethod()
})


setMethod("locNames","genind", function(x, withAlleles=FALSE, ...){
  callNextMethod()
})
setReplaceMethod("locNames", "genind", function(x, value) {
  callNextMethod()
})


###########
# indNames (no gen method)
###########
setGeneric("indNames", function(x,...){
    standardGeneric("indNames")
})

setGeneric("indNames<-", function(x, value){
    standardGeneric("indNames<-")
})

setMethod("indNames","genind", function(x, ...){
    return(rownames(x@tab))
})


setReplaceMethod("indNames","genind",function(x,value) {
    if(length(value) != nrow(x@tab)) stop("Vector length does not match number of individuals")
    rownames(x@tab) <- as.character(value)
    return(x)
})

###########
# popNames (no gen method)
###########
setGeneric("popNames", function(x,...){
  standardGeneric("popNames")
})

setGeneric("popNames<-", function(x, value){
  standardGeneric("popNames<-")
})

setMethod("popNames","genind", function(x, ...){
  return(levels(pop(x)))
})


setReplaceMethod("popNames","genind",function(x, value) {
  value <- as.character(value)
  if(length(value) != length(levels(pop(x)))){
    stop("Vector length does not match number of populations")
  }
  levels(pop(x)) <- value
  return(x)
})

setMethod("popNames","genpop", function(x, ...){
  return(rownames(x@tab))
})


setReplaceMethod("popNames","genpop",function(x, value) {
  value <- as.character(value)
  if (length(value) != nrow(x@tab)){
    stop("Vector length does not match number of populations")
  }
  rownames(x@tab) <- value
  return(x)
})


##########
# alleles
##########
setGeneric("alleles", function(x,...){
    standardGeneric("alleles")
})

setGeneric("alleles<-", function(x, value){
    standardGeneric("alleles<-")
})

setMethod("alleles","gen", function(x, ...){
    return(x@all.names)
})

setReplaceMethod("alleles","gen", function(x, value){
    if(!is.list(value)) stop("replacement value must be a list")
    if(length(value)!=nLoc(x)) stop("replacement list must be of length nLoc(x)")
    if(any(sapply(value, length) != x$loc.n.all)) stop("number of replacement alleles do not match that of the object")
    x@all.names <- value
    names(x@all.names) <- locNames(x)
    return(x)
})

setMethod("alleles","genind", function(x, ...){
  callNextMethod()
})
setReplaceMethod("alleles","genind", function(x, value){
  callNextMethod()
})

setMethod("alleles","genpop", function(x, ...){
    callNextMethod()
})
setReplaceMethod("alleles","genpop", function(x, value){
    callNextMethod()
})



##########
## ploidy (no gen method)
##########
setGeneric("ploidy", function(x,...){
    standardGeneric("ploidy")
})

setGeneric("ploidy<-", function(x, value){
    standardGeneric("ploidy<-")
})

setMethod("ploidy","genind", function(x,...){
    return(x@ploidy)
})


setReplaceMethod("ploidy","genind",function(x,value) {
    value <- as.integer(value)
    if(any(value)<1) stop("Negative or null values provided")
    if(any(is.na(value))) stop("NA values provided")
    if(length(value)!=nInd(x)) value <- rep(value, length=nInd(x))
    slot(x,"ploidy",check=TRUE) <- value
    return(x)
})


setMethod("ploidy","genpop", function(x,...){
    return(x@ploidy)
})


setReplaceMethod("ploidy","genpop",function(x,value) {
    value <- as.integer(value)
    if(any(value)<1) stop("Negative or null values provided")
    if(any(is.na(value))) stop("NA values provided")
    if(length(value)>1) warning("Several ploidy numbers provided; using only the first integer")
    slot(x,"ploidy",check=TRUE) <- value[1]
    return(x)
})






##########
## other
#########
setGeneric("other", function(x,...){
    standardGeneric("other")
})

setGeneric("other<-", function(x, value){
    standardGeneric("other<-")
})

setMethod("other","gen", function(x,...){
    if(length(x@other)==0) return(NULL)
    return(x@other)
})


setReplaceMethod("other","gen",function(x,value) {
    if( !is.null(value) && (!is.list(value) | is.data.frame(value)) ) {
        value <- list(value)
    }
    slot(x,"other",check=TRUE) <- value
    return(x)
})


setMethod("other","genind", function(x,...){
  callNextMethod()
})
setReplaceMethod("other","genind",function(x,value) {
  callNextMethod()
})

setMethod("other","genpop", function(x,...){
    callNextMethod()
})
setReplaceMethod("other","genpop",function(x,value) {
    callNextMethod()
})


