#############
# S4 methods
#############
setAs("genind", "data.frame", function(from, to) {
    return(as.data.frame(from@tab))
})



setAs("genpop", "data.frame", function(from, to) {
    return(as.data.frame(from@tab))
})



setAs("genind", "matrix", function(from, to) {
    return(from@tab)
})



setAs("genpop", "matrix", function(from, to) {
    return(from@tab)
})



setAs("genind", "genpop", function(from, to) {
    if(!is.genind(from)) stop("object is not a valid genind")

    x <- genind2genpop(from, quiet=TRUE)
    warning("You had better use genind2genpop to specify treatment of NAs")

    return(x@tab)
})




setOldClass("ktab")
setAs("genind", "ktab", function(from, to) {
    checkType(from)
    res <- ktab.data.frame(df=as.data.frame(from), blocks=from@loc.nall, rownames=from@ind.names,
                           colnames=unlist(from@all.names), tabnames=from@loc.names)
    return(res)
})




setAs("genpop", "ktab", function(from, to) {
    checkType(from)
    res <- ktab.data.frame(df=as.data.frame(from), blocks=from@loc.nall, rownames=from@pop.names,
                           colnames=unlist(from@all.names), tabnames=from@loc.names)
    return(res)
})




##############
# S3 versions
##############
as.data.frame.genind <- function(x,...){
    return(as(x,"data.frame"))
}



as.data.frame.genpop <- function(x,...){
    return(as(x,"data.frame"))
}



as.matrix.genind <- function(x,...){
    return(as(x,"matrix"))
}



as.matrix.genpop <- function(x,...){
    return(as(x,"matrix"))
}



as.genpop.genind <- function(x,...){
    return(as(x,"genpop"))
}




as.ktab.genind <- function(x,...){
    return(as(x,"ktab"))
}




as.ktab.genpop <- function(x,...){
    return(as(x,"ktab"))
}

