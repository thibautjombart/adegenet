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
    res <- ktab.data.frame(df=as.data.frame(from), blocks=from@loc.n.all, rownames=indNames(from),
                           colnames=unlist(alleles(from)), tabnames=locNames(from))
    return(res)
})




setAs("genpop", "ktab", function(from, to) {
    checkType(from)
    res <- ktab.data.frame(df=as.data.frame(from), blocks=from@loc.n.all, rownames=popNames(from),
                           colnames=unlist(alleles(from)), tabnames=locNames(from))
    return(res)
})




##############
# S3 versions
##############

#' @export
as.data.frame.genind <- function(x,...){
    return(as.data.frame(tab(x, ...)))
}



#' @export
as.data.frame.genpop <- function(x,...){
    return(as.data.frame(tab(x, ...)))
}



#' @method as.matrix genind
#' @export
as.matrix.genind <- function(x,...){
    return(tab(x, ...))
}



#' @method as.matrix genpop
#' @export
as.matrix.genpop <- function(x,...){
    return(tab(x, ...))
}



#' @method as.genpop genind
#' @export
as.genpop.genind <- function(x,...){
    return(as(x,"genpop"))
}




as.ktab.genind <- function(x,...){
    return(as(x,"ktab"))
}




as.ktab.genpop <- function(x,...){
    return(as(x,"ktab"))
}

